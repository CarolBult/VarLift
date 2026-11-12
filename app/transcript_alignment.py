'''
Codon-aware alignment of a query transcript and a homologous transcript, for display.

The two proteins are aligned with the same aligner VarLift uses to compare residue
positions (app/protein_alignment.py), and each transcript's coding sequence is then
threaded through the protein alignment codon by codon, so the bases line up under
the residues they encode. Very long proteins (e.g. titin) are aligned through a
window around the variant only.
'''
import json
import yaml
import logging
from functools import lru_cache
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq

from app.protein_alignment import (_sequences, _species, _full_alignment, _global_aligner,
                                   homologous_position, MAX_FULL_ALIGNMENT_CELLS, SequenceUnavailable)


with open('app/config.yaml') as f:
    config = yaml.safe_load(f)

# Residues either side of the variant when only a window is aligned (very long proteins)
LONG_PROTEIN_WINDOW = 150


@lru_cache(maxsize=2)
def _coding_sequences(species):
    '''Transcript ID -> coding sequence, or None when the parquet hasn't been built.'''
    path = config['paths'][f'{species}_cds_pqt']
    if not Path(path).exists():
        logging.warning(f'{path} not found; run app/data_refresh.py to build it. '
                        'Alignments are shown without codons.')
        return None
    df = pd.read_parquet(path)
    return dict(zip(df['transcript_id'], df['sequence']))


@lru_cache(maxsize=2)
def _domains(species):
    '''Transcript ID -> list of Pfam domains {id, name, start, end}, or {} when the parquet hasn't been built.'''
    path = config['paths'][f'{species}_domain_pqt']
    if not Path(path).exists():
        logging.warning(f'{path} not found; run app/data_refresh.py to build it. Alignments are shown without domains.')
        return {}
    df = pd.read_parquet(path)
    names = _pfam_names()
    df['name'] = [n if isinstance(n, str) and n else names.get(i, i) for i, n in zip(df['domain_id'], df['name'])]
    out = {}
    for tx, i, n, a, b in zip(df['transcript_id'], df['domain_id'], df['name'], df['start'], df['end']):
        out.setdefault(tx, []).append({'id': i, 'name': n, 'start': int(a), 'end': int(b)})
    return out


def _pfam_names():
    '''Pfam accession -> name, from the names VEP lookups have cached.'''
    try:
        with open(config['paths']['pfam_names']) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def _domain_at(domains, residue):
    '''ID of the domain containing a residue (the first, where domains overlap), or None.'''
    if residue is None:
        return None
    return next((d['id'] for d in domains if d['start'] <= residue <= d['end']), None)


def _codons(cds, protein):
    '''
    The codon for each residue of protein, or None when the coding sequence doesn't
    translate to it. Transcripts with an incomplete 5' end start out of frame, so the
    frame that best reproduces the protein is used.
    '''
    if not cds:
        return None
    best = None
    for frame in range(3):
        seq = cds[frame:]
        seq = seq[:len(seq) - len(seq) % 3]
        translated = str(Seq(seq).translate())
        n = min(len(translated), len(protein))
        if n == 0:
            continue
        agree = sum(a == b for a, b in zip(translated[:n], protein[:n]) if b != 'X')
        if best is None or agree > best[0]:
            best = (agree, frame, n)
    agree, frame, n = best
    if agree < 0.95 * n:
        return None
    codons = [cds[frame + 3 * i:frame + 3 * i + 3] for i in range(len(protein))]
    return [c if len(c) == 3 else None for c in codons]


def _segment(seq, centre, half):
    '''(start, end) 0-based half-open around a 1-based centre residue.'''
    start = max(0, centre - 1 - half)
    return start, min(len(seq), centre + half)


def _alignment(seq_q, seq_t, query_residue, query_tx, query_species, target_tx, target_species):
    '''
    Gapped query and target strings and the 0-based offsets they start at. Long
    proteins are aligned only through a window around the query residue.
    '''
    if len(seq_q) * len(seq_t) <= MAX_FULL_ALIGNMENT_CELLS:
        aln = _full_alignment(seq_q, seq_t)
        return aln[0], aln[1], 0, 0, False

    if not query_residue:
        raise SequenceUnavailable('protein too long to align without a variant residue')
    mapped = homologous_position(query_tx, query_species, target_tx, target_species, query_residue) or query_residue
    q0, q1 = _segment(seq_q, query_residue, LONG_PROTEIN_WINDOW)
    t0, t1 = _segment(seq_t, mapped, LONG_PROTEIN_WINDOW)
    aln = _global_aligner.align(seq_q[q0:q1], seq_t[t0:t1])[0]
    return aln[0], aln[1], q0, t0, True


def _gap_runs(columns):
    '''Runs of gap columns: which protein lacks the residues, how many, and where.'''
    runs, start = [], None
    for i, col in enumerate(columns + [None]):
        in_gap = col is not None and (col['q'] == '-' or col['t'] == '-')
        if in_gap and start is None:
            start = i
        elif not in_gap and start is not None:
            run = columns[start:i]
            missing_in = 'target' if run[0]['t'] == '-' else 'query'
            present = 'qp' if missing_in == 'target' else 'tp'
            before = next((columns[k]['qp'] for k in range(start - 1, -1, -1) if columns[k]['qp']), 0)
            runs.append({'missing_in': missing_in, 'length': i - start,
                         'from': run[0][present], 'to': run[-1][present], 'after_query_residue': before})
            start = None
    return runs


def transcript_alignment(query_tx, query_species, target_tx, target_species,
                         query_residue=None, target_residue=None, flank=None):
    '''
    Codon-aware alignment of two transcripts' proteins.

    Parameters:
        query_tx, target_tx: str. Ensembl transcript IDs.
        query_species, target_species: str. 'human' or 'mouse'.
        query_residue, target_residue: int or None. 1-based variant residues, marked in the result.
        flank: int or None. Residues either side of the query residue to return; None for the whole alignment.

    Returns:
        dict. Transcript details, summary statistics, gap runs and the alignment columns, each
        {'q', 't': residue or '-', 'qp', 'tp': residue number or None, 'qc', 'tc': codon or None}.

    Raises:
        SequenceUnavailable: when either protein sequence is missing.
    '''
    q_species, t_species = _species(query_species), _species(target_species)
    seqs_q, seqs_t = _sequences(q_species), _sequences(t_species)
    seq_q = seqs_q.get(query_tx) if seqs_q else None
    seq_t = seqs_t.get(target_tx) if seqs_t else None
    if not seq_q or not seq_t:
        raise SequenceUnavailable(f'{query_tx} or {target_tx}')

    cds_q, cds_t = _coding_sequences(q_species), _coding_sequences(t_species)
    codons_q = _codons(cds_q.get(query_tx) if cds_q else None, seq_q)
    codons_t = _codons(cds_t.get(target_tx) if cds_t else None, seq_t)

    gapped_q, gapped_t, iq, it, windowed = _alignment(seq_q, seq_t, query_residue,
                                                      query_tx, q_species, target_tx, t_species)
    domains_q = _domains(q_species).get(query_tx, [])
    domains_t = _domains(t_species).get(target_tx, [])

    columns = []
    for a, b in zip(gapped_q, gapped_t):
        col = {'q': a, 't': b, 'qp': None, 'tp': None, 'qc': None, 'tc': None}
        if a != '-':
            col['qp'] = iq + 1
            col['qc'] = codons_q[iq] if codons_q else None
            iq += 1
        if b != '-':
            col['tp'] = it + 1
            col['tc'] = codons_t[it] if codons_t else None
            it += 1
        col['qd'] = _domain_at(domains_q, col['qp'])
        col['td'] = _domain_at(domains_t, col['tp'])
        columns.append(col)

    paired = [c for c in columns if c['q'] != '-' and c['t'] != '-']
    identical = sum(c['q'] == c['t'] for c in paired)
    with_codons = [c for c in paired if c['qc'] and c['tc']]
    same_bases = sum(x == y for c in with_codons for x, y in zip(c['qc'], c['tc']))

    aligned_to = next((c['tp'] for c in columns if c['qp'] == query_residue), None) if query_residue else None

    shown = columns
    if flank is not None and query_residue:
        centre = next((i for i, c in enumerate(columns) if c['qp'] == query_residue), None)
        if centre is not None:
            shown = columns[max(0, centre - flank):centre + flank + 1]

    return {
        'query': {'transcript': query_tx, 'species': q_species, 'length': len(seq_q),
                  'residue': query_residue, 'codons': codons_q is not None, 'domains': domains_q},
        'target': {'transcript': target_tx, 'species': t_species, 'length': len(seq_t),
                   'residue': target_residue, 'codons': codons_t is not None, 'domains': domains_t},
        'aligned_to_query_residue': aligned_to,
        'windowed': windowed,
        'stats': {
            'columns': len(columns),
            'aligned_residues': len(paired),
            'identical_residues': identical,
            'protein_identity': round(100 * identical / len(paired), 1) if paired else None,
            'cds_identity': round(100 * same_bases / (3 * len(with_codons)), 1) if with_codons else None,
            'gap_columns': len(columns) - len(paired),
        },
        'gaps': _gap_runs(columns),
        'columns': shown,
        'complete': shown is columns,
    }


def _blocks(gapped_q, gapped_t, q_offset, t_offset):
    '''
    Runs of aligned residues as (query start, query end, target start, target end), 1-based
    inclusive, and the target residues with no query counterpart as (query position they
    follow, target start, target end, count).
    '''
    blocks, inserts = [], []
    iq, it = q_offset, t_offset
    block = insert = None
    for a, b in zip(gapped_q, gapped_t):
        if a != '-':
            iq += 1
        if b != '-':
            it += 1
        if a != '-' and b != '-':
            if insert:
                inserts.append(insert); insert = None
            if block and block[1] == iq - 1 and block[3] == it - 1:
                block[1], block[3] = iq, it
            else:
                if block:
                    blocks.append(block)
                block = [iq, iq, it, it]
        elif b != '-':
            if insert:
                insert['to'] = it; insert['count'] += 1
            else:
                insert = {'after': iq, 'from': it, 'to': it, 'count': 1}
    if block:
        blocks.append(block)
    if insert:
        inserts.append(insert)
    return blocks, inserts


def _query_position_of(blocks, target_residue):
    '''The query residue aligned to a target residue, or None when it has no counterpart.'''
    for q0, _, t0, t1 in blocks:
        if t0 <= target_residue <= t1:
            return q0 + (target_residue - t0)
    return None


def _segments_in_query(domains, blocks):
    '''
    A target transcript's domains in query coordinates: the parts of each domain that align
    to the query protein, as {id, name, start, end} (query residues); parts with no query
    counterpart are left out.
    '''
    segments = []
    for d in domains:
        for q0, _, t0, t1 in blocks:
            lo, hi = max(d['start'], t0), min(d['end'], t1)
            if lo <= hi:
                segments.append({'id': d['id'], 'name': d['name'], 'start': q0 + lo - t0, 'end': q0 + hi - t0})
    return segments


def transcript_overview(query_tx, query_species, query_residue, targets, target_species):
    '''
    Every target transcript aligned to one query transcript, in query coordinates, for
    drawing the transcripts of an allele (or of a human variant) one above the other.

    Parameters:
        query_tx: str. Ensembl transcript ID of the query (submitted variant's) transcript.
        query_species, target_species: str. 'human' or 'mouse'.
        query_residue: int or None. 1-based residue of the submitted variant in the query protein.
        targets: list of (transcript ID, variant residue or None).

    Returns:
        dict. 'query' (transcript, length, residue) and 'targets': for each, its length, aligned
        blocks and insertions, the residue aligned to the query residue, and the query residue
        its own variant residue aligns to. A target whose protein is missing has 'error' set.
    '''
    q_species, t_species = _species(query_species), _species(target_species)
    seqs_q, seqs_t = _sequences(q_species), _sequences(t_species)
    seq_q = seqs_q.get(query_tx) if seqs_q else None
    if not seq_q:
        raise SequenceUnavailable(query_tx)

    out = []
    for tx, residue in targets:
        seq_t = seqs_t.get(tx) if seqs_t else None
        if not seq_t:
            out.append({'transcript': tx, 'error': 'No protein sequence available'})
            continue
        try:
            gapped_q, gapped_t, q0, t0, windowed = _alignment(seq_q, seq_t, query_residue,
                                                              query_tx, q_species, tx, t_species)
        except SequenceUnavailable:
            out.append({'transcript': tx, 'error': 'Too long to align without a variant residue'})
            continue
        blocks, inserts = _blocks(gapped_q, gapped_t, q0, t0)
        paired = sum(b[1] - b[0] + 1 for b in blocks)
        identical = sum(a == b for a, b in zip(gapped_q, gapped_t) if a != '-' and b != '-')
        aligned_to = next((t_start + (query_residue - q_start) for q_start, q_end, t_start, _ in blocks
                           if query_residue and q_start <= query_residue <= q_end), None)
        own_domains = _domains(t_species).get(tx, [])
        out.append({
            'transcript': tx,
            'length': len(seq_t),
            'domains': own_domains,
            'domain_segments': _segments_in_query(own_domains, blocks),
            'residue': residue,
            'aligned_to_query_residue': aligned_to,
            'residue_at_query': _query_position_of(blocks, residue) if residue else None,
            'blocks': blocks,
            'insertions': inserts,
            'protein_identity': round(100 * identical / paired, 1) if paired else None,
            'windowed': windowed,
        })
    return {'query': {'transcript': query_tx, 'species': q_species, 'length': len(seq_q), 'residue': query_residue,
                      'domains': _domains(q_species).get(query_tx, [])},
            'target_species': t_species, 'targets': out}
