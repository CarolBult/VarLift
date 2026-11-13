'''
Structural equivalency between a query variant and its homologous variants.

Each query variant is scored by the class of its molecular consequence, and each
pair gets an Equivalency Tier (1 = most equivalent):

    Coding
        1  Identical change: same amino-acid change at the aligned residue, where
           the alignment around the residue is reliable (scoring.min_local_identity)
        2  Same residue: aligned residue, different change; an identical change
           where the alignment is uncertain; or both truncating
           (stop_gained/frameshift) within scoring.truncation_window residues
        3  Same domain: same consequence, in Pfam domains of the same clan
        4  Same consequence, elsewhere in the protein
        5  Different consequence
    Splice/Intronic
        1  Same splice site and offset (equivalent intron, same HGVS offset)
        2  Same splice site: equivalent intron, both splice-site/region variants on
           the same side of the intron
        3  Same intron
        4  Different intron
        5  Different consequence

Positions are compared through a protein alignment, domains by Pfam clan, and introns
by where they split the coding sequence (also through the alignment), so differences in
residue, exon and intron numbering between orthologs and isoforms don't count against a
pair. The True/False match columns are kept as the evidence behind each tier.
'''
import re
import yaml
import logging
import numpy as np
import pandas as pd
from pathlib import Path

from app.protein_alignment import homologous_position, local_identity, SequenceUnavailable



with open('app/config.yaml') as f:
    config = yaml.safe_load(f)

TRUNCATION_WINDOW = config.get('scoring', {}).get('truncation_window', 10)
ALIGNMENT_FLANK = config.get('scoring', {}).get('alignment_flank', 10)
MIN_LOCAL_IDENTITY = config.get('scoring', {}).get('min_local_identity', 0.3)


CODING = 'Coding'
SPLICE = 'Splice/Intronic'

CODING_TERMS = {'missense_variant', 'synonymous_variant', 'stop_gained', 'stop_lost', 'start_lost',
                'stop_retained_variant', 'start_retained_variant', 'frameshift_variant',
                'inframe_insertion', 'inframe_deletion', 'protein_altering_variant',
                'coding_sequence_variant', 'incomplete_terminal_codon_variant'}

SPLICE_TERMS = {'splice_donor_variant', 'splice_acceptor_variant', 'splice_donor_5th_base_variant',
                'splice_donor_region_variant', 'splice_region_variant',
                'splice_polypyrimidine_tract_variant', 'intron_variant'}

# Criteria that only apply to one class
TRUNCATING_TERMS = {'stop_gained', 'frameshift_variant'}

# Splice-site and splice-region terms (an intron_variant alone is deeper in the intron)
SPLICE_SITE_TERMS = SPLICE_TERMS - {'intron_variant'}

# Criteria that only apply to one class
CLASS_CRITERIA = {CODING: ['AA Match', 'AA Position Match', 'Domain Match'],
                  SPLICE: ['Intron Match', 'Intron Offset Match']}

# Exon number isn't compared: exon numbering differs between orthologs and isoforms, and
# the aligned residue already places a coding variant
STRUCTURAL_COLS = ['Biotype Match', 'Consequence Match', 'AA Match', 'AA Position Match',
                   'Domain Match', 'Intron Match', 'Intron Offset Match']

TIER_BASIS = {
    CODING: {1: 'Identical change', 2: 'Same residue', 3: 'Same domain', 4: 'Same consequence',
             5: 'Different consequence'},
    SPLICE: {1: 'Same splice site and offset', 2: 'Same splice site', 3: 'Same intron',
             4: 'Different intron', 5: 'Different consequence'},
}


def _terms(consequence):
    if isinstance(consequence, list):
        return set(consequence)
    if isinstance(consequence, str):
        return set(consequence.split(','))
    return set()


def scoring_class(consequence):
    '''
    Scoring class of a VEP consequence: Coding when any term changes the coding
    sequence, Splice/Intronic for splice-site and intronic terms, otherwise None.
    '''
    terms = _terms(consequence)
    if terms & CODING_TERMS:
        return CODING
    if terms & SPLICE_TERMS:
        return SPLICE
    return None


def _exons_for(transcripts, species):
    path = config['paths'][f'{species}_exon_pqt']
    if not Path(path).exists():
        logging.warning(f'{path} not found; run app/data_refresh.py to build it. Intron offsets are left blank.')
        return pd.DataFrame(columns=['transcript_id', 'rank', 'start', 'end', 'strand'])
    return pd.read_parquet(path, filters=[('transcript_id', 'in', list(transcripts))])


def add_intron_offsets(vep_df, species):
    '''
    Add intron_offset: the signed distance from an intronic variant to the nearest
    exon boundary of its transcript, written as in HGVS (+5 is five bases into the
    intron after an exon, -2 two bases before the next exon). Exonic rows get None.

    Parameters:
        vep_df: DataFrame. Parsed VEP output with Submission, transcript_id and intron.
        species: str. Either 'homo_sapiens' or 'mus_musculus'.
    '''
    vep_df['intron_offset'] = None
    if vep_df.empty or 'intron' not in vep_df.columns:
        return vep_df

    intronic = vep_df['intron'].notna() & vep_df['transcript_id'].notna()
    if not intronic.any():
        return vep_df

    exons = _exons_for(vep_df.loc[intronic, 'transcript_id'].unique(),
                       'human' if species == 'homo_sapiens' else 'mouse')
    exons = {t: g.set_index('rank') for t, g in exons.groupby('transcript_id')}

    def _offset(row):
        tx = exons.get(row['transcript_id'])
        try:
            rank = int(str(row['intron']).split('/')[0])
            pos = int(row['Submission'].split('\t')[1])
            before, after = tx.loc[rank], tx.loc[rank + 1]
        except (AttributeError, KeyError, ValueError, TypeError, IndexError):
            return None

        if before['strand'] == '+':
            donor, acceptor = pos - before['end'], after['start'] - pos
        else:
            donor, acceptor = before['start'] - pos, pos - after['end']

        if donor <= 0 or acceptor <= 0:
            return None
        # HGVS places the variant relative to the nearer exon, the earlier one on a tie
        return f'+{donor}' if donor <= acceptor else f'-{acceptor}'

    vep_df.loc[intronic, 'intron_offset'] = vep_df[intronic].apply(_offset, axis=1)
    return vep_df


def _match(a, b):
    '''True/False when both have a value, None when neither does, False when only one does.'''
    a_null, b_null = pd.isna(a), pd.isna(b)
    if a_null and b_null:
        return None
    if a_null or b_null:
        return False
    return bool(a == b)


def _aa_change(df, sfx):
    ref, var = df[f'refAA{sfx}'], df[f'varAA{sfx}']
    return [f'{r}/{v}' if pd.notna(r) else None for r, v in zip(ref, var)]


def _residue_position(amino_acids):
    '''Residue number from a change like R145C or RE1812X (frameshift), or None.'''
    match = re.match(r'^[A-Z*-]+(\d+)', amino_acids) if isinstance(amino_acids, str) else None
    return int(match.group(1)) if match else None


def _aligned_positions(df, q, c, q_species, c_species):
    '''
    For each pair, the homolog residue aligned to the query's residue and whether it is
    the homolog's variant residue. Positions are compared by number only when the
    protein sequences aren't available.
    '''
    aligned, matches = [], []
    for q_tx, q_aa, c_tx, c_aa in zip(df[f'Transcript ID{q}'], df[f'Amino Acids{q}'],
                                      df[f'Transcript ID{c}'], df[f'Amino Acids{c}']):
        q_pos, c_pos = _residue_position(q_aa), _residue_position(c_aa)
        if q_pos is None or c_pos is None:
            aligned.append(None)
            matches.append(_match(q_pos, c_pos))
            continue
        try:
            mapped = homologous_position(q_tx, q_species, c_tx, c_species, q_pos)
        except SequenceUnavailable:
            mapped = q_pos
        aligned.append(mapped)
        matches.append(mapped == c_pos)
    return aligned, matches


def _pfam_by_transcript(transcripts, species):
    '''Transcript ID -> list of (start, end, clan or Pfam ID, Pfam ID, name) for its Pfam domains.'''
    path = config['paths'].get(f'{species}_domain_pqt')
    if not path or not Path(path).exists():
        return None
    df = pd.read_parquet(path, columns=['transcript_id', 'domain_id', 'name', 'clan', 'start', 'end'],
                         filters=[('transcript_id', 'in', list(transcripts))])
    # Families outside any clan are compared by their own ID
    key = df['clan'].where(df['clan'].notna() & (df['clan'] != ''), df['domain_id'])
    out = {}
    for tx, a, b, k, i, n in zip(df['transcript_id'], df['start'], df['end'], key, df['domain_id'], df['name']):
        out.setdefault(tx, []).append((int(a), int(b), k, i, n))
    return out


def _clans_at(domains, tx, residue):
    if residue is None:
        return set()
    return {d[2] for d in domains.get(tx, []) if d[0] <= residue <= d[1]}


def _domain_label(domains, tx, residue):
    '''"PF00069 Protein kinase domain" for the Pfam domain(s) at a residue, or ''.'''
    if not domains or residue is None:
        return ''
    found = {d[3]: d[4] for d in domains.get(tx, []) if d[0] <= residue <= d[1]}
    return ', '.join(f'{i} {n}' if isinstance(n, str) and n else i for i, n in sorted(found.items()))


def _domain_matches(df, q, c, q_species, c_species):
    '''
    Whether the two variant residues fall in Pfam domains of the same clan (any of the
    domains at each residue, not only the first VEP reports). Falls back to comparing
    VEP's Pfam IDs when the domain tables haven't been built.
    Also returns, per pair, the query residue's domain(s) for the tier label.
    '''
    q_dom = _pfam_by_transcript(df[f'Transcript ID{q}'].dropna().unique(), q_species)
    c_dom = _pfam_by_transcript(df[f'Transcript ID{c}'].dropna().unique(), c_species)
    if q_dom is None or c_dom is None:
        return [_match(a, b) for a, b in zip(df[f'Pfam Domain ID{q}'], df[f'Pfam Domain ID{c}'])], [''] * len(df)
    out, labels = [], []
    for q_tx, q_aa, c_tx, c_aa in zip(df[f'Transcript ID{q}'], df[f'Amino Acids{q}'],
                                      df[f'Transcript ID{c}'], df[f'Amino Acids{c}']):
        a = _clans_at(q_dom, q_tx, _residue_position(q_aa))
        b = _clans_at(c_dom, c_tx, _residue_position(c_aa))
        out.append(None if not a and not b else bool(a & b))
        labels.append(_domain_label(q_dom, q_tx, _residue_position(q_aa)))
    return out, labels


def _coding_junctions(transcripts, species):
    '''
    Transcript ID -> {intron number: last coding base before it}, for introns inside the
    coding sequence. None when the exon table has no CDS positions.
    '''
    path = config['paths'][f'{species}_exon_pqt']
    if not Path(path).exists():
        return None
    exons = pd.read_parquet(path, filters=[('transcript_id', 'in', list(transcripts))])
    if 'cds_end' not in exons.columns:
        return None
    out = {}
    for tx, g in exons.sort_values('rank').groupby('transcript_id'):
        ends = dict(zip(g['rank'], g['cds_end']))
        starts = dict(zip(g['rank'], g['cds_start']))
        out[tx] = {r: int(e) for r, e in ends.items()
                   if pd.notna(e) and pd.notna(starts.get(r + 1)) and starts[r + 1] == e + 1}
    return out


def _intron_number(intron):
    try:
        return int(str(intron).split('/')[0])
    except (TypeError, ValueError):
        return None


def _intron_matches(df, q, c, q_species, c_species):
    '''
    Whether two intronic variants are in equivalent introns: the query intron's position
    in the coding sequence (residue and codon phase) is carried through the protein
    alignment, and the homolog's intron must split its coding sequence at the aligned
    residue (within one residue) with the same phase. Introns outside the coding sequence
    are compared by number.
    '''
    q_j = _coding_junctions(df[f'Transcript ID{q}'].dropna().unique(), q_species.replace('homo_sapiens', 'human'))
    c_j = _coding_junctions(df[f'Transcript ID{c}'].dropna().unique(), c_species.replace('homo_sapiens', 'human'))
    out = []
    for q_tx, q_in, c_tx, c_in in zip(df[f'Transcript ID{q}'], df[f'Intron Rank{q}'],
                                      df[f'Transcript ID{c}'], df[f'Intron Rank{c}']):
        qn, cn = _intron_number(q_in), _intron_number(c_in)
        if qn is None or cn is None:
            out.append(_match(qn, cn))
            continue
        qj = (q_j or {}).get(q_tx, {}).get(qn)
        cj = (c_j or {}).get(c_tx, {}).get(cn)
        if qj is None or cj is None:
            out.append(qn == cn if qj is None and cj is None else False)
            continue
        try:
            mapped = homologous_position(q_tx, q_species, c_tx, c_species, (qj - 1) // 3 + 1)
        except SequenceUnavailable:
            out.append(qn == cn)
            continue
        out.append(mapped is not None and qj % 3 == cj % 3 and abs((cj - 1) // 3 + 1 - mapped) <= 1)
    return out


def _local_identities(df, q, c, q_species, c_species, position_matches):
    '''Local alignment identity around the query residue, for pairs whose residues align.'''
    out = []
    for q_tx, q_aa, c_tx, matched in zip(df[f'Transcript ID{q}'], df[f'Amino Acids{q}'],
                                         df[f'Transcript ID{c}'], position_matches):
        pos = _residue_position(q_aa)
        if not matched or pos is None:
            out.append(None)
            continue
        try:
            out.append(round(local_identity(q_tx, q_species, c_tx, c_species, pos, ALIGNMENT_FLANK), 2))
        except SequenceUnavailable:
            out.append(None)
    return out


def _main_terms(consequence):
    '''Coding, splice and intronic terms, without modifiers like NMD_transcript_variant.'''
    return _terms(consequence) & (CODING_TERMS | SPLICE_TERMS)


def _tier(cls, r):
    '''Equivalency tier of one pair (see the module docstring), or None for unscored classes.'''
    q_terms, c_terms = r['q_terms'], r['c_terms']
    same_consequence = bool(q_terms & c_terms)
    if cls == CODING:
        if r['AA Position Match'] and r['AA Match']:
            return 2 if r['uncertain'] else 1
        if r['AA Position Match']:
            return 2
        if (q_terms & TRUNCATING_TERMS and c_terms & TRUNCATING_TERMS and r['aligned'] is not None
                and r['c_pos'] is not None and abs(r['aligned'] - r['c_pos']) <= TRUNCATION_WINDOW):
            return 2
        if same_consequence and r['Domain Match']:
            return 3
        return 4 if same_consequence else 5
    if cls == SPLICE:
        if not c_terms & SPLICE_TERMS:
            return 5
        if not r['Intron Match']:
            return 4
        if r['Intron Offset Match']:
            return 1
        q_off, c_off = str(r['q_off'] or ''), str(r['c_off'] or '')
        if q_terms & SPLICE_SITE_TERMS and c_terms & SPLICE_SITE_TERMS and q_off[:1] == c_off[:1]:
            return 2
        return 3
    return None


def structural_matches(df, q, c, q_species, c_species):
    '''
    Compare each query/homolog pair in a merged DataFrame.

    Parameters:
        df: DataFrame. One row per pair, query columns suffixed q and homolog columns suffixed c.
        q: str. Suffix of the query variant's columns.
        c: str. Suffix of the homologous variant's columns.
        q_species, c_species: str. 'human' or 'mouse', for the query and homolog.

    Returns:
        DataFrame. Scoring Class, the match columns in STRUCTURAL_COLS (True, False or None),
        and Aligned Position: the homolog residue that aligns to the query's residue.
    '''
    def pair(col):
        return [_match(a, b) for a, b in zip(df[f'{col}{q}'], df[f'{col}{c}'])]

    out = pd.DataFrame(index=df.index)
    out['Scoring Class'] = df[f'Molecular Consequence{q}'].apply(scoring_class)
    out['Biotype Match'] = pair('Biotype')
    # The main consequence terms are compared, so a modifier such as NMD_transcript_variant
    # on one transcript doesn't make a missense variant a mismatch
    q_terms = df[f'Molecular Consequence{q}'].map(_main_terms)
    c_terms = df[f'Molecular Consequence{c}'].map(_main_terms)
    out['Consequence Match'] = [None if not a and not b else bool(a & b) for a, b in zip(q_terms, c_terms)]
    out['AA Match'] = [_match(a, b) for a, b in zip(_aa_change(df, q), _aa_change(df, c))]
    # Residue numbers differ between orthologs and isoforms, so positions are compared through a protein alignment
    aligned, position_matches = _aligned_positions(df, q, c, q_species, c_species)
    out['AA Position Match'] = position_matches
    out['Domain Match'], domain_labels = _domain_matches(df, q, c, q_species, c_species)
    out['Intron Match'] = _intron_matches(df, q, c, q_species, c_species)
    out['Intron Offset Match'] = pair('Intron Offset')

    # Blank the criteria that belong to the other class
    for cls, criteria in CLASS_CRITERIA.items():
        other = out['Scoring Class'] != cls
        for col in criteria:
            out[col] = out[col].astype(object).where(~other, None)

    out['Aligned Position'] = aligned

    # How well the proteins align around the residue: an identical change at a weakly aligned
    # residue may not be at the equivalent position
    out['Local Alignment Identity'] = _local_identities(df, q, c, q_species, c_species, position_matches)
    uncertain = [v is not None and v < MIN_LOCAL_IDENTITY for v in out['Local Alignment Identity']]

    detail = pd.DataFrame({'q_terms': q_terms.values, 'c_terms': c_terms.values, 'aligned': aligned,
                           'uncertain': uncertain,
                           'c_pos': df[f'Amino Acids{c}'].map(_residue_position).values,
                           'q_off': df[f'Intron Offset{q}'].values, 'c_off': df[f'Intron Offset{c}'].values},
                          index=df.index)
    tiers = [_tier(cls, {**m, **d}) for cls, m, d in zip(out['Scoring Class'], out.to_dict('records'),
                                                          detail.to_dict('records'))]
    out.insert(1, 'Equivalency Tier', tiers)
    out.insert(2, 'Tier Basis', [TIER_BASIS.get(cls, {}).get(t, '') for cls, t in zip(out['Scoring Class'], tiers)])
    # A truncation near, not at, the aligned residue
    near = [t == 2 and not m for t, m in zip(tiers, out['AA Position Match'])]
    out.loc[near, 'Tier Basis'] = f'Truncation within {TRUNCATION_WINDOW} residues'
    # Pairs that rest on the aligned residue say when the alignment there is weak
    weak = [u and t in (1, 2) and m for u, t, m in zip(uncertain, tiers, out['AA Position Match'])]
    out.loc[weak, 'Tier Basis'] = [('Identical change' if a else 'Same residue') + ', alignment uncertain'
                                   for a in out.loc[weak, 'AA Match']]
    # Tier 1 names the domain the change is in
    in_domain = [t == 1 and bool(d) for t, d in zip(tiers, domain_labels)]
    out.loc[in_domain, 'Tier Basis'] = [f'Identical change in {d}' for d, keep in zip(domain_labels, in_domain) if keep]

    # NMD and non-coding transcripts are shown but don't set the allele's tier
    out['Counts Toward Tier'] = (df[f'Biotype{c}'] == 'protein_coding').values
    return out.astype(object)


def add_best_tier(score_df, key):
    '''
    Best Tier: each allele's (or ClinVar variant's) best tier over its protein-coding
    transcripts, or over all its transcripts when none is protein-coding.
    '''
    tier = pd.to_numeric(score_df['Equivalency Tier'], errors='coerce')
    counted = tier.where(score_df['Counts Toward Tier'].astype(bool))
    best = counted.groupby(score_df[key]).transform('min')
    score_df['Best Tier'] = best.fillna(tier.groupby(score_df[key]).transform('min'))
    return score_df


def total_score(score_df, match_cols):
    '''Percent of the applicable criteria (not None) that match.'''
    values = score_df[match_cols].apply(lambda s: s.map(lambda v: np.nan if v is None or pd.isna(v) else float(v)))
    return (values.mean(axis=1) * 100).round(1)
