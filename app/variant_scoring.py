'''
Structural similarity scoring between a query variant and its homologous variants.

Each query variant is scored by the class of its molecular consequence:
    - Coding: amino-acid change, protein position, exon and Pfam domain.
    - Splice/Intronic: intron number and the variant's offset from the nearest
      exon boundary (as in HGVS c.123+5).
Biotype and consequence are compared for both classes. A criterion that does
not apply to the query's class, or that neither variant has a value for
(e.g. neither is in a Pfam domain), is left as None and does not count
towards the Total Score.
'''
import yaml
import logging
import numpy as np
import pandas as pd
from pathlib import Path


with open('app/config.yaml') as f:
    config = yaml.safe_load(f)


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
CLASS_CRITERIA = {CODING: ['AA Match', 'AA Position Match', 'Exon Match', 'Domain Match'],
                  SPLICE: ['Intron Match', 'Intron Offset Match']}

STRUCTURAL_COLS = ['Biotype Match', 'Consequence Match', 'AA Match', 'AA Position Match',
                   'Exon Match', 'Domain Match', 'Intron Match', 'Intron Offset Match']


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


def structural_matches(df, q, c):
    '''
    Compare each query/homolog pair in a merged DataFrame.

    Parameters:
        df: DataFrame. One row per pair, query columns suffixed q and homolog columns suffixed c.
        q: str. Suffix of the query variant's columns.
        c: str. Suffix of the homologous variant's columns.

    Returns:
        DataFrame. Scoring Class and the match columns in STRUCTURAL_COLS (True, False or None).
    '''
    def pair(col):
        return [_match(a, b) for a, b in zip(df[f'{col}{q}'], df[f'{col}{c}'])]

    out = pd.DataFrame(index=df.index)
    out['Scoring Class'] = df[f'Molecular Consequence{q}'].apply(scoring_class)
    out['Biotype Match'] = pair('Biotype')
    out['Consequence Match'] = pair('Molecular Consequence')
    out['AA Match'] = [_match(a, b) for a, b in zip(_aa_change(df, q), _aa_change(df, c))]
    out['AA Position Match'] = pair('Amino Acids')
    out['Exon Match'] = pair('Exon Rank')
    out['Domain Match'] = pair('Pfam Domain ID')
    out['Intron Match'] = pair('Intron Rank')
    out['Intron Offset Match'] = pair('Intron Offset')

    # Blank the criteria that belong to the other class
    for cls, criteria in CLASS_CRITERIA.items():
        other = out['Scoring Class'] != cls
        for col in criteria:
            out[col] = out[col].astype(object).where(~other, None)

    return out.astype(object)


def total_score(score_df, match_cols):
    '''Percent of the applicable criteria (not None) that match.'''
    values = score_df[match_cols].apply(lambda s: s.map(lambda v: np.nan if v is None or pd.isna(v) else float(v)))
    return (values.mean(axis=1) * 100).round(1)
