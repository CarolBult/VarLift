'''
Map a residue position from one protein to its homologous position in another.

Human and mouse orthologs, and the isoforms of either, can differ in length
before a variant (e.g. MECP2 e1: human R145 is mouse R150), so raw residue
numbers can't be compared directly. Proteins are aligned once per transcript
pair with a global alignment; very long proteins (e.g. titin) are aligned
through a window around the position instead.
'''
import yaml
import logging
from functools import lru_cache
from pathlib import Path

import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices


with open('app/config.yaml') as f:
    config = yaml.safe_load(f)

# Proteins longer than this (product of lengths) are aligned through a window
MAX_FULL_ALIGNMENT_CELLS = 25_000_000
WINDOW = 50

_blosum62 = substitution_matrices.load('BLOSUM62')


def _aligner(mode):
    aligner = PairwiseAligner(mode=mode, substitution_matrix=_blosum62,
                              open_gap_score=-10, extend_gap_score=-0.5)
    if mode == 'global':
        # End gaps cost a little more than internal ones, so a gap inside a repeat
        # (e.g. MECP2's N-terminal alanines) isn't pushed past the start methionine
        aligner.open_end_gap_score = -11
        aligner.extend_end_gap_score = -0.5
    return aligner


_global_aligner = _aligner('global')
_local_aligner = _aligner('local')


@lru_cache(maxsize=2)
def _sequences(species):
    '''Transcript ID -> protein sequence, or None when the parquet hasn't been built.'''
    path = config['paths'][f'{species}_protein_pqt']
    if not Path(path).exists():
        logging.warning(f'{path} not found; run app/data_refresh.py to build it. '
                        'Protein positions are compared by residue number instead.')
        return None
    df = pd.read_parquet(path)
    return dict(zip(df['transcript_id'], df['sequence']))


def _species(name):
    return 'human' if name in ('human', 'homo_sapiens') else 'mouse'


def _map_through(alignment, pos0):
    '''Position in the second sequence aligned to 0-based pos0 in the first, or None in a gap.'''
    for (a_start, a_end), (b_start, _) in zip(*alignment.aligned):
        if a_start <= pos0 < a_end:
            return b_start + (pos0 - a_start)
    return None


@lru_cache(maxsize=4096)
def _full_alignment(seq_a, seq_b):
    return _global_aligner.align(seq_a, seq_b)[0]


@lru_cache(maxsize=65536)
def _windowed_position(seq_a, seq_b, pos0):
    start = max(0, pos0 - WINDOW)
    window = seq_a[start:pos0 + WINDOW + 1]
    alignment = _local_aligner.align(window, seq_b)[0]
    return _map_through(alignment, pos0 - start)


class SequenceUnavailable(Exception):
    '''Raised when a transcript's protein sequence isn't available.'''


def homologous_position(tx_a, species_a, tx_b, species_b, pos):
    '''
    The residue in transcript tx_b's protein that aligns to residue pos of tx_a's protein.

    Parameters:
        tx_a, tx_b: str. Ensembl transcript IDs.
        species_a, species_b: str. 'human'/'homo_sapiens' or 'mouse'/'mus_musculus'.
        pos: int. 1-based residue position in tx_a's protein.

    Returns:
        int. 1-based residue position in tx_b's protein, or None when it falls in an alignment gap.

    Raises:
        SequenceUnavailable: when either protein sequence is missing.
    '''
    seqs_a, seqs_b = _sequences(_species(species_a)), _sequences(_species(species_b))
    seq_a = seqs_a.get(tx_a) if seqs_a else None
    seq_b = seqs_b.get(tx_b) if seqs_b else None
    if not seq_a or not seq_b or not 0 < pos <= len(seq_a):
        raise SequenceUnavailable(f'{tx_a} or {tx_b}')

    if len(seq_a) * len(seq_b) <= MAX_FULL_ALIGNMENT_CELLS:
        mapped = _map_through(_full_alignment(seq_a, seq_b), pos - 1)
    else:
        mapped = _windowed_position(seq_a, seq_b, pos - 1)

    return None if mapped is None else mapped + 1
