'''
Lookup HP (Human Phenotype Ontology) terms for human variants from ClinVar's
variant_summary.tsv — a much richer source than CLNDISDB.

Lazy-loaded at first use. Indexed by (Chromosome, Start) since ClinVar often
records ref/alt as 'na' for protein-level variants.
'''
import re
import yaml
import logging
import pandas as pd

with open('app/config.yaml') as f:
    _config = yaml.safe_load(f)

_HP_PATTERN = re.compile(r'Human Phenotype Ontology:(HP:\d+)')

# Lazy singleton: (chrom:str, start:str) -> set of 'HP:...' strings
_HP_INDEX = None


def _build_index():
    '''Load variant_summary.tsv and build a {(chrom, start): set_of_HPs} map.'''
    path = _config['paths']['ncbi_var_summary_grch38']
    logging.info(f'Loading ClinVar variant summary for HP lookup: {path}')

    # Read only the columns we need to keep memory down
    df = pd.read_csv(
        path,
        sep='\t',
        usecols=['Chromosome', 'Start', 'PhenotypeIDS'],
        dtype={'Chromosome': str, 'Start': str, 'PhenotypeIDS': str},
    )

    df = df.dropna(subset=['PhenotypeIDS'])
    df = df[df['PhenotypeIDS'].str.contains('Human Phenotype Ontology:HP:', na=False)]

    index = {}
    for chrom, start, pids in df[['Chromosome', 'Start', 'PhenotypeIDS']].itertuples(index=False):
        hps = set(_HP_PATTERN.findall(pids))
        if not hps:
            continue
        key = (str(chrom), str(start))
        if key in index:
            index[key].update(hps)
        else:
            index[key] = hps

    logging.info(f'ClinVar HP lookup indexed {len(index)} positions')
    return index


def lookup_hp_ids(chromosome, start):
    '''
    Return sorted list of HP IDs for a variant at (chromosome, start).
    Returns an empty list if the variant isn't in the index.
    '''
    global _HP_INDEX
    if _HP_INDEX is None:
        _HP_INDEX = _build_index()
    hps = _HP_INDEX.get((str(chromosome), str(start)), set())
    return sorted(hps)
