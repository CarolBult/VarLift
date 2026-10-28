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
_DISEASE_PATTERN = re.compile(r'(MONDO:MONDO:\d+|OMIM:\d+|Orphanet:\d+)')

# Lazy singleton: (chrom:str, start:str) -> set of 'HP:...' strings

_HP_INDEX = None
_DISEASE_INDEX = None

def _build_indexes():
    '''
    Scan variant_summary.tsv once, build two indexes keyed by (chrom, start):
      - HP_INDEX: direct HP cross-refs (when PhenotypeIDS has 'Human Phenotype Ontology:HP:...')
      - DISEASE_INDEX: disease IDs (MONDO/OMIM/Orphanet) that can be chained through hpoa
    '''
    path = _config['paths']['ncbi_var_summary_grch38']
    logging.info(f'Loading ClinVar variant summary: {path}')

    df = pd.read_csv(
        path,
        sep='\t',
        usecols=['Chromosome', 'Start', 'PhenotypeIDS'],
        dtype={'Chromosome': str, 'Start': str, 'PhenotypeIDS': str},
    )
    df = df.dropna(subset=['PhenotypeIDS'])

    hp_index = {}
    disease_index = {}
    for chrom, start, pids in df[['Chromosome', 'Start', 'PhenotypeIDS']].itertuples(index=False):
        key = (str(chrom), str(start))

        hps = set(_HP_PATTERN.findall(pids))
        if hps:
            if key in hp_index:
                hp_index[key].update(hps)
            else:
                hp_index[key] = hps

        diseases = set()
        for m in _DISEASE_PATTERN.findall(pids):
            # Normalize 'MONDO:MONDO:...' -> 'MONDO:...'
            if m.startswith('MONDO:MONDO:'):
                m = m[6:]
            diseases.add(m)
        if diseases:
            if key in disease_index:
                disease_index[key].update(diseases)
            else:
                disease_index[key] = diseases

    logging.info(f'Index built: {len(hp_index)} HP positions, {len(disease_index)} disease positions')
    return hp_index, disease_index

def _ensure_loaded():
    global _HP_INDEX, _DISEASE_INDEX
    if _HP_INDEX is None or _DISEASE_INDEX is None:
        _HP_INDEX, _DISEASE_INDEX = _build_indexes()

def lookup_hp_ids(chromosome, start):
    '''
    Return sorted list of HP IDs for a variant at (chromosome, start).
    Returns an empty list if the variant isn't in the index.
    '''
    _ensure_loaded()
    hps = _HP_INDEX.get((str(chromosome), str(start)), set())
    return sorted(hps)

def lookup_disease_ids(chromosome, start):
    '''
    Return sorted list of disease IDs (MONDO, OMIM, Orphanet) extracted from
    variant_summary.tsv's PhenotypeIDS for the variant at (chromosome, start).
    Used to chain through phenotype.hpoa for HP lookup.
    '''
    _ensure_loaded()
    diseases = _DISEASE_INDEX.get((str(chromosome), str(start)), set())
    return sorted(diseases)
