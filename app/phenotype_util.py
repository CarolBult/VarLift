import logging
from collections import defaultdict
from pathlib import Path

import yaml
import pandas as pd


# Load Config
with open('app/config.yaml') as f:
    config = yaml.safe_load(f)


# ALL_Phenotype.rpt lists every MP term annotated to any genotype carrying the allele,
# including "normal" (NOT) annotations, with no way to tell them apart: Acvr1<tm1Glh>'s only
# term, "skeleton phenotype", records that its skeleton was examined and found normal.
# MGI_PhenoGenoMP.rpt lists genotype annotations without the normal ones, so an allele's
# term is treated as normal when none of its genotypes is annotated with that term, or with
# a more specific term beneath it, as abnormal. The second test keeps the few terms
# ALL_Phenotype lists above an allele's abnormal annotations (e.g. "abnormal trabecular bone
# morphology" above "increased bone trabecular spacing") from being called normal. It is
# applied only below MP's top-level system terms, which sit above nearly every annotation.
# A term found normal in one genotype and abnormal in another is abnormal.

def _mp_parents():
    '''{MP ID: [parent MP IDs]} from mp.obo, or {} if it is missing.'''
    path = Path(config['paths']['mp_obo'])
    if not path.exists():
        return {}
    parents, term = defaultdict(list), None
    with open(path, encoding='utf-8') as f:
        for line in f:
            if line.startswith('['):
                term = None
            elif line.startswith('id: MP:'):
                term = line[4:].strip()
            elif line.startswith('is_a: MP:') and term:
                parents[term].append(line[6:16])
    return parents


def _abnormal_terms_by_allele():
    '''{allele symbol: MP IDs annotated as abnormal to any genotype carrying it}, or None.'''
    path = Path(config['paths']['mgi_genotype_phenotype'])
    if not path.exists():
        logging.warning('MGI_PhenoGenoMP.rpt not found (python app/data_refresh.py fetches it); '
                        'normal phenotype annotations will be treated as abnormal')
        return None
    geno = pd.read_csv(path, sep='\t', header=None, usecols=[1, 3], dtype=str,
                       names=['Alleles', 'MP'])
    abnormal = defaultdict(set)
    for alleles, mp in geno.itertuples(index=False):
        for symbol in str(alleles).split('|'):
            abnormal[symbol].add(mp)
    return abnormal


def _split_normal(phenotype_df):
    '''Split each allele's Phenotypes into abnormal (Phenotypes) and Normal Phenotypes lists.'''
    phenotypes = phenotype_df['Phenotypes'].fillna('').str.split('|')
    phenotypes = phenotypes.map(lambda terms: [t for t in terms if t])
    abnormal_by_allele = _abnormal_terms_by_allele()
    if abnormal_by_allele is None:
        return phenotypes, phenotypes.map(lambda _: [])

    parents = _mp_parents()
    top_level = {t for t, ps in parents.items() if 'MP:0000001' in ps} | {'MP:0000001'}
    ancestor_cache = {}

    def ancestors(term):
        if term not in ancestor_cache:
            out = {term}
            for p in parents.get(term, ()):
                out |= ancestors(p)
            ancestor_cache[term] = frozenset(out)
        return ancestor_cache[term]

    abnormal_col, normal_col = [], []
    for symbol, terms in zip(phenotype_df['Allele Symbol'], phenotypes):
        annotated = abnormal_by_allele.get(symbol, set())
        above_annotated = set().union(*(ancestors(t) for t in annotated)) - top_level if annotated else set()
        abnormal, normal = [], []
        for t in terms:
            mp = t.split(',', 1)[0]
            (abnormal if mp in annotated or mp in above_annotated else normal).append(t)
        abnormal_col.append(abnormal)
        normal_col.append(normal)
    return (pd.Series(abnormal_col, index=phenotype_df.index),
            pd.Series(normal_col, index=phenotype_df.index))


MGI_phenotype_df = pd.read_csv(config['paths']['mgi_phenotype'], sep='\t')
MGI_phenotype_df['Phenotypes'], MGI_phenotype_df['Normal Phenotypes'] = _split_normal(MGI_phenotype_df)
MGI_phenotype_df = MGI_phenotype_df.set_index('Allele ID')
logging.info(f"MGI phenotypes: {MGI_phenotype_df['Normal Phenotypes'].map(len).sum()} normal annotations "
             f"on {(MGI_phenotype_df['Normal Phenotypes'].map(len) > 0).sum()} alleles kept apart")


def allele_phenotypes(allele_id):
    '''(abnormal, normal) "MP:...,label" lists for one MGI allele ID.'''
    if allele_id in MGI_phenotype_df.index:
        row = MGI_phenotype_df.loc[allele_id]
        return list(row['Phenotypes']), list(row['Normal Phenotypes'])
    return [], []


def allele_phenotype_match(allele_df):
    '''
    Pairs MGI alleles with their detailed list of phenotypes.

    Parameters:
        allele_df. DataFrame. Expects AlleleID, AlleleSymbol, Chromosome, Start, End, Ref, Alt, Molecular Consequence.

    Returns:
        DataFrame. Contains AlleleID, AlleleSymbol, Phenotypes (abnormal) and Normal Phenotypes,
        each a list of "MP:...,label" strings.
    '''
    phenotype_df = allele_df.join(MGI_phenotype_df, on='AlleleID', how='inner')

    phenotype_df.drop(['Allele Symbol', 'Chromosome', 'Start', 'End', 'Ref', 'Alt', 'Molecular Consequence'], inplace=True, axis=1)

    return phenotype_df
