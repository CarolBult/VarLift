'''
Mouse allele coverage: for a gene, list every MGI allele with disease
annotations and indicate whether VarLift's pipeline was able to evaluate it.

Primary source: MGI_DiseaseMouseModel.rpt — the authoritative list of all
mouse models associated with a disease, including knockouts, gene traps, and
transgenics that have no variant-level annotation. VARIANT-ALLELE-MGI.tsv.gz
is used as a secondary source for consequence information when available.
'''
import pandas as pd
import re

from app.gene_util import mus_alleles_df
from app.disease_util import MGI_disease_df


def _load_variant_allele():
    '''VARIANT-ALLELE-MGI table, already loaded by gene_util.'''
    return mus_alleles_df


def _load_disease_models():
    '''
    MGI_DiseaseMouseModel.rpt, already loaded by disease_util, reduced to the
    disease name, allele symbol (e.g. 'Mecp2<tm1Pplt>'), allele MGI ID and
    gene symbol.
    '''
    df = MGI_disease_df[['DOterm', 'AlleleSymbol', 'AlleleID', 'MarkerSymbol']]
    df = df.rename(columns={'DOterm': 'Disease', 'MarkerSymbol': 'GeneSymbol'})
    return df.dropna(subset=['AlleleID', 'GeneSymbol'])


def _unique_consequences(consequence):
    '''
    VARIANT-ALLELE-MGI repeats a consequence once per transcript
    (e.g. 'missense_variant,missense_variant,...'). Keep each one once, in
    order, separated by ', ' so table cells can wrap.
    '''
    if not isinstance(consequence, str):
        return ''
    return ', '.join(dict.fromkeys(c.strip() for c in consequence.split(',') if c.strip()))


def _reason_for(consequence, is_evaluated):
    if is_evaluated:
        return ''
    if not isinstance(consequence, str) or consequence == '' or pd.isna(consequence):
        return 'Structural or transgenic allele (no amino-acid-level annotation)'
    low = consequence.lower()
    if 'missense' in low or 'stop_gained' in low:
        return 'Variant could not be cross-mapped to human'
    return f'Consequence not evaluable: {consequence}'


def allele_coverage_for_gene(gene_symbol, evaluated_allele_ids):
    '''
    List all disease-associated MGI alleles for a mouse gene, with a per-allele
    status indicating whether VarLift's pipeline evaluated it.
    '''
    dm = _load_disease_models()
    va = _load_variant_allele()

    cols = ['Gene Symbol', 'AlleleID', 'AlleleSymbol', 'Consequence',
            'Diseases', 'Evaluated', 'Reason']

    gene_rows = dm[dm['GeneSymbol'] == gene_symbol]
    if gene_rows.empty:
        return pd.DataFrame(columns=cols)

    diseases_by_allele = (
        gene_rows.groupby('AlleleID')['Disease']
        .apply(lambda s: sorted(set(s.dropna())))
        .to_dict()
    )
    symbols_by_allele = (
        gene_rows.drop_duplicates(subset=['AlleleID'])
        .set_index('AlleleID')['AlleleSymbol']
        .to_dict()
    )
    va_by_allele = (
        va[['AlleleId', 'MostSevereConsequenceName']]
        .dropna(subset=['AlleleId']).drop_duplicates(subset=['AlleleId'])
        .set_index('AlleleId')['MostSevereConsequenceName']
        .to_dict()
    )

    evaluated_set = set(evaluated_allele_ids or [])
    rows = []
    for aid in sorted(diseases_by_allele.keys()):
        symbol = symbols_by_allele.get(aid, '')
        consequence = _unique_consequences(va_by_allele.get(aid, ''))
        is_evaluated = aid in evaluated_set
        rows.append({
            'Gene Symbol': gene_symbol,
            'AlleleID': aid,
            'AlleleSymbol': re.sub(r'<([^>]+)>', r'<sup>\1</sup>', symbol) if isinstance(symbol, str) else symbol,
            'Consequence': consequence,
            'Diseases': diseases_by_allele[aid],
            'Evaluated': is_evaluated,
            'Reason': _reason_for(consequence, is_evaluated),
        })
    return pd.DataFrame(rows, columns=cols)
