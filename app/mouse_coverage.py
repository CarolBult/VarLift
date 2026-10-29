'''
Mouse allele coverage: for a gene, list every MGI allele with disease
annotations and indicate whether VarLift's pipeline was able to evaluate it.

Primary source: MGI_DiseaseMouseModel.rpt — the authoritative list of all
mouse models associated with a disease, including knockouts, gene traps, and
transgenics that have no variant-level annotation. VARIANT-ALLELE-MGI.tsv.gz
is used as a secondary source for consequence information when available.
'''
import logging
import pandas as pd

_VARIANT_ALLELE = None
_DISEASE_MODELS = None


def _load_variant_allele():
    global _VARIANT_ALLELE
    if _VARIANT_ALLELE is not None:
        return _VARIANT_ALLELE
    path = 'app/data/MGI/VARIANT-ALLELE-MGI.tsv.gz'
    logging.info(f'Loading MGI variant-allele table: {path}')
    df = pd.read_csv(path, sep='\t', skiprows=14, low_memory=False)
    _VARIANT_ALLELE = df
    return df


def _load_disease_models():
    '''
    Load MGI_DiseaseMouseModel.rpt. The file has no header row; comment
    lines at top start with '#'. Columns we need (0-indexed):
      0: disease name
      5: allele symbol (e.g. 'Mecp2<tm1Pplt>')
      6: allele MGI ID
      10: gene symbol
    '''
    global _DISEASE_MODELS
    if _DISEASE_MODELS is not None:
        return _DISEASE_MODELS
    path = 'app/data/MGI/MGI_DiseaseMouseModel.rpt'
    logging.info(f'Loading MGI disease-model table: {path}')
    df = pd.read_csv(
        path, sep='\t', comment='#', header=None,
        usecols=[0, 5, 6, 10],
        names=['Disease', 'AlleleSymbol', 'AlleleID', 'GeneSymbol'],
        dtype=str, low_memory=False,
    )
    df = df.dropna(subset=['AlleleID', 'GeneSymbol'])
    _DISEASE_MODELS = df
    return df


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
        consequence = va_by_allele.get(aid, '')
        is_evaluated = aid in evaluated_set
        rows.append({
            'Gene Symbol': gene_symbol,
            'AlleleID': aid,
            'AlleleSymbol': symbol,
            'Consequence': consequence if isinstance(consequence, str) else '',
            'Diseases': diseases_by_allele[aid],
            'Evaluated': is_evaluated,
            'Reason': _reason_for(consequence, is_evaluated),
        })
    return pd.DataFrame(rows, columns=cols)
