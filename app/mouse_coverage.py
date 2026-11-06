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
from app.phenotype_util import MGI_phenotype_df


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


def mgi_alleles_at(variants):
    '''
    MGI alleles at exactly the submitted position and alleles, with their
    disease models and phenotypes.
    '''
    va = mus_alleles_df.dropna(subset=['StartPosition'])
    va = va[['AlleleId', 'AlleleSymbol', 'Chromosome', 'StartPosition', 'SequenceOfReference', 'SequenceOfVariant']].copy()
    va['Chromosome'] = va['Chromosome'].astype(str)
    va['StartPosition'] = va['StartPosition'].astype(int)

    keys = variants[['Name', 'Chromosome', 'Start', 'Ref', 'Alt']].copy()
    keys['Chromosome'] = keys['Chromosome'].astype(str)
    keys['Ref'] = keys['Ref'].str.upper()
    keys['Alt'] = keys['Alt'].str.upper()

    matched = keys.merge(va, left_on=['Chromosome', 'Start', 'Ref', 'Alt'],
                         right_on=['Chromosome', 'StartPosition', 'SequenceOfReference', 'SequenceOfVariant'])
    matched = matched[['Name', 'AlleleId', 'AlleleSymbol']].drop_duplicates()
    matched = matched.rename(columns={'Name': 'Input Name', 'AlleleId': 'AlleleID'})

    diseases = (MGI_disease_df.dropna(subset=['AlleleID'])
                .groupby('AlleleID')['DOterm'].apply(lambda s: sorted(set(s.dropna()))).to_dict())
    matched['Associated Diseases'] = matched['AlleleID'].map(lambda a: diseases.get(a, []))

    phenotypes = MGI_phenotype_df['Phenotypes'].str.split('|').to_dict()
    matched['Phenotypes'] = matched['AlleleID'].map(lambda a: phenotypes.get(a, []))

    return matched


def query_allele_coverage(genes, submitted_allele_ids):
    '''
    Coverage table for a mouse query: every disease-associated MGI allele of
    each queried gene, flagging the alleles that match a submitted variant,
    with each allele's phenotypes.
    '''
    frames = [allele_coverage_for_gene(gene, submitted_allele_ids) for gene in genes]
    coverage = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if coverage.empty:
        return coverage

    # Reason explains why VarLift couldn't evaluate an allele, which doesn't apply here
    coverage = coverage.drop(columns=['Reason']).rename(columns={'Evaluated': 'Submitted'})

    phenotypes = MGI_phenotype_df['Phenotypes'].str.split('|').to_dict()
    coverage['Phenotypes'] = coverage['AlleleID'].map(lambda a: phenotypes.get(a, []))
    return coverage


def query_other_alleles(genes, submitted_allele_ids):
    '''
    MGI alleles of each queried gene that have phenotypes but no disease model,
    so they are outside the coverage table.
    '''
    disease_ids = set(_load_disease_models()['AlleleID'])
    phenotypes = MGI_phenotype_df['Phenotypes'].str.split('|').to_dict()

    va = _load_variant_allele()
    va = va[va['AlleleAssociatedGeneSymbol'].isin(genes)]
    va = va[['AlleleAssociatedGeneSymbol', 'AlleleId', 'AlleleSymbol']].drop_duplicates(subset=['AlleleId'])
    va = va[~va['AlleleId'].isin(disease_ids) & va['AlleleId'].isin(phenotypes.keys())]

    other = va.rename(columns={'AlleleAssociatedGeneSymbol': 'Gene Symbol', 'AlleleId': 'AlleleID'})
    other['Submitted'] = other['AlleleID'].isin(set(submitted_allele_ids))
    other['Phenotypes'] = other['AlleleID'].map(phenotypes)
    return other.sort_values(['Submitted', 'AlleleSymbol'], ascending=[False, True]).reset_index(drop=True)
