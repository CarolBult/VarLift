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
from functools import lru_cache

from app.gene_util import mus_alleles_df, mouse_gene_pqt, config
from app.data_util import count_header_lines
from app.disease_util import MGI_disease_df, fetch_mus_doid, map_doids_to_mondo
from app.phenotype_util import allele_phenotypes


@lru_cache(maxsize=1)
def _markers():
    '''MRK_List2.rpt: current MGI markers with GRCm39 coordinates and marker type.'''
    df = pd.read_csv(config['paths']['mgi_markers'], sep='\t', dtype=str)
    df = df[df['Status'] == 'O']
    df = df.rename(columns={'Marker Symbol': 'Symbol', 'Marker Type': 'Type',
                            'genome coordinate start': 'Start', 'genome coordinate end': 'End'})
    df['Start'] = pd.to_numeric(df['Start'], errors='coerce')
    df['End'] = pd.to_numeric(df['End'], errors='coerce')
    return df[['Symbol', 'Chr', 'Start', 'End', 'Type']]


@lru_cache(maxsize=1)
def _phenotypic_alleles():
    '''
    MGI_PhenotypicAllele.rpt: every current MGI allele with its marker and allele
    type. Unlike VARIANT-ALLELE-MGI it includes alleles of gene clusters (e.g.
    Pcdhg) and alleles added since the Alliance release.
    '''
    path = config['paths']['mgi_phenotypic_alleles']
    df = pd.read_csv(path, sep='\t', header=None, dtype=str, quoting=3, skiprows=count_header_lines(path),
                     names=['AlleleID', 'AlleleSymbol', 'AlleleName', 'AlleleType', 'AlleleAttribute', 'PubMedID',
                            'MarkerID', 'MarkerSymbol', 'RefSeq', 'Ensembl', 'MP', 'Synonyms', 'MarkerName'])
    df['Allele Type'] = [', '.join(p for p in [t, a.replace('|', ', ') if isinstance(a, str) else ''] if isinstance(p, str) and p)
                         for t, a in zip(df['AlleleType'], df['AlleleAttribute'])]
    return df[['AlleleID', 'AlleleSymbol', 'MarkerSymbol', 'Allele Type']]


def related_markers(gene, positions):
    '''
    MGI markers other than the gene itself that a submitted variant also falls in:
    genes with an exon at a submitted position (overlapping or nested genes) and
    the gene cluster the gene belongs to. MGI gives most clusters no
    coordinates and publishes no membership file, so a gene is taken to belong to
    a Complex/Cluster/Region marker on its chromosome whose symbol starts its own
    (Pcdhg -> Pcdhga3, Hoxa -> Hoxa1).

    Parameters:
        gene: str. Queried mouse gene symbol.
        positions: iterable of (chromosome, position) of the submitted variants.

    Returns:
        list of str. Marker symbols, clusters first.
    '''
    markers = _markers()
    own = markers[markers['Symbol'] == gene]
    if own.empty:
        return []
    chrom = own['Chr'].iloc[0]

    clusters = markers[(markers['Type'] == 'Complex/Cluster/Region') & (markers['Chr'] == chrom)]
    clusters = [s for s in clusters['Symbol'] if gene.startswith(s) and gene != s]

    # Positions inside the gene itself, so other variants in the batch don't pull in unrelated genes
    s, e = own['Start'].iloc[0], own['End'].iloc[0]
    inside = [p for c, p in positions if str(c) == chrom and pd.notna(s) and s <= p <= e]
    genes = markers[(markers['Type'] == 'Gene') & (markers['Chr'] == chrom) & (markers['Symbol'] != gene)]
    spanning = {sym for p in inside for sym in genes.loc[(genes['Start'] <= p) & (genes['End'] >= p), 'Symbol']}
    # A gene whose span merely contains the position (e.g. in an intron, as Pcdhga1 does for a
    # Pcdhga3 variable-exon variant) is a different gene there; count it only where the
    # position is in one of its exons
    overlapping = sorted(_with_exon_at(spanning, chrom, inside))
    return clusters + overlapping


def _with_exon_at(symbols, chrom, positions):
    '''
    The subset of mouse gene symbols with an exon (of any Ensembl transcript)
    containing one of the positions. Reads the exon table once for all genes.
    '''
    if not symbols or not positions:
        return set()
    genes = mouse_gene_pqt.loc[mouse_gene_pqt['Name'].isin(symbols), ['Name', 'gene_id']]
    if genes.empty:
        return set()
    exons = pd.read_parquet(config['paths']['mouse_exon_pqt'], columns=['seqid', 'start', 'end', 'gene_id'],
                            filters=[('gene_id', 'in', genes['gene_id'].tolist())])
    exons = exons[exons['seqid'].astype(str) == str(chrom)]
    hit = pd.Series(False, index=exons.index)
    for p in positions:
        hit |= (exons['start'] <= p) & (exons['end'] >= p)
    hit_ids = set(exons.loc[hit, 'gene_id'])
    return set(genes.loc[genes['gene_id'].isin(hit_ids), 'Name'])


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
    # One generic reason: MGI's allele type doesn't tell evaluable from unevaluable
    # alleles (engineered and humanized alleles can carry evaluable changes), and the
    # Consequence column already shows the consequence where there is one
    if is_evaluated:
        return ''
    return ('Variant genome coordinates not available or molecular consequence '
            'not evaluable (e.g. 5_prime_UTR_variant).')


def allele_coverage_for_gene(gene_symbol, evaluated_allele_ids, all_alleles=False):
    '''
    List all disease-associated MGI alleles for a mouse gene, with a per-allele
    status indicating whether VarLift's pipeline evaluated it.

    With all_alleles (human queries), alleles with variant data but no disease
    model are listed too, since they are scored as well, and each allele's
    phenotypes are added.
    '''
    dm = _load_disease_models()
    va = _load_variant_allele()

    cols = ['Gene Symbol', 'AlleleID', 'AlleleSymbol', 'Consequence',
            'Diseases', 'Evaluated', 'Reason']

    gene_rows = dm[dm['GeneSymbol'] == gene_symbol]
    gene_va = va[va['AlleleAssociatedGeneSymbol'] == gene_symbol] if all_alleles else va.iloc[0:0]
    if gene_rows.empty and gene_va.empty:
        return pd.DataFrame(columns=cols)

    diseases_by_allele = (
        gene_rows.groupby('AlleleID')['Disease']
        .apply(lambda s: sorted(set(s.dropna())))
        .to_dict()
    )
    symbols_by_allele = (
        gene_va.drop_duplicates(subset=['AlleleId'])
        .set_index('AlleleId')['AlleleSymbol']
        .to_dict()
    )
    symbols_by_allele.update(
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
    for aid in sorted(symbols_by_allele.keys()):
        symbol = symbols_by_allele.get(aid, '')
        consequence = _unique_consequences(va_by_allele.get(aid, ''))
        is_evaluated = aid in evaluated_set
        rows.append({
            'Gene Symbol': gene_symbol,
            'AlleleID': aid,
            # Disease-model symbols use <...>; VARIANT-ALLELE-MGI's already have <sup> tags
            'AlleleSymbol': re.sub(r'<([^>]+)>', r'<sup>\1</sup>', symbol) if isinstance(symbol, str) and '<sup>' not in symbol else symbol,
            'Consequence': consequence,
            'Diseases': diseases_by_allele.get(aid, []),
            'Evaluated': is_evaluated,
            'Reason': _reason_for(consequence, is_evaluated),
        })
    coverage = pd.DataFrame(rows, columns=cols)

    if all_alleles:
        coverage['Phenotypes'] = coverage['AlleleID'].map(lambda a: allele_phenotypes(a)[0])
        coverage['Normal Phenotypes'] = coverage['AlleleID'].map(lambda a: allele_phenotypes(a)[1])
    return coverage


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

    # Disease models as MONDO terms with their IDs, as for the alleles scored against human variants
    doid_map = fetch_mus_doid(matched['AlleleID'].unique())
    mondo_ids, mondo_terms = map_doids_to_mondo(doid_map) if doid_map else ({}, {})
    matched['Associated Diseases'] = matched['AlleleID'].map(lambda a: mondo_terms.get(a, []))
    matched['MONDO'] = matched['AlleleID'].map(lambda a: mondo_ids.get(a, []))

    matched['Phenotypes'] = matched['AlleleID'].map(lambda a: allele_phenotypes(a)[0])
    matched['Normal Phenotypes'] = matched['AlleleID'].map(lambda a: allele_phenotypes(a)[1])

    return matched


def _related_marker_rows(gene, marker, submitted_allele_ids):
    '''
    Coverage rows for the alleles of a cluster or overlapping gene, listed under the
    queried gene. VARIANT-ALLELE-MGI has no cluster alleles, so MGI's phenotypic
    allele report is added to it.
    '''
    rows = allele_coverage_for_gene(marker, submitted_allele_ids, all_alleles=True)
    listed = set(rows['AlleleID'])
    pa = _phenotypic_alleles()
    extra = pa[(pa['MarkerSymbol'] == marker) & ~pa['AlleleID'].isin(listed)]
    if not extra.empty:
        evaluated = set(submitted_allele_ids or [])
        disease = _load_disease_models()
        disease = disease[disease['AlleleID'].isin(extra['AlleleID'])].groupby('AlleleID')['Disease'] \
            .apply(lambda s: sorted(set(s.dropna()))).to_dict()
        added = pd.DataFrame({
            'Gene Symbol': marker,
            'AlleleID': extra['AlleleID'].values,
            'AlleleSymbol': [re.sub(r'<([^>]+)>', r'<sup>\1</sup>', s) for s in extra['AlleleSymbol']],
            'Consequence': '',
            'Diseases': [disease.get(a, []) for a in extra['AlleleID']],
            'Evaluated': [a in evaluated for a in extra['AlleleID']],
            'Reason': '',
        })
        added['Phenotypes'] = added['AlleleID'].map(lambda a: allele_phenotypes(a)[0])
        added['Normal Phenotypes'] = added['AlleleID'].map(lambda a: allele_phenotypes(a)[1])
        rows = pd.concat([rows, added], ignore_index=True)
    rows['Marker'] = marker
    rows['Gene Symbol'] = gene
    return rows


def query_allele_coverage(genes, submitted_allele_ids, positions=()):
    '''
    Coverage table for a mouse query: every MGI allele of each queried gene
    (disease models and alleles with variant data), flagging the alleles that
    match a submitted variant, with each allele's phenotypes.

    Alleles of MGI markers the variant also falls in (the gene's cluster, such as
    Pcdhg, and overlapping genes) are listed under the queried gene too, with their
    marker in 'Marker'. Many of them, such as large cluster deletions, have no
    variant coordinates in MGI, so the variant can't match them by position, but
    the user should still see that they exist.

    Parameters:
        genes: list of str. Queried mouse gene symbols.
        submitted_allele_ids: set of str. MGI alleles at the submitted coordinates.
        positions: iterable of (chromosome, position) of the submitted variants.
    '''
    positions = list(positions)
    frames = []
    for gene in genes:
        own = allele_coverage_for_gene(gene, submitted_allele_ids, all_alleles=True)
        own['Marker'] = gene
        frames.append(own)
        frames += [_related_marker_rows(gene, marker, submitted_allele_ids)
                   for marker in related_markers(gene, positions)]
    frames = [f for f in frames if not f.empty]
    coverage = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if coverage.empty:
        return coverage

    types = dict(zip(_phenotypic_alleles()['AlleleID'], _phenotypic_alleles()['Allele Type']))
    coverage['Allele Type'] = coverage['AlleleID'].map(types).fillna('')
    coverage = coverage.drop_duplicates(subset=['Gene Symbol', 'AlleleID'])

    # Reason explains why VarLift couldn't evaluate an allele, which doesn't apply here
    cols = ['Gene Symbol', 'Marker', 'AlleleID', 'AlleleSymbol', 'Allele Type', 'Consequence',
            'Diseases', 'Evaluated', 'Phenotypes', 'Normal Phenotypes']
    return coverage[[c for c in cols if c in coverage.columns]]
