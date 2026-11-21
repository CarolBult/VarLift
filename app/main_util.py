from app.log_util import log_batch_query
from app.processing_util import prepare_submission, process_human_fetch
from app.gene_util import fetch_gene_info, fetch_mus_alleles, fetch_homologous_gene, process_mus_alleles
from app.vep_util import prepare_vep_output, run_vep_scorable
from app.variant_scoring import structural_matches, total_score, add_best_tier, STRUCTURAL_COLS
from app.phenotype_util import allele_phenotype_match
from app.disease_util import assign_clinvar, fetch_mus_doid, clndisdb_to_mondo, map_doids_to_mondo, extract_hp_ids

import time
import yaml
import logging
import pandas as pd
from pathlib import Path

# Load config
with open('app/config.yaml') as f:
    config = yaml.safe_load(f)


class NoCodingVariantsError(Exception):
    '''Raised when VEP finds no gene consequences for any submitted variant.'''



def human_variant_diseases(variants):
    '''
    Collect ClinVar disease associations and HP terms for human variants.

    Parameters
        variants: DataFrame. Expects Name, Chromosome, Start, Stop, Ref, Alt.

    Returns
        DataFrame. Contains Name, Associated Diseases, MONDO and HP_IDs_variant.
    '''
    from app.clinvar_summary import lookup_hp_ids, lookup_disease_ids
    from app.hpoa import lookup_hps_for_mondo, current_hp_ids

    disease_df = assign_clinvar(variants)
    disease_df = extract_hp_ids(disease_df)        # HP IDs embedded directly in CLNDISDB
    disease_df = clndisdb_to_mondo(disease_df)     # resolves MONDO list per variant

    # The variant's own ClinVar records (exact ref/alt match in the VCF); variant_summary
    # is read for these records only, never for other changes at the same position
    allele_ids = dict(zip(disease_df['Name'], disease_df['ClinVar AlleleIDs']))

    # Variant-level HP terms: union of
    #   - direct HP cross-refs in CLNDISDB (extract_hp_ids result)
    #   - HPs from variant_summary.tsv's PhenotypeIDS column (richer)
    #   - HPs looked up by MONDO disease via phenotype.hpoa
    coord_to_vs_hp = {}
    for _, v in variants.iterrows():
        hps = lookup_hp_ids(v['Chromosome'], v['Start'], allele_ids.get(v['Name'], []))
        if hps:
            coord_to_vs_hp[v['Name']] = hps

    # Pre-fetch disease IDs from variant_summary.tsv for each input variant
    coord_to_vs_diseases = {}
    for _, v in variants.iterrows():
        dids = lookup_disease_ids(v['Chromosome'], v['Start'], allele_ids.get(v['Name'], []))
        if dids:
            coord_to_vs_diseases[v['Name']] = dids

    # Backfill Associated Diseases and MONDO from variant_summary.tsv when ClinVar's
    # VCF (CLNDISDB) didn't yield any MONDO disease associations: the same records'
    # variant_summary rows also list traits the VCF leaves out (e.g. somatic ones)
    from app.disease_util import mondo_xref_map, mondo_term_map

    def _backfill(row):
        if row.get('MONDO') and len(row['MONDO']) > 0:
            return row  # CLNDISDB already provided; keep it
        dids = coord_to_vs_diseases.get(row['Name'], [])
        if not dids:
            return row
        # OMIM and Orphanet IDs resolve to MONDO, so each disease is named once,
        # by its MONDO term, and the two lists stay in step as with CLNDISDB
        mondos = {d if d.startswith('MONDO:') else mondo_xref_map.get(d) for d in dids} - {None}
        pairs = sorted((mondo_term_map.get(m, m), m) for m in mondos)
        row['MONDO'] = [m for _, m in pairs]
        row['Associated Diseases'] = [name for name, _ in pairs]
        return row

    disease_df = disease_df.apply(_backfill, axis=1)
    disease_df = disease_df.drop(columns=['ClinVar AlleleIDs'])


    def _build_variant_hps(row):
        hps = set()
        # 1. Direct HP cross-refs embedded in CLNDISDB (via extract_hp_ids)
        clndisdb = row.get('HP_IDs', [])
        if isinstance(clndisdb, list):
            hps.update(clndisdb)
        # 2. Direct HPs from variant_summary.tsv
        hps.update(coord_to_vs_hp.get(row['Name'], []))
        # 3. HPs chained via MONDO disease IDs (from ClinVar's CLNDISDB)
        mondos = row.get('MONDO', [])
        if isinstance(mondos, list):
            for mondo in mondos:
                hps.update(lookup_hps_for_mondo(mondo))
        # 4. HPs chained via MONDO disease IDs (from variant_summary.tsv)
        for did in coord_to_vs_diseases.get(row['Name'], []):
            if did.startswith('MONDO:'):
                hps.update(lookup_hps_for_mondo(did))
            # OMIM and Orphanet feed hpoa directly, bypassing MONDO xref
            # (we'd need a parallel lookup in hpoa.py — do that if needed)
        # ClinVar cites some retired and obsolete HP IDs; count each phenotype once, by its current term
        return current_hp_ids(hps)

    disease_df['HP_IDs_variant'] = disease_df.apply(_build_variant_hps, axis=1)
    disease_df = disease_df.drop(columns=['HP_IDs'])

    return disease_df


# ---- HUMAN QUERY ----
def hvar_query(variants, assembly='GRCh38'):
    '''
    Annotate a batch of human variant using local gene metadata, Ensembl VEP, and ClinVar.

    Parameters

    Returns
        gene_df: pandas.DataFrame. Gene-level metadata and input variant summary.
        protein_df: pandas.DataFrame. Transcript-level information.
        input_gene_df: pandas.DataFrame. Gene mapping to input.
        unevaluated_df: pandas.DataFrame. Variants without a coding consequence.
    '''
    s = time.time()
    # Check assembly
    if assembly != 'GRCh38':
        logging.info(f'Assembly Error: {assembly} is not supported. Only GRCh38 is supported.')
        raise ValueError('Assembly must be GRCh38')

    # Log the query
    log_batch_query(variants)
    
    # Prepare common notation names
    variants, submission_names_map = prepare_submission(variants)


    # Query VEP for coding, splice-site and intronic consequences
    vep_df, noncoding_vep_df = run_vep_scorable(variants, 'homo_sapiens', query=True)

    # Variants that can't be scored are reported separately
    from app.unevaluated import unevaluated_variants
    unevaluated_df = unevaluated_variants(variants, vep_df.get('Submission', []), 'homo_sapiens', noncoding_vep_df)

    if vep_df.empty:
        if unevaluated_df.empty:
            raise NoCodingVariantsError('No protein-coding consequences found on GRCh38.')
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), unevaluated_df

    # Assign input field
    vep_df['Name'] = vep_df['Submission'].map(submission_names_map)

    # Clean VEP output into protein DataFrame
    protein_df = prepare_vep_output(vep_df)


    # Create a human gene to input mapping DataFrame
    input_gene_df = (protein_df[['Name', 'Gene Symbol']].drop_duplicates()
                     .rename(columns={'Gene Symbol': 'Hum Gene',
                                      'Name': 'Input Name'}).copy())

    # Build gene table
    genes = protein_df['Gene Symbol'].unique()
    gene_df = fetch_gene_info(genes, species='human')


    # Disease Associations
    from app.hpoa import lookup_hps_for_gene

    disease_df = human_variant_diseases(variants)

    protein_df = protein_df.merge(disease_df, on='Name', how='left')

    # Gene-level HP terms: looked up once per unique gene symbol in protein_df
    gene_to_hp = {g: lookup_hps_for_gene(g) for g in protein_df['Gene Symbol'].dropna().unique()}
    protein_df['HP_IDs_gene'] = protein_df['Gene Symbol'].map(gene_to_hp)

    # Save results
    Path("app/results").mkdir(parents=True, exist_ok=True)
    gene_df.to_csv('app/results/hum_gene_df.csv', index=False)
    protein_df.to_csv('app/results/hum_protein_df.csv', index=False)

    e = time.time()

    print(f'Human Time: {e-s}')

    return gene_df, protein_df, input_gene_df, unevaluated_df


def mvar_fetch(input_gene_df, assembly='GRCm39'):
    '''
    Annotate a mouse variant using local gene metadata, Ensembl VEP, and MouseMine.

    This function performs the following steps:
        1. Fetch mouse gene metadata.
        2. Identify all variants for the gene in MGI local database.
        3. Call Ensembl VEP for each variant → transcript-level annotations.
        4. Retrieve domain names from InterPro.
        5. Query MouseMine for ontology annotations.

    Parameters
        hum_gene_df: df. Human gene DF 
        assembly : string, optional. Genome assembly. Only "GRCm39" is supported.

    Returns
        gene_df: pandas.DataFrame. Gene-level metadata and input variant summary.
        protein_df: pandas.DataFrame. Transcript-level information.
        phenotype_df: pandas.DataFrame. Allele phenotype annotation.
        gene_input_df: pandas.DataFrame. Gene mapping to input and homolog.
    '''
    s = time.time()
    # Check assembly
    if assembly != 'GRCm39':
        raise ValueError('Assembly must be GRCm39')
    

    # Extract Mus gene symbol
    gene_input_df = fetch_homologous_gene(input_gene_df, 'mouse')

    # Build orthologous gene table
    genes = gene_input_df['Mus Gene'].unique()
    mouse_gene_df = fetch_gene_info(genes, species='mouse')


    # Extract gene ids
    MGI_gene_ids = mouse_gene_df['Accession'].unique()

    # Fetch mouse alleles
    mouse_allele_df = fetch_mus_alleles(MGI_gene_ids)
    mouse_prt_df = process_mus_alleles(mouse_allele_df)


    if mouse_prt_df.empty:
        logging.info('No mouse alleles found.')
        return mouse_gene_df, mouse_prt_df, pd.DataFrame(), gene_input_df, pd.DataFrame()

    # Extract phenotypes
    phenotype_df = allele_phenotype_match(mouse_allele_df)


    # Query MGD for ontology associations. Alleles without a disease model are
    # still scored: the human variant may match an allele no disease is annotated to
    doid_map = fetch_mus_doid(mouse_prt_df['AlleleID'].unique())

    # Prepare common notation names
    mouse_prt_df, _ = prepare_submission(mouse_prt_df)
    allele_map = mouse_prt_df[['Submission', 'Gene Symbol', 'AlleleID', 'AlleleSymbol']]


    # Query VEP for each variant
    variant_vep, _ = run_vep_scorable(mouse_prt_df, 'mus_musculus', query=False)

    if variant_vep.empty or 'consequence_terms' not in variant_vep.columns:
        logging.info('No mouse alleles with a coding, splice-site or intronic consequence.')
        return mouse_gene_df, pd.DataFrame(), phenotype_df, gene_input_df, pd.DataFrame()

    # Clean VEP output
    mouse_prt_df = prepare_vep_output(variant_vep)


    # Drop transcripts not in homologous gene
    mouse_prt_df = mouse_prt_df[mouse_prt_df['Gene Symbol'].isin(genes)]

    # Map Alleles back to VEP output
    allele_map = allele_map.rename(columns={'Gene Symbol': 'Gene Symbol Allele'})
    mouse_prt_df = mouse_prt_df.merge(allele_map, on='Submission', how='left')

    # Only keep Allele Gene Symbol
    mouse_prt_df.drop('Gene Symbol', inplace=True, axis=1)
    mouse_prt_df.rename(columns={'Gene Symbol Allele': 'Gene Symbol'}, inplace=True)


    # Compute MONDO mapping
    mondo_id_map, mondo_term_map = map_doids_to_mondo(doid_map) if doid_map else ({}, {})

    # Insert into the variant df (empty lists for alleles without a disease model)
    mouse_prt_df['MONDO'] = mouse_prt_df['AlleleID'].map(lambda a: mondo_id_map.get(a, []))
    mouse_prt_df['Associated Diseases'] = mouse_prt_df['AlleleID'].map(lambda a: mondo_term_map.get(a, []))


    # Rearrange columns
    mouse_prt_df = mouse_prt_df[['Submission', 'Gene Symbol', 'AlleleID', 'AlleleSymbol', 'Transcript ID', 'Biotype', 
                                 'Exon Rank', 'Intron Rank', 'Intron Offset', 'Pfam Domain ID', 'Molecular Consequence', 
                                 'Sequence Change', 'Amino Acids', 'refAA', 'varAA', 'Associated Diseases', 'MONDO']]


    # Replace all null with None
    mouse_gene_df = mouse_gene_df.where(pd.notnull(mouse_gene_df), None)
    mouse_prt_df = mouse_prt_df.where(pd.notnull(mouse_prt_df), None)
    phenotype_df = phenotype_df.where(pd.notnull(phenotype_df), None)

    # Save debugging tables to CSVs
    mouse_gene_df.to_csv('app/results/mouse_gene_df.csv', index=False)
    mouse_prt_df.to_csv('app/results/mouse_protein_df.csv', index=False)
    phenotype_df.to_csv('app/results/mouse_phenotype_df.csv', index=False)
    gene_input_df.to_csv('app/results/gene_input_df.csv', index=False)

    e = time.time()
    print(f'Mouse Time: {e-s}')

    # Allele coverage: every MGI allele of each output gene (disease models and
    # alleles with variant data), flagging which VarLift could evaluate
    from app.mouse_coverage import allele_coverage_for_gene
    evaluated_ids = set(mouse_prt_df['AlleleID'].unique())
    coverage_frames = []
    for gene in mouse_gene_df['Gene Symbol'].dropna().unique():
        coverage_frames.append(allele_coverage_for_gene(gene, evaluated_ids, all_alleles=True))
    allele_coverage_df = pd.concat(coverage_frames, ignore_index=True) if coverage_frames else pd.DataFrame()

    return mouse_gene_df, mouse_prt_df, phenotype_df, gene_input_df, allele_coverage_df

def hvar_query_score(hum_prt_df, mouse_prt_df, gene_inputs, mouse_phenotype_df=None):
    '''
    Creates score DF for human and mouse results.

    Parameters:
        hum_prt: pandas.DataFrame. Human variant dataframe.
        mouse_prt: pandas.DataFrame. Mouse model dataframe.
        gene_inputs: pandas.DataFrame. Gene-homolog-input mapping dataframe.

    Returns:
        pandas.DataFrame. contains score information for each variant-model pair.
    '''
    s = time.time()

    if mouse_prt_df.empty or hum_prt_df.empty:
        logging.info('No mouse alleles found.')
        return pd.DataFrame()


    hum_prt_df['MONDO_set'] = hum_prt_df['MONDO'].apply(set)
    mouse_prt_df['MONDO_set'] = mouse_prt_df['MONDO'].apply(set)

    expanded_df = gene_inputs.merge(hum_prt_df, left_on='Input Name', right_on='Name', suffixes=('', '_human'))
    expanded_df = expanded_df.merge(mouse_prt_df, left_on='Mus Gene', right_on='Gene Symbol', suffixes=('_human', '_mouse'))
    # None of the variants' homologs have mouse alleles (e.g. a variant VEP places only in an overlapping gene)
    if expanded_df.empty:
        logging.info('No mouse alleles for the variants\' homologs.')
        return pd.DataFrame()

    score_df = pd.DataFrame()
    score_df['Input Name'] = expanded_df['Input Name']
    score_df['AlleleID'] = expanded_df['AlleleID']
    score_df['Allele Symbol'] = expanded_df['AlleleSymbol']
    score_df['Transcript ID'] = expanded_df['Transcript ID_mouse']
    # The human transcript the mouse transcript was compared with (shown in the alignment view)
    score_df['Query Transcript ID'] = expanded_df['Transcript ID_human']
    score_df = pd.concat([score_df, structural_matches(expanded_df, '_human', '_mouse', 'human', 'mouse')], axis=1)
    # Mouse alleles without a disease model leave Disease Match out of the Total Score (None)
    score_df['Disease Match'] = expanded_df.apply(
        lambda r: bool(r['MONDO_set_human'] & r['MONDO_set_mouse']) if r['MONDO_set_mouse'] else None, axis=1)

    # columns which can be attributed to score
    match_cols = STRUCTURAL_COLS + ['Disease Match']

    # Percent of the criteria that apply to the variant's class. Superseded by the
    # Equivalency Tier; kept (hidden in the table) for the download and clinvar_batch
    score_df['Total Score'] = total_score(score_df, match_cols)
    score_df = add_best_tier(score_df, 'AlleleID')


    score_df = score_df.where(pd.notnull(score_df), None)

    e = time.time()
    print(f'Score Time: {e-s}')

    # Phenotype similarity, computed once per (variant, allele) and joined onto every
    # transcript row of that allele; see phenotype_similarity() below
    pheno_cols = ['Matched Disease(s)', 'Phenotype Score (variant)',
                  'Phenotype Score (gene)', 'MP Terms Mapped', 'MP Terms Total']
    pheno_df = phenotype_similarity(hum_prt_df, score_df[['Input Name', 'AlleleID']], mouse_phenotype_df)
    if pheno_df.empty:
        for c in pheno_cols:
            score_df[c] = None
    else:
        score_df = score_df.merge(pheno_df[['Input Name', 'AlleleID'] + pheno_cols],
                                  on=['Input Name', 'AlleleID'], how='left')

    # Belt-and-suspenders NaN handling: new columns from .apply() can get
    # promoted to float64 with NaN which doesn't JSON-serialize cleanly.
    for col in ['Matched Disease(s)', 'Phenotype Score (variant)',
                'Phenotype Score (gene)', 'MP Terms Mapped', 'MP Terms Total']:

        if col in score_df.columns:
            score_df[col] = score_df[col].astype(object).where(score_df[col].notna(), None)
    score_df = score_df.astype(object).where(pd.notnull(score_df), None)
    e = time.time()
    print(f'Score Time: {e-s}')

    score_df.to_csv('app/results/score_df.csv', index=False)

    return score_df

def phenotype_similarity(hum_prt_df, pairs, mouse_phenotype_df, detail=False):
    '''
    Phenotype similarity of mouse alleles to human variants, one row per (variant, allele).

    Two scores, both the asymmetric coverage in app.phenotype_similarity (the share of
    the allele's translated HP terms found in a human reference set):
    - Variant score: against each disease the variant is associated with, one at a time.
      The top score is reported with every disease that reaches it ("Matched Disease(s)");
      a tie has no single best disease.
    - Gene score: against the human gene's HP terms (genes_to_phenotype.txt).

    Parameters
        hum_prt_df: DataFrame. Human variants with Name, MONDO and HP_IDs_gene.
        pairs: DataFrame. Input Name and AlleleID of each (variant, allele) to score.
        mouse_phenotype_df: DataFrame. AlleleID, Phenotypes (abnormal) and optionally Normal
            Phenotypes, as lists of "MP:...,label" strings. Normal annotations are not scored.
        detail: bool. Also return the term mapping behind each score, for display:
            'Variant Mapping' and 'Gene Mapping' (see phenotype_mapping_detail).

    Returns
        DataFrame. Input Name, AlleleID, Matched Disease(s), Phenotype Score (variant),
        Phenotype Score (gene), MP Terms Mapped, MP Terms Total (and the mappings).
    '''
    from app.phenotype_similarity import phenotype_mapping_detail, scored_mp_terms
    from app.hpoa import lookup_hps_for_mondo, lookup_name_for_mondo, hp_names

    if (mouse_phenotype_df is None or mouse_phenotype_df.empty or pairs.empty
            or not ({'MONDO', 'HP_IDs_gene'} & set(hum_prt_df.columns))):
        return pd.DataFrame()

    # AlleleID -> set of MP IDs, and MP ID -> label
    allele_to_mp, mp_names = {}, {}
    for aid, phenos in mouse_phenotype_df[['AlleleID', 'Phenotypes']].itertuples(index=False):
        if isinstance(phenos, list):
            for p in phenos:
                if isinstance(p, str) and p.startswith('MP:'):
                    mp, _, label = p.partition(',')
                    allele_to_mp.setdefault(aid, set()).add(mp)
                    mp_names[mp] = label.strip()

    # Normal (NOT) annotations are kept out of the scores; the mapping view lists them
    allele_to_normal = {}
    if detail and 'Normal Phenotypes' in mouse_phenotype_df.columns:
        for aid, phenos in mouse_phenotype_df[['AlleleID', 'Normal Phenotypes']].itertuples(index=False):
            for p in phenos if isinstance(phenos, list) else []:
                mp, _, label = p.partition(',')
                allele_to_normal.setdefault(aid, {})[mp] = label.strip()

    variants = hum_prt_df.drop_duplicates(subset=['Name']).set_index('Name')
    def _list(name, col):
        v = variants[col].get(name) if col in variants.columns else None
        return list(v) if isinstance(v, (list, tuple)) else []

    disease_hps = {}
    def _disease_hps(m):
        if m not in disease_hps:
            disease_hps[m] = set(lookup_hps_for_mondo(m))
        return disease_hps[m]

    rows = []
    for name, aid in pairs.drop_duplicates().itertuples(index=False):
        # Top-level system terms and "no abnormal phenotype detected" are not scored
        mp_terms, unscored = scored_mp_terms(allele_to_mp.get(aid, set()))

        # Variant score: each of the variant's diseases with HPO annotations, separately
        per_disease = []
        for m in _list(name, 'MONDO'):
            if _disease_hps(m):
                d = phenotype_mapping_detail(mp_terms, _disease_hps(m), mp_names)
                if d['score'] is not None:
                    per_disease.append((d['score'], m, d))
        top = [x for x in per_disease if x[0] == max(s for s, _, _ in per_disease)] if per_disease else []
        named = {m: lookup_name_for_mondo(m) for _, m, _ in top}
        top_names = {m: n or m for m, n in named.items()}   # the mapping view falls back to the ID

        gene = phenotype_mapping_detail(mp_terms, set(_list(name, 'HP_IDs_gene')), mp_names)
        row = {'Input Name': name, 'AlleleID': aid,
               'Matched Disease(s)': '; '.join(sorted({n for n in named.values() if n})),
               'Phenotype Score (variant)': top[0][0] if top else None,
               'Phenotype Score (gene)': gene['score'],
               'MP Terms Mapped': gene['mapped_mp_count'], 'MP Terms Total': gene['total_mp_count']}

        if detail:
            # Each mapping lists one section per human reference: the variant's tied
            # top disease(s), each with its own arithmetic, or the human gene
            not_scored = [{'mp': mp, 'mp_name': mp_names.get(mp)} for mp in sorted(unscored)]
            normal = [{'mp': mp, 'mp_name': label}
                      for mp, label in sorted(allele_to_normal.get(aid, {}).items(), key=lambda x: x[1])]

            def _mapping(score, sections):
                base = sections[0][1] if sections else phenotype_mapping_detail(mp_terms, set(), mp_names)
                return {'score': score, 'total_mp_count': base['total_mp_count'],
                        'unmapped': base['unmapped'], 'not_scored': not_scored, 'normal': normal,
                        'sections': [dict(d, reference=ref) for ref, d in sections]}

            seen, variant_sections = set(), []
            for _, m, d in top:
                if top_names[m] not in seen:        # two MONDO IDs can share a name
                    seen.add(top_names[m])
                    variant_sections.append((top_names[m], d))
            variant = _mapping(top[0][0] if top else None, variant_sections)
            gene_map = _mapping(gene['score'], [('the human gene', gene)] if gene['score'] is not None else [])

            terms = {h for mp in (variant, gene_map) for sec in mp['sections'] for p in sec['pairs']
                     for h in (p['hp'], p['human_hp'], p['mica']) if h}
            names = hp_names(terms)
            for mp in (variant, gene_map):
                for sec in mp['sections']:
                    for p in sec['pairs']:
                        p['hp_name'] = names.get(p['hp'])
                        p['human_hp_name'] = names.get(p['human_hp'])
                        p['mica_name'] = names.get(p['mica'])
            row['Variant Mapping'] = variant
            row['Gene Mapping'] = gene_map
        rows.append(row)

    return pd.DataFrame(rows)


def gene_phenotype_similarity(allele_coverage_df, gene_inputs):
    '''
    Pooled gene-level phenotype similarity for a mouse query. A submitted mouse variant
    may be novel and match no MGI allele, so it has no phenotypes of its own. Instead the
    abnormal MP terms of every MGI allele of the queried mouse gene are pooled into one
    set and scored against the HPO terms of the human homolog (genes_to_phenotype.txt),
    giving one score per submitted variant.

    Parameters
        allele_coverage_df: DataFrame. Mouse Allele Coverage table (query_allele_coverage):
            Marker, AlleleID, Phenotypes and Normal Phenotypes. Alleles of related markers
            (gene clusters, overlapping genes) are not pooled.
        gene_inputs: DataFrame. Input Name, Mus Gene and Hum Gene.

    Returns
        DataFrame. One row per Input Name with the columns of phenotype_similarity(detail=True)
        plus Mus Gene, Hum Gene, Alleles Pooled, Alleles Total and Human HP Terms. AlleleID is the pooled
        set's label; Matched Disease(s) and Phenotype Score (variant) are empty.
    '''
    from app.hpoa import lookup_hps_for_gene

    if (allele_coverage_df is None or allele_coverage_df.empty or gene_inputs.empty
            or not {'Marker', 'AlleleID', 'Phenotypes'} <= set(allele_coverage_df.columns)
            or not {'Input Name', 'Mus Gene', 'Hum Gene'} <= set(gene_inputs.columns)):
        return pd.DataFrame()

    def _terms(values):
        return {p for v in values if isinstance(v, list) for p in v if isinstance(p, str) and p.startswith('MP:')}

    alleles = allele_coverage_df.drop_duplicates(subset=['Marker', 'AlleleID'])
    mapping = gene_inputs[['Input Name', 'Mus Gene', 'Hum Gene']].dropna().drop_duplicates()

    # One pooled phenotype set per mouse gene. A term abnormal in any allele counts as
    # abnormal; Normal lists only terms no allele of the gene has as abnormal
    pooled, counts = [], {}
    for gene in mapping['Mus Gene'].unique():
        own = alleles[alleles['Marker'] == gene]
        abnormal = _terms(own['Phenotypes'])
        if not abnormal:
            continue
        abnormal_ids = {p.partition(',')[0] for p in abnormal}
        normal = {p for p in _terms(own.get('Normal Phenotypes', pd.Series(dtype=object)))
                  if p.partition(',')[0] not in abnormal_ids}
        with_phenotypes = own['Phenotypes'].map(lambda p: isinstance(p, list) and len(p) > 0).sum()
        counts[gene] = (int(with_phenotypes), len(own))
        pooled.append({'AlleleID': f'{gene} (all alleles)', 'Phenotypes': sorted(abnormal),
                       'Normal Phenotypes': sorted(normal)})
    if not pooled:
        return pd.DataFrame()

    pairs = mapping[mapping['Mus Gene'].isin(counts)][['Input Name', 'Mus Gene']].drop_duplicates()
    pairs = pairs.assign(AlleleID=pairs['Mus Gene'] + ' (all alleles)')[['Input Name', 'AlleleID']]

    # Human reference per input: the HPO terms of its human homolog(s)
    gene_hps = {g: set(lookup_hps_for_gene(g) or []) for g in mapping['Hum Gene'].unique()}
    reference = (mapping.groupby('Input Name')['Hum Gene']
                 .agg(lambda genes: sorted(set().union(*(gene_hps[g] for g in genes))))
                 .rename('HP_IDs_gene').rename_axis('Name').reset_index())

    result = phenotype_similarity(reference, pairs, pd.DataFrame(pooled), detail=True)
    if result.empty:
        return result
    genes = mapping.groupby('Input Name').agg({'Mus Gene': ', '.join, 'Hum Gene': lambda g: ', '.join(sorted(set(g)))})
    result['Mus Gene'] = result['AlleleID'].str.replace(' (all alleles)', '', regex=False)
    result['Hum Gene'] = result['Input Name'].map(genes['Hum Gene'])
    result['Alleles Pooled'] = result['Mus Gene'].map(lambda g: counts[g][0])
    result['Alleles Total'] = result['Mus Gene'].map(lambda g: counts[g][1])
    result['Human HP Terms'] = result['Input Name'].map(reference.set_index('Name')['HP_IDs_gene'].map(len))
    return result

# ---- MOUSE QUERY ----
def mvar_query(variants, assembly='GRCm39'):
    '''
    '''
    s = time.time()
    # Check assembly
    if assembly != 'GRCm39':
        raise ValueError('Assembly must be GRCm39')
    
    # Log the query
    log_batch_query(variants)
    
    # Prepare HGVS notation
    variants, submission_name_map = prepare_submission(variants)


    # Query VEP for coding, splice-site and intronic consequences
    vep_df, noncoding_vep_df = run_vep_scorable(variants, 'mus_musculus', query=True)

    # Variants that can't be scored are reported separately
    from app.unevaluated import unevaluated_variants
    unevaluated_df = unevaluated_variants(variants, vep_df.get('Submission', []), 'mus_musculus', noncoding_vep_df)

    # MGI alleles at the submitted coordinates
    from app.mouse_coverage import mgi_alleles_at, query_allele_coverage
    alleles = mgi_alleles_at(variants)
    submitted_ids = set(alleles['AlleleID'])
    # Submitted positions, for the clusters and overlapping genes listed in the coverage table
    positions = list(zip(variants['Chromosome'].astype(str), variants['Start'].astype(int)))
    unevaluated_genes = set(unevaluated_df.get('Gene Symbol', pd.Series(dtype=object)).dropna())

    if vep_df.empty:
        if unevaluated_df.empty:
            raise NoCodingVariantsError('No protein-coding consequences found on GRCm39.')
        allele_coverage_df = query_allele_coverage(sorted(unevaluated_genes), submitted_ids, positions)
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), unevaluated_df, allele_coverage_df

    # Assign input field
    vep_df['Name'] = vep_df['Submission'].map(submission_name_map)


    # Clean VEP output into protein DataFrame
    protein_df = prepare_vep_output(vep_df)


    # Add the matching MGI alleles and their disease models to each variant
    allele_cols = (alleles.groupby('Input Name')
                   .agg(**{'AlleleID': ('AlleleID', ', '.join),
                           'AlleleSymbol': ('AlleleSymbol', ', '.join)})
                   .reset_index().rename(columns={'Input Name': 'Name'}))

    # Disease models of all the matched alleles, as (term, MONDO ID) pairs kept in step
    disease_pairs = alleles.groupby('Input Name').apply(
        lambda g: sorted({p for ds, ms in zip(g['Associated Diseases'], g['MONDO']) for p in zip(ds, ms)}))
    allele_cols['Associated Diseases'] = allele_cols['Name'].map(lambda n: [d for d, _ in disease_pairs.get(n, [])])
    allele_cols['MONDO'] = allele_cols['Name'].map(lambda n: [m for _, m in disease_pairs.get(n, [])])
    protein_df = protein_df.merge(allele_cols, on='Name', how='left')


    # Create a gene to input mapping DataFrame
    input_gene_df = (protein_df[['Name', 'Gene Symbol']].drop_duplicates()
                     .rename(columns={'Gene Symbol': 'Mus Gene',
                                      'Name': 'Input Name'}).copy())

    # Build gene table
    genes = protein_df['Gene Symbol'].unique()
    gene_df = fetch_gene_info(genes, species='mouse')

    # Replace all null with None
    gene_df = gene_df.where(pd.notnull(gene_df), None)
    protein_df = protein_df.where(pd.notnull(protein_df), None)


    # Save results
    Path("app/results").mkdir(parents=True, exist_ok=True)
    gene_df.to_csv('app/results/mouse_gene_df.csv', index=False)
    protein_df.to_csv('app/results/mouse_protein_df.csv', index=False)

    e = time.time()
    print(f'Mouse Time: {e-s}')

    # Every MGI allele of the queried genes, flagging the submitted ones
    genes = sorted(set(protein_df['Gene Symbol'].dropna()) | unevaluated_genes)
    allele_coverage_df = query_allele_coverage(genes, submitted_ids, positions)

    return gene_df, protein_df, input_gene_df, unevaluated_df, allele_coverage_df


def hvar_fetch(mouse_prt_df, input_mapping_df, assembly='GRCh38'):
    '''
    '''
    s = time.time()

    # Check assembly
    if assembly != 'GRCh38':
        logging.info(f'Assembly Error: {assembly} is not supported. Only GRCh38 is supported.')
        raise ValueError('Assembly must be GRCh38')
        
    
    # Extract Hum gene symbol
    gene_input_df = fetch_homologous_gene(input_mapping_df, 'human')


    # Build orthologous gene table
    genes = gene_input_df['Hum Gene'].unique()
    hum_gene_df = fetch_gene_info(genes, species='human')

    # Extract NCBI variants
    homologous_variants, phenotype_df = process_human_fetch(hum_gene_df, mouse_prt_df, gene_input_df)

    if homologous_variants.empty:
        logging.info('No homologous variants found.')
        return hum_gene_df, pd.DataFrame(), pd.DataFrame(), gene_input_df


    # Query VEP for each variant
    variants, submission_HGVS_map = prepare_submission(homologous_variants)
    
    variant_vep, _ = run_vep_scorable(variants, 'homo_sapiens', query=False)


    # Assign name field
    variant_vep['Name'] = variant_vep['Submission'].map(submission_HGVS_map)

    # Clean VEP output
    hum_protein_df = prepare_vep_output(variant_vep)


    # Disease Associations
    disease_df = assign_clinvar(variants)
    disease_df = clndisdb_to_mondo(disease_df).drop(columns=['ClinVar AlleleIDs'])
    hum_protein_df = hum_protein_df.merge(disease_df, on='Name', how='left')

    hum_gene_df = hum_gene_df.where(pd.notnull(hum_gene_df), None)
    hum_protein_df = hum_protein_df.where(pd.notnull(hum_protein_df), None)

    hum_gene_df.to_csv('app/results/hum_gene_df.csv', index=False)
    hum_protein_df.to_csv('app/results/hum_protein_df.csv', index=False)
    phenotype_df.to_csv('app/results/hum_phenotype_df.csv', index=False)
    gene_input_df.to_csv('app/results/gene_input_df.csv', index=False)

    e = time.time()

    print(f'Human Time: {e-s}')


    return hum_gene_df, hum_protein_df, phenotype_df, gene_input_df


def mvar_query_score(mouse_prt_df, hum_prt_df, gene_inputs):
    '''
    Creates score df for mouse and human results.

    Parameters:
        mouse_prt: pandas.DataFrame. Mouse variant dataframe.
        hum_prt: pandas.DataFrame. Human variant dataframe.
        gene_inputs: pandas.DataFrame. Gene-homolog-input mapping dataframe.

    Returns:
        pandas.DataFrame. contains score information for each variant-model pair.
    '''
    s = time.time()

    if mouse_prt_df.empty or hum_prt_df.empty:
        logging.info('No homologous variants found.')
        return pd.DataFrame()


    expanded_df = gene_inputs.merge(hum_prt_df, left_on='Hum Gene', right_on='Gene Symbol', suffixes=('', '_human'))
    expanded_df = expanded_df.merge(mouse_prt_df, left_on='Input Name', right_on='Name', suffixes=('', '_mouse'))

    score_df = pd.DataFrame()
    score_df['Input Name'] = expanded_df['Input Name']
    score_df['Name'] = expanded_df['Name']
    score_df['Transcript ID'] = expanded_df['Transcript ID']
    # The mouse transcript the human transcript was compared with (shown in the alignment view)
    score_df['Query Transcript ID'] = expanded_df['Transcript ID_mouse']
    score_df = pd.concat([score_df, structural_matches(expanded_df, '_mouse', '', 'mouse', 'human')], axis=1)

    # columns which can be attributed to score
    match_cols = STRUCTURAL_COLS

    # Percent of the criteria that apply to the variant's class. Superseded by the
    # Equivalency Tier; kept (hidden in the table) for the download
    score_df['Total Score'] = total_score(score_df, match_cols)
    score_df = add_best_tier(score_df, 'Name')


    score_df = score_df.astype(object).where(pd.notnull(score_df), None)

    e = time.time()
    print(f'Score Time: {e-s}')

    score_df.to_csv('app/results/score_df.csv', index=False)

    return score_df
