from app.log_util import log_batch_query
from app.processing_util import prepare_submission, process_human_fetch
from app.gene_util import fetch_gene_info, fetch_mus_alleles, fetch_homologous_gene, process_mus_alleles
from app.vep_util import prepare_vep_output, run_vep
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


# Load homology
homology_df = pd.read_csv(config['paths']['mgi_homology'],
                          sep='\t', names=['MGI_ID', 'MusGeneSymbol', 'MusEntrezGeneID', 
                                           'MusHGNC_ID', 'HumGeneSymbol', 'HumEntrezGeneID'])
homology_dict = dict(zip(homology_df['HumGeneSymbol'], homology_df['MusGeneSymbol']))



# ---- HUMAN QUERY ----
def hvar_query(variants, assembly='GRCh38'):
    '''
    Annotate a batch of human variant using local gene metadata, Ensembl VEP, and ClinVar.

    Parameters

    Returns
        gene_df: pandas.DataFrame. Gene-level metadata and input variant summary.
        protein_df: pandas.DataFrame. Transcript-level information.
        input_gene_df: pandas.DataFrame. Gene mapping to input.
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


    # Query VEP
    vep_df = run_vep(variants, 'homo_sapiens', query=True)


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
    from app.clinvar_summary import lookup_hp_ids, lookup_disease_ids
    from app.hpoa import lookup_hps_for_mondo, lookup_hps_for_gene

    disease_df = assign_clinvar(variants)
    disease_df = extract_hp_ids(disease_df)        # HP IDs embedded directly in CLNDISDB
    disease_df = clndisdb_to_mondo(disease_df)     # resolves MONDO list per variant

    # Variant-level HP terms: union of
    #   - direct HP cross-refs in CLNDISDB (extract_hp_ids result)
    #   - HPs from variant_summary.tsv's PhenotypeIDS column (richer, by chrom/pos)
    #   - HPs looked up by MONDO disease via phenotype.hpoa
    coord_to_vs_hp = {}
    for _, v in variants.iterrows():
        hps = lookup_hp_ids(v['Chromosome'], v['Start'])
        if hps:
            coord_to_vs_hp[v['Name']] = hps

    # Pre-fetch disease IDs from variant_summary.tsv for each input variant
    coord_to_vs_diseases = {}
    for _, v in variants.iterrows():
        dids = lookup_disease_ids(v['Chromosome'], v['Start'])
        if dids:
            coord_to_vs_diseases[v['Name']] = dids

    # Backfill Associated Diseases and MONDO from variant_summary.tsv when
    # ClinVar's VCF (CLNDISDB) didn't yield any disease associations. Many
    # protein-level ClinVar records are absent from the VCF but still carry
    # disease cross-references in variant_summary.tsv.
    from app.hpoa import lookup_name_for_disease

    name_to_diseases = {}
    for _, v in variants.iterrows():
        dids = lookup_disease_ids(v['Chromosome'], v['Start'])
        if dids:
            name_to_diseases[v['Name']] = dids

    def _backfill(row):
        if row.get('MONDO') and len(row['MONDO']) > 0:
            return row  # CLNDISDB already provided; keep it
        dids = name_to_diseases.get(row['Name'], [])
        if not dids:
            return row
        mondos = sorted({d for d in dids if d.startswith('MONDO:')})
        # Build names: prefer MONDO names via OMIM/Orphanet from hpoa
        names = []
        seen_names = set()
        for d in dids:
            n = lookup_name_for_disease(d)
            if n and n not in seen_names:
                names.append(n)
                seen_names.add(n)
        row['MONDO'] = mondos
        row['Associated Diseases'] = names
        return row

    disease_df = disease_df.apply(_backfill, axis=1)




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
        return sorted(hps)

    disease_df['HP_IDs_variant'] = disease_df.apply(_build_variant_hps, axis=1)
    disease_df = disease_df.drop(columns=['HP_IDs'])

    protein_df = protein_df.merge(disease_df, on='Name', how='left')

    # Gene-level HP terms: looked up once per unique gene symbol in protein_df
    gene_to_hp = {g: lookup_hps_for_gene(g) for g in protein_df['Gene Symbol'].dropna().unique()}
    protein_df['HP_IDs_gene'] = protein_df['Gene Symbol'].map(gene_to_hp)

    # Save results
    Path("app/results").mkdir(parents=True, exist_ok=True)
    gene_df.to_csv('app/results/hum_gene_df.csv', index=False)
    protein_df.to_csv('app/results/hum_protein_df.csv', index=False)

    # Print resulting tables
    # print(gene_df)
    # print('----------')
    # print(protein_df)
    # print('----------')
    # print(input_gene_df)
    # print('----------')
    e = time.time()

    print(f'Human Time: {e-s}')

    return gene_df, protein_df, input_gene_df


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
        return mouse_gene_df, mouse_prt_df, pd.DataFrame(), gene_input_df


    # Extract phenotypes
    phenotype_df = allele_phenotype_match(mouse_allele_df)


    # Query MGD for ontology associations
    try:
        doid_map = fetch_mus_doid(mouse_prt_df['AlleleID'].unique())
    except:
        logging.info('No Alleles found with disease associations.')
        return mouse_gene_df, mouse_prt_df, phenotype_df, gene_input_df

    # Drop alleles without disease associations
    mouse_prt_df = mouse_prt_df[mouse_prt_df['AlleleID'].isin(doid_map.keys())]


    if mouse_prt_df.empty:
        logging.info('No disease associated mouse alleles found.')
        return mouse_gene_df, mouse_prt_df, pd.DataFrame(), gene_input_df


    # Prepare common notation names
    mouse_prt_df, _ = prepare_submission(mouse_prt_df)
    allele_map = mouse_prt_df[['Submission', 'Gene Symbol', 'AlleleID', 'AlleleSymbol']]


    # Query VEP for each variant
    variant_vep = run_vep(mouse_prt_df, 'mus_musculus', query=False)

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
    mondo_id_map, mondo_term_map = map_doids_to_mondo(doid_map)

    # Insert into the variant df
    mouse_prt_df['MONDO'] = mouse_prt_df['AlleleID'].map(mondo_id_map)
    mouse_prt_df['Associated Diseases'] = mouse_prt_df['AlleleID'].map(mondo_term_map)


    # Rearrange columns
    mouse_prt_df = mouse_prt_df[['Submission', 'Gene Symbol', 'AlleleID', 'AlleleSymbol', 'Transcript ID', 'Biotype', 
                                 'Exon Rank', 'Pfam Domain ID', 'Molecular Consequence', 
                                 'Codon Switch', 'Amino Acids', 'refAA', 'varAA', 'Associated Diseases', 'MONDO']]


    # Replace all null with None
    mouse_gene_df = mouse_gene_df.where(pd.notnull(mouse_gene_df), None)
    mouse_prt_df = mouse_prt_df.where(pd.notnull(mouse_prt_df), None)
    phenotype_df = phenotype_df.where(pd.notnull(phenotype_df), None)

    # Save debugging tables to CSVs
    mouse_gene_df.to_csv('app/results/mouse_gene_df.csv', index=False)
    mouse_prt_df.to_csv('app/results/mouse_protein_df.csv', index=False)
    phenotype_df.to_csv('app/results/mouse_phenotype_df.csv', index=False)
    gene_input_df.to_csv('app/results/gene_input_df.csv', index=False)

    # Print resulting tables
    # print(mouse_gene_df)
    # print('----------')
    # print(mouse_prt_df)
    # print('----------')

    e = time.time()
    print(f'Mouse Time: {e-s}')

    return mouse_gene_df, mouse_prt_df, phenotype_df, gene_input_df

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

    score_df = pd.DataFrame()
    score_df['Input Name'] = expanded_df['Input Name']
    score_df['AlleleID'] = expanded_df['AlleleID']
    score_df['Allele Symbol'] = expanded_df['AlleleSymbol']
    score_df['Transcript ID'] = expanded_df['Transcript ID_mouse']
    score_df['Biotype Match'] = expanded_df['Biotype_human'] == expanded_df['Biotype_mouse']
    score_df['Consequence Match'] = expanded_df['Molecular Consequence_human'] == expanded_df['Molecular Consequence_mouse']
    score_df['AA Match'] = (expanded_df['refAA_human'] == expanded_df['refAA_mouse']) & (expanded_df['varAA_human'] == expanded_df['varAA_mouse'])
    score_df['AA Position Match'] = expanded_df['Amino Acids_human'] == expanded_df['Amino Acids_mouse']
    score_df['Exon Match'] = expanded_df['Exon Rank_human'] == expanded_df['Exon Rank_mouse']
    score_df['Domain Match'] = expanded_df['Pfam Domain ID_human'] == expanded_df['Pfam Domain ID_mouse']
    score_df['Disease Match'] = expanded_df.apply(lambda r: bool(r['MONDO_set_human'] & r['MONDO_set_mouse']), axis=1)

    # columns which can be attributed to score
    match_cols = ['Biotype Match', 'Consequence Match', 'AA Match', 'AA Position Match', 'Exon Match', 'Domain Match', 'Disease Match']

    # Calculate precentage of hits
    score_df['Total Score'] = (score_df[match_cols].sum(axis=1) / score_df[match_cols].notna().sum(axis=1) * 100).round(1)


    score_df = score_df.where(pd.notnull(score_df), None)

    e = time.time()
    print(f'Score Time: {e-s}')

    # Calculate precentage of hits
    score_df['Total Score'] = (score_df[match_cols].sum(axis=1) / score_df[match_cols].notna().sum(axis=1) * 100).round(1)

    # Phenotype equivalency scoring (two tracks):
    # - Variant score: asymmetric coverage of the mouse allele's translated
    #   HPs against each disease the input variant is associated with. We
    #   report the top-scoring disease(s) and the shared top score — so a
    #   researcher sees which disease the match is for, not an aggregate
    #   across diseases.
    # - Gene score: asymmetric coverage against the gene's full HP landscape
    #   (independent of disease choice).
    from app.phenotype_similarity import phenotype_coverage_detailed
    from app.hpoa import lookup_hps_for_mondo, lookup_name_for_mondo

    cols = ['Matched Disease(s)', 'Phenotype Score (variant)',
            'Phenotype Score (gene)', 'MP Terms Mapped', 'MP Terms Total']

    have_mouse = mouse_phenotype_df is not None and not mouse_phenotype_df.empty
    have_mondo = 'MONDO' in hum_prt_df.columns
    have_hp_gene = 'HP_IDs_gene' in hum_prt_df.columns

    if have_mouse and (have_mondo or have_hp_gene):
        # AlleleID -> set of MP IDs (strip ",label" suffix)
        allele_to_mp = {}
        for aid, phenos in mouse_phenotype_df[['AlleleID', 'Phenotypes']].itertuples(index=False):
            if isinstance(phenos, list):
                mps = [p.split(',')[0] for p in phenos if isinstance(p, str) and p.startswith('MP:')]
                allele_to_mp.setdefault(aid, set()).update(mps)

        # Input Name -> list of MONDO disease IDs for the variant
        name_to_mondos = {}
        if have_mondo:
            for name, mondos in hum_prt_df[['Name', 'MONDO']].drop_duplicates(subset=['Name']).itertuples(index=False):
                name_to_mondos[name] = list(mondos) if isinstance(mondos, list) else []

        # Input Name -> gene HP set (invariant across the variant's transcripts)
        name_to_hp_gene = {}
        if have_hp_gene:
            for name, hps in hum_prt_df[['Name', 'HP_IDs_gene']].drop_duplicates(subset=['Name']).itertuples(index=False):
                name_to_hp_gene[name] = set(hps) if isinstance(hps, list) else set()

        def _score_row(r):
            mp_terms = allele_to_mp.get(r['AlleleID'], set())
            gene_hps = name_to_hp_gene.get(r['Input Name'], set())

            # Per-disease variant scoring
            mondos = name_to_mondos.get(r['Input Name'], [])
            per_disease = []  # list of (score, mondo_id)
            for m in mondos:
                disease_hps = set(lookup_hps_for_mondo(m))
                if not disease_hps:
                    continue
                d = phenotype_coverage_detailed(mp_terms, disease_hps)
                if d['score'] is not None:
                    per_disease.append((d['score'], m))

            if per_disease:
                top_score = max(s for s, _ in per_disease)
                top_mondos = [m for s, m in per_disease if s == top_score]
                names = sorted({n for n in (lookup_name_for_mondo(m) for m in top_mondos) if n})
                matched_names = '; '.join(names)
                variant_score = top_score
            else:
                matched_names = ''
                variant_score = None

            # Gene score (asymmetric). Its mapped/total counts describe the
            # mouse MP set and are the same regardless of which human reference
            # we score against, so we use them for the two counter columns.
            d_gene = phenotype_coverage_detailed(mp_terms, gene_hps)

            return pd.Series(
                [matched_names, variant_score, d_gene['score'],
                 d_gene['mapped_mp_count'], d_gene['total_mp_count']],
                index=cols,
            )

        score_df[cols] = score_df.apply(_score_row, axis=1)
    else:
        for c in cols:
            score_df[c] = None

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


    # Query VEP
    vep_df = run_vep(variants, 'mus_musculus', query=True)


    # Assign input field
    vep_df['Name'] = vep_df['Submission'].map(submission_name_map)


    # Clean VEP output into protein DataFrame
    protein_df = prepare_vep_output(vep_df)


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

    # # Print resulting tables
    # print(gene_df)
    # print('----------')
    # print(protein_df)
    # print('----------')

    e = time.time()
    print(f'Mouse Time: {e-s}')

    return gene_df, protein_df, input_gene_df


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
    
    variant_vep = run_vep(variants, 'homo_sapiens', query=False)


    # Assign name field
    variant_vep['Name'] = variant_vep['Submission'].map(submission_HGVS_map)

    # Clean VEP output
    hum_protein_df = prepare_vep_output(variant_vep)


    # Disease Associations
    disease_df = assign_clinvar(variants)
    disease_df = clndisdb_to_mondo(disease_df)
    hum_protein_df = hum_protein_df.merge(disease_df, on='Name', how='left')

    hum_gene_df = hum_gene_df.where(pd.notnull(hum_gene_df), None)
    hum_protein_df = hum_protein_df.where(pd.notnull(hum_protein_df), None)

    hum_gene_df.to_csv('app/results/hum_gene_df.csv', index=False)
    hum_protein_df.to_csv('app/results/hum_protein_df.csv', index=False)
    phenotype_df.to_csv('app/results/hum_phenotype_df.csv', index=False)
    gene_input_df.to_csv('app/results/gene_input_df.csv', index=False)

    # # Print resulting tables
    # print(hum_gene_df)
    # print('----------')
    # print(hum_protein_df)
    # print('----------')


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

    expanded_df.to_csv('app/testing_results/score_test.csv', index=False)

    score_df = pd.DataFrame()
    score_df['Input Name'] = expanded_df['Input Name']
    score_df['Name'] = expanded_df['Name']
    score_df['Transcript ID'] = expanded_df['Transcript ID']
    score_df['Biotype Match'] = expanded_df['Biotype_mouse'] == expanded_df['Biotype']
    score_df['Consequence Match'] = expanded_df['Molecular Consequence_mouse'] == expanded_df['Molecular Consequence']
    score_df['AA Match'] = (expanded_df['refAA_mouse'] == expanded_df['refAA']) & (expanded_df['varAA_mouse'] == expanded_df['varAA'])
    score_df['AA Position Match'] = expanded_df['Amino Acids_mouse'] == expanded_df['Amino Acids']
    score_df['Exon Match'] = expanded_df['Exon Rank_mouse'] == expanded_df['Exon Rank']
    score_df['Domain Match'] = expanded_df['Pfam Domain ID_mouse'] == expanded_df['Pfam Domain ID']

    # columns which can be attributed to score
    match_cols = ['Biotype Match', 'Consequence Match', 'AA Match', 'AA Position Match', 'Exon Match', 'Domain Match']

    # Calculate precentage of hits
    score_df['Total Score'] = (score_df[match_cols].sum(axis=1) / score_df[match_cols].notna().sum(axis=1) * 100).round(1)


    score_df = score_df.astype(object).where(pd.notnull(score_df), None)

    e = time.time()
    print(f'Score Time: {e-s}')

    score_df.to_csv('app/results/score_df.csv', index=False)

    return score_df
