import yaml
import json
import logging
import pandas as pd
from cyvcf2 import VCF


# Load Config
with open('app/config.yaml') as f:
    config = yaml.safe_load(f)


clinvar_vcf = VCF(config['paths']['clinvar_vcf'])


MGI_disease_df = pd.read_csv(config['paths']['mgi_disease'],
                             sep='\t', comment='#', header=None,
                             names=['DOterm', 'DOID', 'NOTmodel', 'AllelePairs', 'StrainBackground',
                                    'AlleleSymbol', 'AlleleID', 'NumReferences', 'AlleleRepositoryID',
                                    'AlleleRRID', 'MarkerSymbol', 'MarkerMGIid', 'GeneRepositoryID', ''])


mondo_xref_map = json.load(open(config['paths']['mondo_xref_map']))
mondo_term_map = json.load(open(config['paths']['mondo_term_map']))



def assign_clinvar(variants):
    '''
    Fetch ClinVar disease associations for a list of variants.

    Parameters:
        variants: DataFrame. Expects Name, Chromosome, Start, Stop, Ref, Alt.
    
    Returns:
        DataFrame. Contains Name and CLINDISDB.
    '''
    results = []

    # Precompute HGVS list
    hgvs_list = variants['Name'].tolist()

    # Cache to avoid redundant queries
    cache = {}

    for idx, row in variants.iterrows():
        chrom, start, stop, ref, alt = row.Chromosome, row.Start, row.Stop, row.Ref, row.Alt
        key = (chrom, start, stop, ref, alt)

        if key not in cache:
            cache[key] = search_clinvar(chrom, start, stop, ref, alt)

        results.append({'Name': hgvs_list[idx], 'CLNDISDB': cache[key]})

    return pd.DataFrame(results)

def search_clinvar(chrom, start, end, ref, alt):
    '''
    Searches ClinVar VCF for disease associations with a given variant.

    Parameters:
        chrom: str. Chromosome of variant.
        start: int. Genomic location of variant start.
        end: int. Genomic location of variant end.
        ref: string. Reference allele.
        alt: string. Alternate allele.
    
    Returns:
        list. containing: chromosome, start, end, ref, alt, diseases
    '''
    # Craft region string
    region = f'{chrom}:{start}-{end}'
    diseases = None

    try:
        # Filter by region and search for variants
        for record in clinvar_vcf(region):
            if record.REF == ref and alt in record.ALT:

                # Extract disease names
                disease_info = record.INFO.get('CLNDISDB')

                if disease_info: 
                    diseases = [d.split(',') for d in disease_info.split('|')]


    except Exception as e:
        logging.warning(f'ClinVar search failed for {region}: {e}')
            
    return diseases


def fetch_mus_doid(allele_ids):
    '''
    Fetches all Mouse disease associations with human models for all alleles from list of MGI_gene_ids.
    
    Parameters:
        MGI_gene_ids: list of str. A list of all MGI_gene_ids to return strings for.
    
    Returns:
        dictionary. Mapping of AlleleID to set of DOIDs.
    '''

    filtered_disease_df = MGI_disease_df[MGI_disease_df['AlleleID'].isin(allele_ids)]

    doid_allele_dict  = (filtered_disease_df.groupby('AlleleID')['DOID'].apply(set).to_dict())
        
    return doid_allele_dict

def map_doids_to_mondo(doid_map):
    '''
    Convert {AlleleID: {DOID}} to a long table, merge MONDO info, then collapse back into per-allele lists.

    Parameters:
        doid_map: dict. Mapping of AlleleID to DOID.

    Returns:
        tuple of dicts. Mapping of AlleleID to MONDO ID and mapping of AlleleID to Disease Association terms.
    '''
    # Normalize to long df
    rows = [{'AlleleID': allele, 'DOID': doid}
            for allele, doids in doid_map.items()
            for doid in doids]
    
    df = pd.DataFrame(rows)

    # Map MONDO ID and MONDO term
    df['MONDO'] = df['DOID'].map(mondo_xref_map)
    df['Disease Association'] = df['MONDO'].map(mondo_term_map)

    # Collapse back to allele level
    mondo_by_allele = df.groupby('AlleleID')['MONDO'].apply(lambda s: sorted(set(s.dropna())))
    term_by_allele  = df.groupby('AlleleID')['Disease Association'].apply(lambda s: sorted(set(s.dropna())))

    return mondo_by_allele.to_dict(), term_by_allele.to_dict()


def clndisdb_to_mondo(clinvar_results):
    '''
    Extracts MONDO IDs and terms from ClinVar CLNDISDB entries.

    Parameters:
        clinvar_results: DataFrame. containing 'CLNDISDB' column with disease entries.
    
    Returns:
        DataFrame. with added 'MONDO' and 'Associated Diseases' columns.
    '''
    mondos_col = []
    terms_col = []

    # Iterate rows
    for _, disease_list in enumerate(clinvar_results['CLNDISDB'].values):

        mondos = []
        terms = []

        # Variants not in ClinVar come through as None; preserve column alignment
        if disease_list:
            for disease_tokens in disease_list:

                # Extract MONDO and term
                mondo, term = resolve_disease_tokens(disease_tokens)

                if mondo:
                    mondos.append(mondo)
                    terms.append(term)

        mondos_col.append(mondos)
        terms_col.append(terms)
    # Assign results to HGVS DataFrame
    clinvar_results['Associated Diseases'] = terms_col
    clinvar_results['MONDO'] = mondos_col

    clinvar_results.drop('CLNDISDB', inplace=True, axis=1)

    return clinvar_results

def resolve_disease_tokens(tokens):
    '''
    Given a list of tokens from a CLNDISDB disease entry, extracts MONDO term.

    Parameters:
        tokens: list of str. Tokens from a CLNDISDB entry (disease IDs).

    Returns:
        tuple: MONDO ID (str) and MONDO term (str). If none found, (None, None).
    '''
    for token in tokens:

        # Direct MONDO entry
        if token.startswith("MONDO:"):
            _, mondo = token.split(":", 1)
            return mondo, mondo_term_map.get(mondo)

        # Normalize MedGen to UMLS
        if token.startswith("MedGen:"):
            token = "UMLS:" + token.split(":")[1]

        # Xref lookup
        if token in mondo_xref_map:
            mondo = mondo_xref_map[token]
            return mondo, mondo_term_map.get(mondo)

    # No valid MONDO mapping
    return None, None

def extract_hp_ids(clinvar_results):
    '''
    Extract HP (Human Phenotype Ontology) IDs from the CLNDISDB field.

    CLNDISDB entries look like:
      [['MedGen:C...', 'OMIM:...', 'Human_Phenotype_Ontology:HP:0001234'], ...]

    Adds a 'HP_IDs' column (list of 'HP:...' strings per variant).
    '''
    hp_col = []
    for disease_list in clinvar_results['CLNDISDB'].values:
        hp_ids = set()
        if disease_list:
            for disease_entry in disease_list:
                for xref in disease_entry:
                    if xref.startswith('Human_Phenotype_Ontology:HP:'):
                        hp_ids.add(xref.split(':', 1)[1])  # keep "HP:0001234"
        hp_col.append(sorted(hp_ids))
    clinvar_results['HP_IDs'] = hp_col
    return clinvar_results
