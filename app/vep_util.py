import yaml
import logging
import requests
import numpy as np
import pandas as pd
import subprocess
import json
import os
from pathlib import Path


# Load Config
with open('app/config.yaml') as f:
    config = yaml.safe_load(f)

# Load APIs
ensembl_base_url = config['api']['ensembl_base_url']
interpro_base_url = config['api']['interpro_base_url']

# Load VEP paths
apptainer_main = config['paths']['apptainer_main']
vep_output_file = config['paths']['vep_output']
vep_input_file = config['paths']['vep_input']


def get_domain_name(pfam_domain_id):
    '''
    Query InterPro Pfam database to get domain name.

    Parameters:
        pfam_domain_id: string. Pfam domain ID to extract full name for.

    Returns:
        string. The corresponding long name for given domain ID.
    '''
    # Construct URL
    url = interpro_base_url + pfam_domain_id
    headers = {'Accept': 'application/json'}

    # Log InterPro request
    logging.info(f'InterPro Request: {url}')

    try:
        # Make request to InterPro REST API
        r = requests.get(url, headers=headers, timeout=10)
        r.raise_for_status()
        json = r.json()

        # Extract domain long name
        metadata = json['metadata']
        long_name = metadata['name'].get('name') or metadata['name'].get('value')

        # Return domain long name
        return long_name
    
    # Return None if the request fails
    except Exception:
        logging.error(f'InterPro Request Failed: {url}')
        return None


def extract_pfam_name(domain_list):
    '''
    Extract first Pfam domain name from a VEP domains list.
    
    Parameters:
        domain_list: list of dicts. Output from VEP.
    
    Returns:
        str. Pfam domain ID or None.
    '''
    if not isinstance(domain_list, list):
        return None
    for d in domain_list:
        if isinstance(d, dict) and d.get('db') == 'Pfam':
            return d.get('name')
    return None


def vep_input(variants):
    '''
    Creates a VCFv4.0 of unique variants in given DataFrame.

    Parameters:
        variants: DataFrame. Contains Chromosome, Start, Ref, Alt
    
    Creates file processing/vep/input.vcf
    '''
    unique_variants = variants.drop_duplicates(subset='Submission')

    print(unique_variants)

    Path(vep_input_file).parent.mkdir(parents=True, exist_ok=True)

    with open(vep_input_file, 'w') as f:
        # VCF meta-information
        f.write('##fileformat=VCFv4.0\n')

        # Header line
        f.write('#CHROM\tPOS\tID\tREF\tALT\n')

        # Write variants
        for _, row in unique_variants.iterrows():
            f.write(f'{row['Submission']}\n')
            

def apptainer_vep(species, query=True, input_file='input.vcf', output_file='output.json'):
    '''
    Runs local VEP on input file for species.

    Parameters:
        species: str. Either 'homo_sapiens' or 'mus_musculus'.
        input_file: str. Input file in processing dictionary. Defaults to 'input.vcf'
        output_file: str. Output file in processing dictionary. Defaults to 'output.json'
    '''

    host_main = os.path.expanduser('~/VarLift/app/vep')

    cmd = ['apptainer', 'exec',
           '--bind', f'{host_main}:/main',
           f'{host_main}/ensembl-vep.sif', 'vep',
           '--dir_cache', '/main/vep_cache',
           '--input_file', f'/main/processing/{input_file}',
           '--output_file', f'/main/processing/{output_file}',
           '--json',
           '--cache',
           '--offline',
           '--force_overwrite',
           '--species', species,
           '--domains',
           '--biotype',
           '--symbol',
           '--numbers',
           '--coding_only']
    
    if species == 'homo_sapiens':
        cmd += ['--polyphen', 'b',
                '--pick', '--pick_order', 
                'mane_select,length']
    
    if species == 'mus_musculus' and query:
        cmd += ['--pick',
                '--pick_order', 'canonical,length']

    subprocess.run(cmd, check=True)


def parse_vep_json(input_file=vep_output_file):
    '''
    Parse VEP json output (Docker output) into a DataFrame.

    Parameters:
        inpute_file. str. path to json for parsing.
    
    Returns:
        DataFrame. VEP Output.
    '''

    all_dfs = []

    with open(input_file) as f:
        for line in f:
            if not line.strip():
                continue

            variant = json.loads(line)

            consequences = variant.get('transcript_consequences', [])
            if not consequences:
                logging.warning(f"No transcript consequences for {variant.get('input')}")
                continue

            # Flatten transcript consequences
            df = pd.json_normalize(consequences)

            # Keep only columns you care about
            keep = ['gene_symbol', 'transcript_id',
                    'polyphen_prediction', 'polyphen_score',
                    'amino_acids', 'protein_start',
                    'consequence_terms', 'exon',
                    'domains', 'codons',
                    'impact', 'biotype',
                    'hgvsp', 'hgvsc']
            df = df[[c for c in keep if c in df.columns]]

            # Extract Pfam domain
            if 'domains' in df.columns:
                df['Domain'] = df['domains'].apply(extract_pfam_name)
                df.drop(columns=['domains'], inplace=True)
            else:
                df['Domain'] = None

            # Attach submitted input (optional)
            df['Submission'] = variant.get('input')

            all_dfs.append(df)

    vep_df = pd.concat(all_dfs, ignore_index=True) if all_dfs else pd.DataFrame()
    vep_df = vep_df.replace({np.nan: None})

    return vep_df


def prepare_vep_output(vep_df):
    '''
    Cleans VEP output for desired protein information and format.

    Parameters:
        vep_df: DataFrame. Vep output from fetch_vep_data function.
    
    Returns:
        DataFrame. Containing Transcript ID, Gene Symbol, HGVS, Biotype, Exon Rank, Pfam Domain ID, Pfam Domain Name, Polyphen Prediction
        Polyphen Score, Molecular Consequence, Codon Switch, Amino Acids, refAA, varAA.
    '''
    # Split amino_acids to REFAA and VARAA if in X/Y format
    vep_df['protein_start'] = vep_df['protein_start'].astype('Int64')


    # Extract ref/varAA
    aa = vep_df['amino_acids'].astype('string').str.strip().str.extract(r"^(?P<refAA>[A-Z\*])(?:/(?P<varAA>[A-Z\*]))?$")

    vep_df['refAA'] = aa['refAA']
    vep_df['varAA'] = aa['varAA']

    # Handle synoymous variants
    vep_df["varAA"] = vep_df["varAA"].fillna("=")


    # Build the change string
    ref_str = vep_df['refAA'].astype('string')
    pos_str = vep_df['protein_start'].astype('string')
    var_str = vep_df['varAA'].astype('string')

    vep_df['amino_acids'] = ref_str.str.cat(pos_str).str.cat(var_str)


    # Merge consequence terms into a single string
    vep_df['consequence_terms'] = vep_df['consequence_terms'].apply(lambda x: ','.join(x) if isinstance(x, list) else x)


    # # Fetch domain names
    # Fetch Pfam domain names from InterPro (one API call per unique domain)
    if 'Domain' in vep_df.columns:
        unique_domains = [d for d in vep_df['Domain'].unique() if pd.notna(d)]
        name_map = {d: get_domain_name(d) for d in unique_domains}
        vep_df['domain_name'] = vep_df['Domain'].map(name_map)
    # for domain in vep_df['Domain'].unique():
    #     if domain == None:
    #         dom_name = None
    #     else:
    #         dom_name = get_domain_name(domain)
    #     vep_df.loc[vep_df['Domain'] == domain, 'domain_name'] = dom_name

    
    # Create protein DataFrame
#    keep = ['Name', 'Submission', 'gene_symbol', 'transcript_id', 'HGVS', 'biotype', 'exon', 'Domain', 'polyphen_prediction', 
#           'polyphen_score', 'consequence_terms', 'codons', 'amino_acids', 'refAA', 'protein_start', 'varAA']
    keep = ['Name', 'Submission', 'gene_symbol', 'transcript_id', 'HGVS', 'biotype', 'exon', 'Domain', 'domain_name', 'polyphen_prediction',
            'polyphen_score', 'consequence_terms', 'codons', 'amino_acids', 'refAA', 'protein_start', 'varAA']
    cols_present = [col for col in keep if col in vep_df.columns]
    protein_df = vep_df[cols_present].copy()

    # Update column names
    protein_df.rename(columns={'transcript_id': 'Transcript ID', 'gene_symbol': 'Gene Symbol', 'biotype': 'Biotype',
                       'exon': 'Exon Rank', 'Domain': 'Pfam Domain ID', 'domain_name': 'Pfam Domain Name',
                       'polyphen_prediction': 'Polyphen Prediction', 'polyphen_score': 'Polyphen Score',
                       'consequence_terms': 'Molecular Consequence', 'codons': 'Codon Switch', 'amino_acids': 'Amino Acids'},
                       inplace=True,  errors='ignore')

#    protein_df.rename(columns={'transcript_id': 'Transcript ID', 'gene_symbol': 'Gene Symbol', 'biotype': 'Biotype', 
#                       'exon': 'Exon Rank', 'Domain': 'Pfam Domain ID', 
#                       'polyphen_prediction': 'Polyphen Prediction', 'polyphen_score': 'Polyphen Score', 
#                       'consequence_terms': 'Molecular Consequence', 'codons': 'Codon Switch', 'amino_acids': 'Amino Acids'}, 
#                       inplace=True,  errors='ignore')
    
    return protein_df


def run_vep(variants, species, query=True):
    '''
    Runs VEP using help functions.

    Parameters:
        variants: DataFrame. Contains Chromosome, Start, Ref, Alt
        species: str. Either 'homo_sapiens' or 'mus_musculus'.
    
    Returns:
        DataFrame. VEP Output.
    '''

    vep_input(variants)

    try:
        apptainer_vep(species, query=query)
        
        vep_df = parse_vep_json()

    except:
        logging.error('VEP Request Failed')
        vep_df = pd.DataFrame(columns=['gene_symbol','transcript_id','amino_acids','protein_start',
                                       'consequence_terms','exon','codons','impact','biotype',
                                       'Domain','Submission','Name'])
    
    return vep_df

