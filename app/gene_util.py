import re
import yaml
import logging
import requests
import pandas as pd

from app.data_util import count_header_lines

# Load Config
with open('app/config.yaml') as f:
    config = yaml.safe_load(f)

# Load Ensembl URL
ensembl_base_url = config['api']['ensembl_base_url']

# Load Gene Parquets
mouse_gene_pqt = pd.read_parquet(config['paths']['mouse_gene_pqt'])
human_gene_pqt = pd.read_parquet(config['paths']['human_gene_pqt'])

# Load homology
homology_df = pd.read_csv(config['paths']['mgi_homology'],
                          sep='\t', names=['MGI_ID', 'MusGeneSymbol', 'MusEntrezGeneID', 
                                           'Mus_HGNC_ID', 'HumGeneSymbol', 'HumEntrezGeneID'])

# Load Mus alleles
mus_alleles_df = pd.read_csv(config['paths']['mgi_alleles'], sep='\t', compression='gzip', low_memory=False,
                             skiprows=count_header_lines(config['paths']['mgi_alleles']))


def get_gene_info(gene, species='human'):
    '''
    Query Ensembl REST API and extract info for a given gene.

    Parameters:
        gene: string. Gene symbol to extract info for.
        species: string. Species to extract gene from. Defaults to 'human'.

    Returns:
        dictionary. Containins HGNC ID, Esnembl ID, name, strand, biotype, and description.
    '''
    # Construct URL
    url = f'{ensembl_base_url}/lookup/symbol/{species}/{gene}?content-type=application/json'

    # Log Ensembl request
    logging.info(f'Ensemble Gene Info Request: {url}')

    try:
        # Make request to Ensemble REST API
        request = requests.get(url, headers={'Content-Type': 'application/json'})
        request.raise_for_status()
        data = request.json()

        # Extract relevant data
        ens_id = data.get('id')
        name = data.get('display_name')
        biotype = data.get('biotype')
        strand = data.get('strand')
        match = re.match(r'^(.*?) \[Source:(.+) Symbol;Acc:(.+)\]$', data.get('description')) # Split string to access description and HGNC ID
        description = match.group(1).strip()
        
        # Different MGI IDs for mouse and HGNC IDs for human
        if species == 'human':
            hgnc_id = match.group(3).strip()

            # Create gene dictionary
            gene_info = {'HGNC_id': hgnc_id, 'ENSMBL_id': ens_id, 'name': name, 'strand': strand, 'biotype': biotype, 'description': description}
        
        if species == 'mouse':
            hgnc_id = match.group(3).strip()

            # Create gene dictionary
            gene_info = {'MGI_id': hgnc_id, 'ENSMBL_id': ens_id, 'name': name, 'strand': strand, 'biotype': biotype, 'description': description}
        
    except Exception as e:
        logging.error(f'Ensembl Gene Info Request Failed: {e}')
        ValueError(f'Ensembl Gene Info Request Failed: {e}')

    return gene_info


def fetch_gene_info(genes, species='human'):
    '''
    Fetch gene info from saved parquets in /data/Ensembl/ for a lsit of gene symbols.

    Parameters:
        genes: list of strings. Gene symbols to extract info for.
        species: string. 'human' or 'mouse'. Defaults to 'human'.

    Returns:
        DataFrame. Containins HGNC ID, Ensembl ID, name, strand, biotype, and description.
    '''
    # Assign correct parquet based on species
    if species == 'human': pqt = human_gene_pqt
    elif species == 'mouse': pqt = mouse_gene_pqt

    # Subset parquet to requested genes
    subset = pqt[pqt['Name'].isin(genes)]

    # Preserve order of input genes
    subset = subset.set_index('Name').reindex(list(genes)).reset_index()

    # Extract description and accession
    subset[['Description', 'Accession']] = subset['description'].str.extract(r'^(.*?) \[Source:.+ Symbol.+Acc:(.+)\]$')

    # Clean DataFrame
    subset.drop('description', axis=1, inplace=True)
    
    subset.rename(columns={'Name': 'Gene Symbol', 'biotype': 'Biotype', 'seqid': 'Chromosome', 'start': 'Start', 
                           'end': 'End', 'strand': 'Strand', 'gene_id': 'Ensembl_ID'}, inplace=True)
    
    subset = subset[['Gene Symbol', 'Description', 'Biotype', 'Chromosome', 'Start', 'End', 'Strand', 'Ensembl_ID', 'Accession']]

    return subset


def fetch_homologous_gene(input_mapping_df, species):
    '''
    '''
    if species == 'mouse':
        input_mapping_df = homology_df.merge(input_mapping_df, left_on='HumGeneSymbol', right_on='Hum Gene', how='inner')\
                                    .rename(columns={'MusGeneSymbol': 'Mus Gene'}).copy()
        input_mapping_df.drop(['HumGeneSymbol', 'Mus_HGNC_ID', 'HumEntrezGeneID', 'MusEntrezGeneID'], axis=1, inplace=True)

    elif species == 'human':
        input_mapping_df = homology_df.merge(input_mapping_df, left_on='MusGeneSymbol', right_on='Mus Gene', how='inner')\
                                    .rename(columns={'HumGeneSymbol': 'Hum Gene'}).copy()
        input_mapping_df.drop(['MusGeneSymbol', 'Mus_HGNC_ID', 'HumEntrezGeneID', 'MusEntrezGeneID'], axis=1, inplace=True)

    return input_mapping_df



def fetch_mus_alleles(MGI_gene_ids):
    '''
    Fetches all Mouse alleles with human models for all genes from list of MGI_gene_ids.
    
    Parameters:
        MGI_gene_ids: list of str. A list of all MGI_gene_ids to retur strings for.
    
    Returns:
        DataFrame. Containing Gene, AlleleID, AlleleSymbol, Chromosome, Start, End, Ref, Alt, Molecular Consequence.
    '''
    # Only retain variants with logged genomic information
    mouse_select = mus_alleles_df[(mus_alleles_df['AlleleAssociatedGeneId'].isin(MGI_gene_ids))]

    # Create allele DataFrame
    mouse_allele_df = pd.DataFrame()
    old_cols = ['AlleleAssociatedGeneSymbol', 'AlleleId', 'AlleleSymbol', 'Chromosome', 'StartPosition', 
                'EndPosition', 'SequenceOfReference', 'SequenceOfVariant'] # Column names in mouse_select
    new_cols = ['Gene Symbol', 'AlleleID', 'AlleleSymbol', 'Chromosome', 'Start', 'End', 'Ref', 'Alt'] # Corresponding column names for mouse_allele_df
    mouse_allele_df[new_cols] = mouse_select[old_cols] # Update column names

    # Extract primary Molecular Consequence term
    mouse_allele_df['Molecular Consequence'] = mouse_select['MostSevereConsequenceName'].str.split(',').str[0]

    return mouse_allele_df

def process_mus_alleles(mouse_allele_df):
    '''
    Filters alleles to one with genomic coordinates listed and makes this integers.

    Parameters:
        mouse_allele_df: DataFrame. MGI allele information.
    
    Returns:
        DataFrame. Contains mouse allele information to continue through pipeline.
    '''
    
    mouse_prt_df = mouse_allele_df[mouse_allele_df['Start'].notnull()]

    # Convert coordinates to integers
    mouse_prt_df['Start'] = mouse_prt_df['Start'].astype(int)
    mouse_prt_df['End'] = mouse_prt_df['End'].astype(int)

    return mouse_prt_df