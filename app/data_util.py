import json
import re
import yaml
import obonet
import requests
import pandas as pd
from urllib.parse import unquote
from pathlib import Path


# Load Config
with open('app/config.yaml') as f:
    config = yaml.safe_load(f)

# Create AA map
aa_multi_to_single = {'Ala':'A', 'Arg':'R', 'Asn':'N', 'Asp':'D',
                      'Cys':'C', 'Glu':'E', 'Gln':'Q', 'Gly':'G',
                      'His':'H', 'Ile':'I', 'Leu':'L', 'Lys':'K',
                      'Met':'M', 'Phe':'F', 'Pro':'P', 'Ser':'S',
                      'Thr':'T', 'Trp':'W', 'Tyr':'Y', 'Val':'V',
                      'Ter':'*', '=':'='}


def download_file(input_url, output_path):
    '''
    Downloads a file to a given location from a url.

    Paramters:
        input_url: str. URL to request from.
        output_path: str. File to save URL to.
    '''
    try:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        
        with requests.get(input_url, stream=True) as r:
            r.raise_for_status()
            with open(output_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
        print(f"Downloaded '{output_path}' successfully.")

    except requests.exceptions.RequestException as e:
        print(f'Download failed: {e}')


def parse_mondo_obo(mondo_path):
    '''
    Parses the MONDO OBO file to create a mapping of disease synonyms to their primary names.
    Creates a JSON file 'mondo_term_map.json' mapping MONDO ID to term name and mondo_xref_map
    mapping ref ids to mondo ids in the 'data/MONDO/' directory.

    Parameters:
        mondo_path. str. Path to mondo file parse.
    '''
    graph = obonet.read_obo(mondo_path)

    mondo_term_map = {}
    mondo_xref_map = {}
    for term, data in graph.nodes(data=True):
        name = data.get("name")

        id = term

        mondo_term_map[id] = name

        for x_id in data.get('xref', []):
            x_ref = x_id.split(' ')[0]
            mondo_xref_map[x_ref] = id

    json.dump(mondo_term_map, open(config['paths']['mondo_term_map'], 'w'), indent=4)
    json.dump(mondo_xref_map, open(config['paths']['mondo_xref_map'], 'w'), indent=4)


def parse_ensembl_gff3(gff3_path, species):
    '''
    Parses Ensembl GFF3 file to extract gene information and saves it as a parquet file.

    Parameters:
        gff3_path. str. Path to gff3 file.
        species. str. 'human' or 'mouse'.
    '''
    cols = ['seqid','source','type','start','end','score','strand','phase','attributes'] # Current gff3 formatting
    gene_df = pd.read_csv(gff3_path, sep='\t', comment='#', compression='gzip', header=None, names=cols, low_memory=False)
    gene_df = gene_df[(gene_df['type'] == 'gene')]

    # Extract fields using regex
    gene_df['Name'] = gene_df['attributes'].str.extract(r'Name=([^;]+)')
    gene_df['description'] = gene_df['attributes'].str.extract(r'description=([^;]+)')[0].apply(
    lambda x: unquote(x) if isinstance(x, str) else x)
    #gene_df['description'] = gene_df['attributes'].str.extract(r'description=([^;]+)')
    gene_df['biotype'] = gene_df['attributes'].str.extract(r'biotype=([^;]+)')
    gene_df['gene_id'] = gene_df['attributes'].str.extract(r'ID=gene:([^;]+)')

    # Keep only needed columns
    gene_df = gene_df[['Name', 'description', 'biotype', 'seqid', 'start', 'end', 'strand', 'gene_id']]

    gene_df.to_parquet(config['paths'][f'{species}_gene_pqt'])


def parse_ncbi_variant_summary(summary_path, output_path):
    '''
    Parses NCBI TSV file as a source of all human variants and 
    extracts relevation information and saves it as a csv.

    Parameters:
        summary_path. str. Path to TSV file.
        output_path. str. Path to final, cleaned TSV file.
    '''
    cols = ['#AlleleID', 'Type', 'Name', 'GeneSymbol', 'HGNC_ID', 
        'ClinSigSimple', 'PhenotypeIDS', 'PhenotypeList', 
        'Assembly', 'Chromosome', 'Start', 'Stop', 
        'ReferenceAllele', 'AlternateAllele']

    var = pd.read_csv(summary_path, sep='\t', usecols=cols)

    var_grch38 = var[(var['Assembly'] == 'GRCh38')] #  & (var['Type'] == 'single nucleotide variant')

    var_grch38.rename(columns={'#AlleleID': 'AlleleID'}, inplace=True)


    prot_pattern = re.compile(r'p\.([A-Za-z]{3})\d+([A-Za-z]{3}|=)')

    # Extract REF/ALT proteins as new columns
    var_grch38[['refAA', 'varAA']] = var_grch38['Name'].str.extract(prot_pattern)


    var_grch38['refAA'] = var_grch38['refAA'].map(aa_multi_to_single)
    var_grch38['varAA'] = var_grch38['varAA'].map(aa_multi_to_single)

    var_grch38.to_csv(output_path, sep='\t', index=False)
