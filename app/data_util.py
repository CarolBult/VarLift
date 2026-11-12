import gzip
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


def count_header_lines(path):
    '''
    Counts the '#' comment lines at the top of a file, for use as pandas skiprows.
    comment='#' would also cut data lines at any '#', such as transgene symbols
    like Tg(Lck-Notch3)#Issc, losing every column after it.

    Parameters:
        path: str. Path to file, optionally gzipped.

    Returns:
        int. Number of leading lines that start with '#'.
    '''
    opener = gzip.open if str(path).endswith('.gz') else open

    count = 0
    with opener(path, 'rt') as f:
        for line in f:
            if not line.startswith('#'):
                break
            count += 1

    return count


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


def parse_ensembl_exons(gff3_path, species):
    '''
    Parses Ensembl GFF3 file to extract exon coordinates per transcript and saves them
    as a parquet file. Used to place intronic variants relative to the nearest exon.

    Parameters:
        gff3_path. str. Path to gff3 file.
        species. str. 'human' or 'mouse'.
    '''
    cols = ['seqid','source','type','start','end','score','strand','phase','attributes']
    gff_df = pd.read_csv(gff3_path, sep='\t', comment='#', compression='gzip', header=None, names=cols,
                         usecols=['seqid', 'type', 'start', 'end', 'strand', 'attributes'], low_memory=False,
                         dtype={'seqid': str})

    # Transcript -> gene, so a gene's exons can be found (e.g. whether a variant lies in an exon of an overlapping gene)
    transcript_rows = gff_df.loc[gff_df['attributes'].str.startswith('ID=transcript:', na=False), 'attributes']
    transcripts = transcript_rows.str.extract(r'^ID=transcript:([^;]+);Parent=gene:([^;]+)').dropna()
    transcript_gene = dict(zip(transcripts[0], transcripts[1]))

    exon_df = gff_df[gff_df['type'] == 'exon'].copy()

    exon_df['transcript_id'] = exon_df['attributes'].str.extract(r'Parent=transcript:([^;]+)', expand=False)
    exon_df['rank'] = exon_df['attributes'].str.extract(r'rank=(\d+)', expand=False).astype(int)
    exon_df['gene_id'] = exon_df['transcript_id'].map(transcript_gene)

    exon_df = exon_df[['transcript_id', 'rank', 'start', 'end', 'strand', 'seqid', 'gene_id']].dropna(subset=['transcript_id'])

    exon_df.to_parquet(config['paths'][f'{species}_exon_pqt'], index=False)


def _read_ensembl_fasta(fasta_path, id_pattern):
    '''(transcript ID, sequence) for each record of a gzipped Ensembl FASTA, the ID taken from the header by id_pattern.'''
    rows = []
    transcript, sequence = None, []
    with gzip.open(fasta_path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if transcript:
                    rows.append((transcript, ''.join(sequence)))
                match = re.search(id_pattern, line)
                transcript, sequence = (match.group(1) if match else None), []
            else:
                sequence.append(line.strip())
        if transcript:
            rows.append((transcript, ''.join(sequence)))
    return pd.DataFrame(rows, columns=['transcript_id', 'sequence']).drop_duplicates('transcript_id')


def parse_ensembl_peptides(pep_path, species):
    '''
    Parses an Ensembl peptide FASTA to extract each transcript's protein sequence and
    saves it as a parquet file. Used to align homologous proteins so residue positions
    can be compared across species.

    Parameters:
        pep_path. str. Path to gzipped pep.all FASTA file.
        species. str. 'human' or 'mouse'.
    '''
    protein_df = _read_ensembl_fasta(pep_path, r'transcript:(\S+?)(?:\.\d+)?(?:\s|$)')
    protein_df.to_parquet(config['paths'][f'{species}_protein_pqt'], index=False)


def parse_ensembl_cds(cds_path, species):
    '''
    Parses an Ensembl CDS FASTA (headers start with the versioned transcript ID) to
    extract each transcript's coding sequence and saves it as a parquet file. Used to
    show the codons beneath a protein alignment.

    Parameters:
        cds_path. str. Path to gzipped cds.all FASTA file.
        species. str. 'human' or 'mouse'.
    '''
    cds_df = _read_ensembl_fasta(cds_path, r'^>(\S+?)(?:\.\d+)?(?:\s|$)')
    cds_df.to_parquet(config['paths'][f'{species}_cds_pqt'], index=False)


def parse_ensembl_domains(mysql_dir, species):
    '''
    Parses the Pfam protein features of an Ensembl core database dump (the
    analysis, translation, transcript and protein_feature tables) into each
    transcript's Pfam domains with their residue ranges, saved as a parquet file.
    Used to draw domains on the transcript alignment views.

    Parameters:
        mysql_dir. str. Directory holding the four gzipped table dumps.
        species. str. 'human' or 'mouse'.
    '''
    def table(name, columns, usecols):
        return pd.read_csv(Path(mysql_dir) / f'{name}.txt.gz', sep='\t', header=None, names=columns,
                           usecols=usecols, quoting=3, na_values='\\N', keep_default_na=False, low_memory=False)

    analysis = table('analysis', ['analysis_id', 'created', 'logic_name', 'db', 'db_version', 'db_file', 'program',
                                  'program_version', 'program_file', 'parameters', 'module', 'module_version',
                                  'gff_source', 'gff_feature'], ['analysis_id', 'logic_name'])
    pfam_ids = set(analysis.loc[analysis['logic_name'].str.lower() == 'pfam', 'analysis_id'])

    features = table('protein_feature', ['protein_feature_id', 'translation_id', 'seq_start', 'seq_end', 'hit_start',
                                         'hit_end', 'hit_name', 'analysis_id', 'score', 'evalue', 'perc_ident',
                                         'external_data', 'hit_description', 'cigar_line', 'align_type'],
                     ['translation_id', 'seq_start', 'seq_end', 'hit_name', 'analysis_id', 'hit_description'])
    features = features[features['analysis_id'].isin(pfam_ids)]

    translations = table('translation', ['translation_id', 'transcript_id', 'seq_start', 'start_exon_id', 'seq_end',
                                         'end_exon_id', 'stable_id', 'version', 'created_date', 'modified_date'],
                         ['translation_id', 'transcript_id'])
    transcripts = table('transcript', ['transcript_id', 'gene_id', 'analysis_id', 'seq_region_id', 'seq_region_start',
                                       'seq_region_end', 'seq_region_strand', 'display_xref_id', 'source', 'biotype',
                                       'description', 'is_current', 'canonical_translation_id', 'stable_id', 'version',
                                       'created_date', 'modified_date'], ['transcript_id', 'stable_id'])

    domains = (features.merge(translations, on='translation_id')
                       .merge(transcripts.rename(columns={'transcript_id': 'tid', 'stable_id': 'transcript_id'}),
                              left_on='transcript_id', right_on='tid', suffixes=('_internal', '')))
    domains = domains.rename(columns={'hit_name': 'domain_id', 'seq_start': 'start', 'seq_end': 'end'})

    # Pfam's descriptive family names (Ensembl stores the short IDs, e.g. PK_Tyr_Ser-Thr) and clans,
    # which group related families such as the two protein kinase models PF00069 and PF07714
    clans = pd.read_csv(config['paths']['pfam_clans'], sep='\t', header=None, compression='gzip', quoting=3,
                        names=['domain_id', 'clan', 'clan_name', 'short_name', 'description'])
    domains = domains.merge(clans[['domain_id', 'clan', 'description']], on='domain_id', how='left')
    domains['name'] = domains['description'].fillna(domains['hit_description'])
    domains = (domains[['transcript_id', 'domain_id', 'name', 'clan', 'start', 'end']]
               .drop_duplicates().sort_values(['transcript_id', 'start']))
    domains.to_parquet(config['paths'][f'{species}_domain_pqt'], index=False)


def parse_ncbi_variant_summary(summary_path, output_path):
    '''
    Parses NCBI TSV file as a source of all human variants and
    extracts relevation information and saves it as a parquet.

    Parameters:
        summary_path. str. Path to TSV file.
        output_path. str. Path to final, cleaned parquet file.
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

    write_variant_summary_parquet(var_grch38, output_path)


def write_variant_summary_parquet(var_df, output_path):
    '''
    Saves the parsed variant summary as a parquet sorted by position, so lookups
    by position or gene only read the few row groups that hold matching rows.

    Parameters:
        var_df. DataFrame. Parsed GRCh38 variant summary.
        output_path. str. Path to parquet file.
    '''
    var_df = var_df.copy()

    # Chromosome mixes numbers and letters; store as text so parquet accepts it
    var_df['Chromosome'] = var_df['Chromosome'].map(lambda c: c if pd.isna(c) else str(c))

    # Remember the original file order so results can be returned in it
    var_df['_order'] = range(len(var_df))

    var_df = var_df.sort_values(['Chromosome', 'Start'], kind='stable')
    var_df.to_parquet(output_path, index=False, row_group_size=20000, compression='zstd')
