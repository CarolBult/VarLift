import re
import pandas as pd

from app.clinvar_summary import variants_for_genes


def process_batch_query(input):
    '''
    Processes a string resulting from the batch input of variants with formatting:
    chromosome:start-end:ref/allele

    Parameters:
        input: string. Multiline string containg all variants to process.

    Returns: DataFrame. Containing Chromosome, Start, Ref, Alt, for each variant.
    '''
    # Regex pattern to match each line: gene:chromosome:start-end:ref/allele
    pattern = r'^(\w+):(\d+)-(\d+):([ACGTacgt]+)/([ACGTacgt]+)$'

    variants = []
    # Go through each line
    for line in input.splitlines():
        if line.strip() == '':
            continue
        
        # Extract inputted variant fields
        match = re.match(pattern, line)
        if match:
            variants.append({'Name': match.group(0),
                             'Chromosome': match.group(1),
                             'Start': int(match.group(2)),
                             'Stop': int(match.group(3)),
                             'Ref': match.group(4),
                             'Alt': match.group(5)})
        
    # Create pandas dataframe
    variants = pd.DataFrame(variants)

    return variants


def prepare_submission(variants):
    '''
    Converts a DataFrame of variants into a list of HGHVS notations for VEP querying
    
    Parameters:
        variants: DataFrame. containing Chromosome, Start, Ref, Alt.

    Returns:
        tuple. list of strings and dictionary. HGVS notations (chrom:g.startRef>Alt) and mapping of HGVS to names.
    '''

    variants['Submission'] = (variants['Chromosome'].astype(str) + '\t' + variants['Start'].astype(str) + '\t.\t' + variants['Ref'] + '\t' + variants['Alt'])

    try: submission_HGVS_map = dict(zip(variants['Submission'], variants['Name']))

    except: submission_HGVS_map = None

    return variants, submission_HGVS_map


def process_human_fetch(hum_gene_df, mouse_protein_df, gene_input_df):

    merged_df = mouse_protein_df.merge(gene_input_df, left_on='Gene Symbol', right_on='Mus Gene')

    # Only the homologous genes' rows can match the merge below
    variant_summary_grch38 = variants_for_genes(merged_df['Hum Gene'].unique().tolist())

    homologous_variants = variant_summary_grch38.merge(
        merged_df[['Hum Gene', 'refAA', 'varAA']],
        left_on=['GeneSymbol', 'refAA', 'varAA'],
        right_on=['Hum Gene', 'refAA', 'varAA'],
        how='inner'
    )

    homologous_variants = homologous_variants.merge(
        hum_gene_df[['Gene Symbol', 'Strand']],
        left_on=['GeneSymbol'],
        right_on=['Gene Symbol'],
        how='inner'
    )

    genomic_pattern = re.compile(r'(?i)\b[crgm]\.\d+(?:[_+-]\d+)?(?P<ref>[ACGT])>(?P<alt>[ACGT])')

    # Extract REF/ALT as new columns
    homologous_variants[['Ref', 'Alt']] = homologous_variants['Name'].str.extract(genomic_pattern)
    homologous_variants.drop(columns=['ReferenceAllele', 'AlternateAllele'], inplace=True)


    # Map alleles to correct strand
    compliment_map = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    mask = homologous_variants['Strand'] == '-'

    homologous_variants.loc[mask, 'Ref'] = (homologous_variants.loc[mask, 'Ref'].map(compliment_map))
    homologous_variants.loc[mask, 'Alt'] = (homologous_variants.loc[mask, 'Alt'].map(compliment_map))


    # Extract phenotypes
    phenotype_df = homologous_variants[['Name', 'PhenotypeIDS', 'PhenotypeList']].copy()
    phenotype_df['PhenotypeIDS'] = phenotype_df['PhenotypeIDS'].str.split('|')
    phenotype_df['PhenotypeList'] = phenotype_df['PhenotypeList'].str.split('|')

    phenotype_df = phenotype_df.explode(['PhenotypeIDS', 'PhenotypeList'], ignore_index=True)
    phenotype_df['HP_ID'] = phenotype_df['PhenotypeIDS'].str.extract(r'(HP:\d+)', expand=False)
    phenotype_df = phenotype_df.dropna(subset=['HP_ID'])
    phenotype_df['Phenotypes'] = phenotype_df['HP_ID'] + ',' + phenotype_df['PhenotypeList']


    phenotype_df = (
        phenotype_df
        .groupby('Name', sort=False)['Phenotypes']
        .apply(list)
        .reset_index()
    )

    return homologous_variants, phenotype_df