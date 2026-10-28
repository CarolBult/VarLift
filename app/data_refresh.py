import yaml
from data_util import download_file, parse_mondo_obo, parse_ensembl_gff3, parse_ncbi_variant_summary


# Load Config
with open('app/config.yaml') as f:
    config = yaml.safe_load(f)


def refresh_mondo():
    '''
    Downloads the mondo OBO file and parse it to create a mapping of Mondo IDs to names.
    '''
    url = config['data']['mondo_obo']
    output_path = config['paths']['mondo_term_map']

    download_file(url, output_path)
    parse_mondo_obo(output_path)


def refresh_mgi():
    '''
    Downloads MGI data files.
    '''
    download_file(config['data']['mgi_alleles'], config['paths']['mgi_alleles'])
    download_file(config['data']['mgi_homology'], config['paths']['mgi_homology'])
    download_file(config['data']['mgi_disease'], config['paths']['mgi_disease'])
    download_file(config['data']['mgi_phenotype'], config['paths']['mgi_phenotype'])

def refresh_mp_hp_mapping():
    '''
    Downloads the Monarch mapping-commons MP↔HP SSSOM mapping file.
    Used for mouse-to-human phenotype equivalency scoring.
    '''
    download_file(config['data']['mp_hp_mapping'], config['paths']['mp_hp_mapping'])


def refresh_mp_hp_mapping():
    '''
    Downloads the Monarch mapping-commons MP↔HP SSSOM mapping file.
    Used for mouse-to-human phenotype equivalency scoring.
    '''
    download_file(config['data']['mp_hp_mapping'], config['paths']['mp_hp_mapping'])


def refresh_clinvar():
    '''
    Downloads ClinVar files.
    '''
    url_vcf = config['data']['clinvar_vcf']
    url_tbi = config['data']['clinvar_tbi']
    out_vcf = config['paths']['clinvar_vcf']
    out_tbi = config['paths']['clinvar_tbi']

    download_file(url_vcf, out_vcf)
    download_file(url_tbi, out_tbi)


def refresh_ensembl():
    '''
    Downloads Ensembl GFF3 files.
    * Note that VEP Caches are not included here.
    '''
    version = config['data']['ensemble_version']

    url_hum = config['data']['ensembl_hum_gff3'].format(ensemble_version=version)
    url_mus = config['data']['ensembl_mus_gff3'].format(ensemble_version=version)
    out_hum = config['paths']['ensembl_hum_gff3'].format(ensemble_version=version)
    out_mus = config['paths']['ensembl_mus_gff3'].format(ensemble_version=version)

    download_file(url_hum, out_hum)
    download_file(url_mus, out_mus)

    parse_ensembl_gff3(out_hum, 'human')
    parse_ensembl_gff3(out_mus, 'mouse')


def refresh_ncbi():
    '''
    Download and parse NCBI variant summary.
    '''
    url = config['data']['ncbi_var_summary']
    output_path = config['paths']['ncbi_var_summary']

    download_file(url, output_path)

    input_path = config['paths']['ncbi_var_summary']
    output_path = config['paths']['ncbi_var_summary_grch38']

    parse_ncbi_variant_summary(input_path, output_path)
    

def refresh_data():
    '''
    Refreshes all data with updated downloads. Both downloads and parses all data sources.
    '''
    refresh_mondo()
    refresh_mgi()
    refresh_mp_hp_mapping()
    refresh_clinvar()
    refresh_ensembl()
    refresh_ncbi()



if __name__ == '__main__':
    refresh_data()
