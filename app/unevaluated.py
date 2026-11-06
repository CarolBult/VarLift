'''
Report submitted variants that VarLift cannot evaluate.

VarLift matches variants by amino-acid position, so variants with no coding
consequence (intronic, splice-site, UTR, intergenic) never enter the scoring
pipeline. They can still have disease associations and mouse alleles, so this
module annotates them with their non-coding consequence and reports what is
known about them.
'''
import pandas as pd

from app.vep_util import run_vep
from app.gene_util import homology_df
from app.mouse_coverage import mgi_alleles_at


def _reason_for(consequence):
    if not isinstance(consequence, str) or consequence == '':
        return 'VEP could not annotate this variant'
    if 'intergenic' in consequence:
        return 'Not within a gene, so there is no amino-acid change to match'
    return f'No amino-acid change ({consequence}), so the variant cannot be matched by protein position'


def _annotate(variants, species):
    '''
    One row per variant: the gene, transcript and consequence VEP picks without
    the coding-only filter. Variants VEP could not annotate keep empty fields.
    '''
    vep_df = run_vep(variants, species, query=True, coding_only=False)

    cols = ['Submission', 'gene_symbol', 'transcript_id', 'biotype', 'consequence_terms', 'exon', 'intron']
    vep_df = vep_df.reindex(columns=cols).drop_duplicates(subset=['Submission'])
    vep_df['consequence_terms'] = vep_df['consequence_terms'].apply(
        lambda x: ','.join(x) if isinstance(x, list) else x)

    annotated = variants[['Name', 'Submission']].merge(vep_df, on='Submission', how='left')
    annotated = annotated.rename(columns={'Name': 'Input Name', 'gene_symbol': 'Gene Symbol',
                                          'transcript_id': 'Transcript ID', 'biotype': 'Biotype',
                                          'consequence_terms': 'Molecular Consequence',
                                          'exon': 'Exon', 'intron': 'Intron'})
    annotated['Reason'] = annotated['Molecular Consequence'].apply(_reason_for)
    return annotated


def _homologs(genes, species):
    '''Map gene symbols to their homologs in the other species, joined by ", ".'''
    src, dst = ('HumGeneSymbol', 'MusGeneSymbol') if species == 'homo_sapiens' else ('MusGeneSymbol', 'HumGeneSymbol')
    hom = homology_df[homology_df[src].isin(genes)]
    return hom.groupby(src)[dst].apply(lambda s: ', '.join(sorted(set(s.dropna())))).to_dict()


def unevaluated_variants(variants, coding_submissions, species):
    '''
    Annotate the submitted variants that have no coding consequence.

    Parameters:
        variants: DataFrame. Submitted variants with Name, Submission, Chromosome, Start, Stop, Ref, Alt.
        coding_submissions: iterable of str. Submissions that VEP found coding consequences for.
        species: str. Either 'homo_sapiens' or 'mus_musculus'.

    Returns:
        DataFrame. One row per variant (per matching MGI allele for mouse) with its gene,
        consequence, the reason it was not evaluated, and known disease associations.
    '''
    # assign_clinvar pairs rows with names by position, so the index must start at 0
    rest = variants[~variants['Submission'].isin(set(coding_submissions))].reset_index(drop=True)
    if rest.empty:
        return pd.DataFrame()

    report = _annotate(rest, species)

    homolog_col = 'Mouse Homolog' if species == 'homo_sapiens' else 'Human Homolog'
    homologs = _homologs(report['Gene Symbol'].dropna().unique(), species)
    report.insert(report.columns.get_loc('Gene Symbol') + 1, homolog_col, report['Gene Symbol'].map(homologs))

    if species == 'homo_sapiens':
        from app.main_util import human_variant_diseases
        diseases = human_variant_diseases(rest)[['Name', 'Associated Diseases', 'MONDO']]
        report = report.merge(diseases.rename(columns={'Name': 'Input Name'}), on='Input Name', how='left')
    else:
        alleles = mgi_alleles_at(rest)
        report = report.merge(alleles, on='Input Name', how='left')

    report = report.drop(columns=['Submission'])
    return report.astype(object).where(pd.notnull(report), None)
