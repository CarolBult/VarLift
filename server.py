from fastapi import FastAPI
from pydantic import BaseModel
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.encoders import jsonable_encoder

import logging

from app.main_util import hvar_query, mvar_fetch, hvar_query_score, mvar_query, hvar_fetch, mvar_query_score, NoCodingVariantsError, phenotype_similarity
from app.processing_util import process_batch_query
from app.transcript_alignment import transcript_alignment, transcript_overview
from app.protein_alignment import SequenceUnavailable

ORGANISM_LABELS = {'human': 'Human (GRCh38)', 'mouse': 'Mouse (GRCm39)'}

def sort_phenotypes_by_score(phenotype_df, score_df):
    '''
    Reorder phenotype_df so the most equivalent alleles appear first: by each allele's
    Best Tier, then by Disease Match. Alleles with no score go to the bottom.
    '''
    if phenotype_df.empty or score_df.empty or 'AlleleID' not in phenotype_df.columns or 'AlleleID' not in score_df.columns:
        return phenotype_df

    per_allele = score_df.assign(_disease=score_df['Disease Match'].eq(True)).groupby('AlleleID').agg(
        tier=('Best Tier', 'min'), disease=('_disease', 'max'))
    allele_order = per_allele.sort_values(['tier', 'disease'], ascending=[True, False], kind='stable').index.tolist()
    rank = {aid: i for i, aid in enumerate(allele_order)}

    phenotype_df = phenotype_df.copy()
    phenotype_df['_sort'] = phenotype_df['AlleleID'].map(rank).fillna(len(allele_order))
    phenotype_df = (
        phenotype_df
        .sort_values('_sort', kind='stable')
        .drop(columns='_sort')
        .reset_index(drop=True)
    )
    return phenotype_df

app = FastAPI()

app.mount('/static', StaticFiles(directory='static'), name='static')
# User guide (docs/USER_GUIDE.md and its screenshots), shown at /help
app.mount('/guide', StaticFiles(directory='docs'), name='guide')

@app.get('/')
def serve_dashboard():
    return FileResponse('templates/dashboard.html')

@app.get('/help')
def serve_help():
    return FileResponse('templates/help.html')

@app.get('/alignment')
def alignment(query_tx: str, query_species: str, target_tx: str, target_species: str,
              query_residue: int | None = None, target_residue: int | None = None, flank: int | None = None):
    '''
    Codon-aware alignment of a query transcript and a homologous transcript, for the
    Align view of the Variants Evaluated table. flank limits the columns returned to
    that many residues either side of query_residue; omit it for the whole alignment.
    '''
    try:
        return transcript_alignment(query_tx, query_species, target_tx, target_species,
                                    query_residue, target_residue, flank)
    except SequenceUnavailable:
        return JSONResponse(status_code=404, content={
            'error': f'No protein sequence is available for {query_tx} or {target_tx}.'})
    except Exception:
        logging.exception('Alignment failed')
        return JSONResponse(status_code=500, content={'error': 'The alignment could not be made.'})

class OverviewInput(BaseModel):
    query_tx: str
    query_species: str
    query_residue: int | None = None
    target_species: str
    targets: list[tuple[str, int | None]]

@app.post('/transcript_overview')
def overview(data: OverviewInput):
    '''
    All of an allele's transcripts aligned to the query transcript, for the transcript
    overview opened from the Aligned Position column of the scoring table.
    '''
    try:
        return transcript_overview(data.query_tx, data.query_species, data.query_residue,
                                   data.targets, data.target_species)
    except SequenceUnavailable:
        return JSONResponse(status_code=404, content={'error': f'No protein sequence is available for {data.query_tx}.'})
    except Exception:
        logging.exception('Transcript overview failed')
        return JSONResponse(status_code=500, content={'error': 'The transcripts could not be aligned.'})

class VariantInput(BaseModel):
    variants: str
    organism: str

@app.post('/run_variants')
def run_variants(data: VariantInput):
    '''
    Run the pipeline, returning a readable error message instead of a bare 500
    when the input is unusable or the analysis fails.
    '''
    organism = ORGANISM_LABELS.get(data.organism, data.organism)

    try:
        return run_pipeline(data)

    except NoCodingVariantsError:
        message = (f'VEP could not annotate any of the submitted variants on the {organism} genome. '
                   f'Please confirm that the Organism setting matches your coordinates.')
        return JSONResponse(status_code=422, content={'error': message})

    except Exception:
        logging.exception('Analysis failed')
        message = (f'The analysis could not be completed. Please confirm that the Organism setting ({organism}) '
                   f'matches your coordinates and that each line uses the format chromosome:start-end:ref/alt, '
                   f'optionally followed by a tab and a label.')
        return JSONResponse(status_code=500, content={'error': message})


def run_pipeline(data):

    print(data.variants)
    print(data.organism)

    variants = process_batch_query(data.variants)

    if variants.empty:
        message = 'No variants were recognized. Enter one variant per line in the format chromosome:start-end:ref/alt, optionally followed by a tab and a label.'
        return JSONResponse(status_code=422, content={'error': message})

    # Labels are only displayed with the submitted variants; the analysis never sees them
    labels = (variants[variants['Label'] != ''].groupby('Name', sort=False)['Label']
              .agg(lambda s: '; '.join(dict.fromkeys(s))).to_dict())
    variants = variants.drop(columns='Label')

    import pandas as pd
    allele_coverage_df = pd.DataFrame()

    # Stay empty when none of the variants has a coding consequence
    output_gene_df, output_prt_df, output_phenotype_df, score_df = (pd.DataFrame() for _ in range(4))
    phenotype_similarity_df = pd.DataFrame()   # human queries only

    if data.organism == 'human':

        query_gene_df, query_prt_df, input_genes, unevaluated_df = hvar_query(variants)
        if not query_prt_df.empty:
            output_gene_df, output_prt_df, output_phenotype_df, input_genes, allele_coverage_df = mvar_fetch(input_genes)
            score_df = hvar_query_score(query_prt_df, output_prt_df, input_genes, mouse_phenotype_df=output_phenotype_df)
            output_phenotype_df = sort_phenotypes_by_score(output_phenotype_df, score_df)
            # One row per evaluated allele, with the term mapping behind each phenotype score
            if not score_df.empty:
                phenotype_similarity_df = phenotype_similarity(
                    query_prt_df, score_df[['Input Name', 'AlleleID']], output_phenotype_df, detail=True)

    if data.organism == 'mouse':

        query_gene_df, query_prt_df, input_genes, unevaluated_df, allele_coverage_df = mvar_query(variants)
        if not query_prt_df.empty:
            output_gene_df, output_prt_df, output_phenotype_df, input_genes = hvar_fetch(query_prt_df, input_genes)
            score_df = mvar_query_score(query_prt_df, output_prt_df, input_genes)
            output_phenotype_df = sort_phenotypes_by_score(output_phenotype_df, score_df)

    results = {'query_genes': query_gene_df.to_dict(orient='records'),
            'query_variants': query_prt_df.to_dict(orient='records'),
            'output_genes': output_gene_df.to_dict(orient='records'),
            'output_variants': output_prt_df.to_dict(orient='records'),
            'scores': score_df.to_dict(orient='records'),
            'phenotypes': output_phenotype_df.to_dict(orient='records'),
            'gene_mapping': input_genes.to_dict(orient='records'),
            'allele_coverage': allele_coverage_df.to_dict(orient='records'),
            'unevaluated_variants': unevaluated_df.to_dict(orient='records'),
            'phenotype_similarity': phenotype_similarity_df.to_dict(orient='records')}

    # Names for the HP term lists (HP_IDs_variant, HP_IDs_gene), shown beside each ID
    from app.hpoa import hp_names
    hp_ids = {hp for table in results.values() for row in table for col, value in row.items()
              if col.startswith('HP_IDs') and isinstance(value, list) for hp in value}
    results['hp_names'] = hp_names(hp_ids)
    results['labels'] = labels

    safe_result = jsonable_encoder(results)

    import math
    def _scrub(x):
        if isinstance(x, float) and math.isnan(x):
            return None
        if isinstance(x, list):
            return [_scrub(v) for v in x]
        if isinstance(x, dict):
            return {k: _scrub(v) for k, v in x.items()}
        return x

    safe_result = _scrub(safe_result)

    return JSONResponse(content=safe_result)

