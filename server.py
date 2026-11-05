from fastapi import FastAPI
from pydantic import BaseModel
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.encoders import jsonable_encoder

import math
import logging

from app.main_util import hvar_query, mvar_fetch, hvar_query_score, mvar_query, hvar_fetch, mvar_query_score, NoCodingVariantsError
from app.processing_util import process_batch_query

ORGANISM_LABELS = {'human': 'Human (GRCh38)', 'mouse': 'Mouse (GRCm39)'}

def sort_phenotypes_by_score(phenotype_df, score_df):
    '''
    Reorder phenotype_df so alleles with the highest Total Score appear first.
    Uses each allele's best (max) transcript score. Alleles with no score go to the bottom.
    '''
    if phenotype_df.empty or score_df.empty or 'AlleleID' not in phenotype_df.columns:
        return phenotype_df

    allele_order = (
        score_df.groupby('AlleleID')['Total Score']
        .max()
        .sort_values(ascending=False)
        .index.tolist()
    )
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


def scrub_nan(x):
    '''
    Replace NaN floats with None throughout nested lists and dicts, since
    JSON has no NaN.
    '''
    if isinstance(x, float) and math.isnan(x):
        return None
    if isinstance(x, list):
        return [scrub_nan(v) for v in x]
    if isinstance(x, dict):
        return {k: scrub_nan(v) for k, v in x.items()}
    return x


app = FastAPI()

app.mount('/static', StaticFiles(directory='static'), name='static')

@app.get('/')
def serve_dashboard():
    return FileResponse('templates/dashboard.html')

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
        message = (f'None of the submitted variants fall in a protein-coding gene on the {organism} genome. '
                   f'Please confirm that the Organism setting matches your coordinates.')
        return JSONResponse(status_code=422, content={'error': message})

    except Exception:
        logging.exception('Analysis failed')
        message = (f'The analysis could not be completed. Please confirm that the Organism setting ({organism}) '
                   f'matches your coordinates and that each line uses the format chromosome:start-end:ref/alt.')
        return JSONResponse(status_code=500, content={'error': message})


def run_pipeline(data):

    print(data.variants)
    print(data.organism)

    variants = process_batch_query(data.variants)

    if variants.empty:
        message = 'No variants were recognized. Enter one variant per line in the format chromosome:start-end:ref/alt.'
        return JSONResponse(status_code=422, content={'error': message})

    if data.organism == 'human':

        query_gene_df, query_prt_df, input_genes = hvar_query(variants)
        output_gene_df, output_prt_df, output_phenotype_df, input_genes = mvar_fetch(input_genes)
        score_df = hvar_query_score(query_prt_df, output_prt_df, input_genes, mouse_phenotype_df=output_phenotype_df)
        output_phenotype_df = sort_phenotypes_by_score(output_phenotype_df, score_df)

    if data.organism == 'mouse':

        query_gene_df, query_prt_df, input_genes = mvar_query(variants)
        output_gene_df, output_prt_df, output_phenotype_df, input_genes = hvar_fetch(query_prt_df, input_genes)
        score_df = mvar_query_score(query_prt_df, output_prt_df, input_genes)
        output_phenotype_df = sort_phenotypes_by_score(output_phenotype_df, score_df)

    results = {'query_genes': query_gene_df.to_dict(orient='records'),
            'query_variants': query_prt_df.to_dict(orient='records'),
            'output_genes': output_gene_df.to_dict(orient='records'),
            'output_variants': output_prt_df.to_dict(orient='records'),
            'scores': score_df.to_dict(orient='records'),
            'phenotypes': output_phenotype_df.to_dict(orient='records'),
            'gene_mapping': input_genes.to_dict(orient='records')}

    safe_result = scrub_nan(jsonable_encoder(results))

    return JSONResponse(content=safe_result)
