'''
Score ClinVar's human variants against mouse alleles, one gene at a time.

Runs the same pipeline as a human query in the web app, outside the server,
and writes one set of files per gene so a long run can be stopped and resumed.

Usage (from the repository root):
    python -m app.clinvar_batch --genes MECP2
    python -m app.clinvar_batch --all --forks 6

Output, under --out (default app/results/clinvar_batch):
    variants/<GENE>.parquet     ClinVar records submitted for the gene
    scores/<GENE>.parquet       every ClinVar variant / mouse allele transcript pair
    best/<GENE>.parquet         per ClinVar variant, the top Total Score and its allele(s)
    unevaluated/<GENE>.parquet  variants that can't be scored (UTR, non-coding, ...)
    summary.tsv                 one line per gene: counts and run time
'''
import argparse
import logging
import time
from pathlib import Path

import pandas as pd
import pyarrow as pa

from app import vep_util
from app.disease_util import clinvar_vcf
from app.gene_util import mus_alleles_df, homology_df
from app.main_util import hvar_query, mvar_fetch, hvar_query_score, NoCodingVariantsError
from app.vep_util import config


human_genes_df = pd.read_parquet(config['paths']['human_gene_pqt'])


def genes_with_mouse_models():
    '''Human genes whose mouse homolog has a disease-associated MGI allele with coordinates.'''
    alleles = mus_alleles_df[(mus_alleles_df['HasDiseaseAnnotations'] == 'yes')
                             & mus_alleles_df['StartPosition'].notna()]
    mouse_genes = set(alleles['AlleleAssociatedGeneSymbol'])
    return sorted(set(homology_df.loc[homology_df['MusGeneSymbol'].isin(mouse_genes), 'HumGeneSymbol'].dropna()))


def clinvar_variants_for_gene(gene):
    '''
    ClinVar VCF records assigned to the gene (GENEINFO), in the web app's input format.

    Returns:
        DataFrame. Name, Chromosome, Start, Stop, Ref, Alt and the ClinVar fields
        ClinVar AlleleID, Clinical Significance, Review Status and ClinVar Consequence.
    '''
    gene_rows = human_genes_df[human_genes_df['Name'] == gene]
    if gene_rows.empty:
        return pd.DataFrame()
    g = gene_rows.iloc[0]

    rows = []
    for record in clinvar_vcf(f"{g['seqid']}:{g['start']}-{g['end']}"):
        geneinfo = record.INFO.get('GENEINFO') or ''
        if gene not in [entry.split(':')[0] for entry in geneinfo.split('|')]:
            continue
        ref = record.REF
        for alt in record.ALT:
            # Symbolic and missing alleles can't be annotated from coordinates
            if not alt or set(ref + alt) - set('ACGT'):
                continue
            start = record.POS
            stop = start + len(ref) - 1
            rows.append({'Name': f'{record.CHROM}:{start}-{stop}:{ref}/{alt}',
                         'Chromosome': str(record.CHROM), 'Start': start, 'Stop': stop,
                         'Ref': ref, 'Alt': alt,
                         'ClinVar AlleleID': record.INFO.get('ALLELEID'),
                         'Clinical Significance': record.INFO.get('CLNSIG'),
                         'Review Status': record.INFO.get('CLNREVSTAT'),
                         'ClinVar Consequence': record.INFO.get('MC')})

    return pd.DataFrame(rows).drop_duplicates(subset='Name').reset_index(drop=True) if rows else pd.DataFrame()


def best_matches(scores):
    '''
    Per ClinVar variant, the top Total Score and every allele that reaches it.
    Ties are all listed: there is no single best match when scores are equal.
    '''
    if scores.empty:
        return pd.DataFrame()
    scores = scores.assign(**{'Total Score': pd.to_numeric(scores['Total Score'])})
    top = scores[scores['Total Score'] == scores.groupby('Input Name')['Total Score'].transform('max')]
    return (top.groupby('Input Name')
            .agg(**{'Scoring Class': ('Scoring Class', 'first'),
                    'Best Total Score': ('Total Score', 'first'),
                    'Best AlleleID(s)': ('AlleleID', lambda s: '; '.join(sorted(set(s)))),
                    'Best Allele Symbol(s)': ('Allele Symbol', lambda s: '; '.join(sorted(set(s))))})
            .reset_index())


SUMMARY_COLS = ['Gene', 'Status', 'ClinVar Variants', 'Scored Variants', 'Unevaluated Variants',
                'Mouse Alleles', 'Score Rows', 'Seconds']


def _write(df, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df = df.copy()
    # Columns mixing value types (e.g. text and numbers) are stored as text so parquet accepts them
    for col in df.columns:
        try:
            pa.array(df[col])
        except (pa.ArrowInvalid, pa.ArrowTypeError):
            df[col] = df[col].map(lambda v: None if v is None else str(v))
    df.to_parquet(path, index=False)


def run_gene(gene, out_dir):
    '''Score one gene's ClinVar variants and write its files. Returns a summary dict.'''
    # Clear results from an earlier run so a rerun can't leave stale files behind
    for old in out_dir.glob(f'*/{gene}.parquet'):
        old.unlink()

    summary = {'Gene': gene, 'Status': 'ok', 'ClinVar Variants': 0, 'Scored Variants': 0,
               'Unevaluated Variants': 0, 'Mouse Alleles': 0, 'Score Rows': 0}

    variants = clinvar_variants_for_gene(gene)
    summary['ClinVar Variants'] = len(variants)
    if variants.empty:
        summary['Status'] = 'no ClinVar variants'
        return summary

    clinvar_cols = ['Name', 'ClinVar AlleleID', 'Clinical Significance', 'Review Status', 'ClinVar Consequence']
    clinvar_info = variants[clinvar_cols].rename(columns={'Name': 'Input Name'})
    _write(variants, out_dir / 'variants' / f'{gene}.parquet')

    try:
        _, protein_df, input_genes, unevaluated_df = hvar_query(variants[['Name', 'Chromosome', 'Start', 'Stop', 'Ref', 'Alt']].copy())
    except NoCodingVariantsError:
        summary['Status'] = 'VEP annotated no variants'
        return summary

    if not unevaluated_df.empty:
        _write(unevaluated_df, out_dir / 'unevaluated' / f'{gene}.parquet')
    summary['Unevaluated Variants'] = len(unevaluated_df)
    summary['Scored Variants'] = protein_df['Name'].nunique() if not protein_df.empty else 0

    scores = pd.DataFrame()
    if not protein_df.empty:
        _, mouse_prt_df, mouse_phenotype_df, input_genes, _ = mvar_fetch(input_genes)
        summary['Mouse Alleles'] = mouse_prt_df['AlleleID'].nunique() if not mouse_prt_df.empty else 0
        scores = hvar_query_score(protein_df, mouse_prt_df, input_genes, mouse_phenotype_df=mouse_phenotype_df)

    if scores.empty:
        summary['Status'] = 'no mouse alleles to score against'
    else:
        scores = clinvar_info.merge(scores, on='Input Name', how='right')
        _write(scores, out_dir / 'scores' / f'{gene}.parquet')
        _write(clinvar_info.merge(best_matches(scores), on='Input Name', how='inner'),
               out_dir / 'best' / f'{gene}.parquet')
        summary['Score Rows'] = len(scores)

    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    which = parser.add_mutually_exclusive_group(required=True)
    which.add_argument('--genes', nargs='+', help='Human gene symbols to run')
    which.add_argument('--all', action='store_true', help='Every human gene with a disease-associated mouse allele')
    parser.add_argument('--out', default='app/results/clinvar_batch', help='Output directory')
    parser.add_argument('--forks', type=int, default=None, help='VEP worker processes')
    parser.add_argument('--rerun', action='store_true', help='Rerun genes that already have results')
    args = parser.parse_args()

    vep_util.vep_forks = args.forks
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / 'summary.tsv'

    genes = args.genes if args.genes else genes_with_mouse_models()

    summary_df = (pd.read_csv(summary_path, sep='\t') if summary_path.exists()
                  else pd.DataFrame(columns=SUMMARY_COLS))

    # Genes that failed with an error are run again; finished genes are skipped
    done = set() if args.rerun else set(summary_df.loc[~summary_df['Status'].str.startswith('error'), 'Gene'])

    todo = [g for g in genes if g not in done]
    print(f'{len(todo)} gene(s) to run, {len(genes) - len(todo)} already done')

    for i, gene in enumerate(todo, 1):
        s = time.time()
        try:
            summary = run_gene(gene, out_dir)
        except Exception as e:
            logging.exception(f'ClinVar batch failed for {gene}')
            summary = {'Gene': gene, 'Status': f'error: {e}'}
        summary['Seconds'] = round(time.time() - s, 1)

        # One line per gene: a rerun replaces the gene's earlier line
        summary_df = pd.concat([summary_df[summary_df['Gene'] != gene],
                                pd.DataFrame([summary]).reindex(columns=SUMMARY_COLS)], ignore_index=True)
        summary_df.astype({c: 'Int64' for c in SUMMARY_COLS[2:-1]}).to_csv(summary_path, sep='\t', index=False)
        print(f"[{i}/{len(todo)}] {gene}: {summary.get('Status')} "
              f"({summary.get('ClinVar Variants', 0)} variants, {summary.get('Score Rows', 0)} score rows, "
              f"{summary.get('Seconds', '?')} s)")


if __name__ == '__main__':
    main()
