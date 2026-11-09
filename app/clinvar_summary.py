'''
Lookups against ClinVar's variant_summary, stored as a position-sorted parquet
(see data_util.write_variant_summary_parquet).

- HP terms and disease IDs for human variants, by (Chromosome, Start) and the
  ALLELEIDs of the ClinVar VCF records with the variant's exact ref/alt. Only
  the variant's own records count: another change at the same position (R206H
  for R206L) says nothing about the variant. A much richer HP source than CLNDISDB.
- All variants for a set of genes, for finding homologous human variants.

Each lookup reads only the parquet row groups that can hold matching rows,
so nothing is held in memory between searches beyond a small result cache.
'''
import re
import yaml
from functools import lru_cache

import pyarrow.parquet as pq

with open('app/config.yaml') as f:
    _config = yaml.safe_load(f)

_PATH = _config['paths']['ncbi_var_summary_pqt']

_HP_PATTERN = re.compile(r'Human Phenotype Ontology:(HP:\d+)')
_DISEASE_PATTERN = re.compile(r'(MONDO:MONDO:\d+|OMIM:\d+|Orphanet:\d+)')


@lru_cache(maxsize=4096)
def _phenotype_ids_at(chromosome, start, allele_ids):
    '''
    Collect HP IDs and disease IDs (MONDO/OMIM/Orphanet) from PhenotypeIDS
    for the ClinVar records at (chromosome, start) with the given AlleleIDs.
    '''
    if not allele_ids:
        return frozenset(), frozenset()

    table = pq.read_table(_PATH, columns=['PhenotypeIDS'],
                          filters=[('Chromosome', '=', str(chromosome)), ('Start', '=', int(start)),
                                   ('AlleleID', 'in', list(allele_ids))])

    hps = set()
    diseases = set()
    for pids in table.column('PhenotypeIDS').to_pylist():
        if pids is None:
            continue

        hps.update(_HP_PATTERN.findall(pids))

        for m in _DISEASE_PATTERN.findall(pids):
            # Normalize 'MONDO:MONDO:...' -> 'MONDO:...'
            if m.startswith('MONDO:MONDO:'):
                m = m[6:]
            diseases.add(m)

    return frozenset(hps), frozenset(diseases)


def lookup_hp_ids(chromosome, start, allele_ids):
    '''
    Return sorted list of HP IDs for a variant at (chromosome, start), from its
    ClinVar records (allele_ids, from the VCF). Returns an empty list if the
    variant isn't in ClinVar.
    '''
    hps, _ = _phenotype_ids_at(str(chromosome), int(start), tuple(sorted(allele_ids)))
    return sorted(hps)


def lookup_disease_ids(chromosome, start, allele_ids):
    '''
    Return sorted list of disease IDs (MONDO, OMIM, Orphanet) extracted from
    variant_summary's PhenotypeIDS for the variant at (chromosome, start), from
    its ClinVar records (allele_ids, from the VCF).
    Used to chain through phenotype.hpoa for HP lookup.
    '''
    _, diseases = _phenotype_ids_at(str(chromosome), int(start), tuple(sorted(allele_ids)))
    return sorted(diseases)


def variants_for_genes(genes):
    '''
    Return all variant_summary rows for the given gene symbols, in the
    original file order.
    '''
    genes = [g for g in genes if isinstance(g, str)]
    parquet = pq.ParquetFile(_PATH)

    # The file is sorted by position, so a gene's rows sit in a few row groups.
    # Find them from the GeneSymbol column alone, then read only those groups.
    row_groups = [i for i in range(parquet.num_row_groups)
                  if genes and parquet.read_row_group(i, columns=['GeneSymbol']).column(0).to_pandas().isin(genes).any()]

    df = parquet.read_row_groups(row_groups).to_pandas() if row_groups else parquet.schema_arrow.empty_table().to_pandas()
    df = df[df['GeneSymbol'].isin(genes)]
    return df.sort_values('_order').drop(columns='_order').reset_index(drop=True)
