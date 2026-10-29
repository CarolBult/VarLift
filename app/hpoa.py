'''
Human Phenotype Ontology annotation lookups.

Loads two HPO files lazily at first call:
- phenotype.hpoa: disease (OMIM/ORPHA/DECIPHER) -> HP terms, for variant-level
  HP lookup via the MONDO -> OMIM/Orphanet chain.
- genes_to_phenotype.txt: gene symbol -> HP terms, for gene-level HP lookup.

Also builds a reverse index of mondo_xref_map.json so a MONDO ID can be walked
to its OMIM/Orphanet/other equivalents.
'''
import json
import yaml
import logging
import pandas as pd
from collections import defaultdict

with open('app/config.yaml') as f:
    _config = yaml.safe_load(f)

# Lazy singletons
_DISEASE_TO_HP = None       # {OMIM:/ORPHA:/DECIPHER: id -> set of HP:... strings}
_GENE_TO_HP = None          # {gene_symbol -> set of HP:... strings}
_DISEASE_NAMES = None     # {OMIM:/ORPHA:/DECIPHER: id -> disease name string}
_MONDO_TO_XREFS = None      # {MONDO:id -> set of other disease IDs}

def _build_mondo_xrefs():
    '''Invert mondo_xref_map (xref -> MONDO) into MONDO -> set of xrefs.'''
    with open(_config['paths']['mondo_xref_map']) as f:
        xref = json.load(f)
    inv = defaultdict(set)
    for xref_id, mondo in xref.items():
        inv[mondo].add(xref_id)
    logging.info(f'MONDO->xref index: {len(inv)} MONDO IDs')
    return dict(inv)

def _build_disease_to_hp():
    '''Load phenotype.hpoa into a {disease_id: set(HP:...)} map.'''
    path = _config['paths']['hpoa']
    logging.info(f'Loading HPO annotations: {path}')
    df = pd.read_csv(
        path, sep='\t', comment='#', dtype=str,
    )
    # Skip NOT annotations (the disease explicitly does NOT have this phenotype)
    df = df[df['qualifier'].fillna('') != 'NOT']
    index = defaultdict(set)
    for did, hp in df[['database_id', 'hpo_id']].itertuples(index=False):
        if isinstance(did, str) and isinstance(hp, str) and hp.startswith('HP:'):
            index[did].add(hp)
    logging.info(f'phenotype.hpoa indexed {len(index)} diseases')
    return dict(index)


def _build_gene_to_hp():
    '''Load genes_to_phenotype.txt into a {gene_symbol: set(HP:...)} map.'''
    path = _config['paths']['genes_to_phenotype']
    logging.info(f'Loading gene->HP annotations: {path}')
    df = pd.read_csv(
        path, sep='\t', comment='#', dtype=str,
    )
    index = defaultdict(set)
    for gene, hp in df[['gene_symbol', 'hpo_id']].itertuples(index=False):
        if isinstance(gene, str) and isinstance(hp, str) and hp.startswith('HP:'):
            index[gene].add(hp)
    logging.info(f'genes_to_phenotype indexed {len(index)} genes')
    return dict(index)


def _normalize_disease_id(did):
    '''
    Normalize a disease ID so it matches phenotype.hpoa's conventions.
    mondo_xref_map uses 'Orphanet:778'; hpoa uses 'ORPHA:778'.
    '''
    if did.startswith('Orphanet:'):
        return 'ORPHA:' + did.split(':', 1)[1]
    return did


def lookup_hps_for_mondo(mondo_id):
    '''
    Return a sorted list of HP IDs for a MONDO disease by walking:
      MONDO -> xref IDs (OMIM/ORPHA/DECIPHER) -> phenotype.hpoa -> HP terms.
    Returns an empty list if the MONDO isn't recognised or has no HP annotations.
    '''
    global _MONDO_TO_XREFS, _DISEASE_TO_HP
    if _MONDO_TO_XREFS is None:
        _MONDO_TO_XREFS = _build_mondo_xrefs()
    if _DISEASE_TO_HP is None:
        _DISEASE_TO_HP = _build_disease_to_hp()

    xrefs = _MONDO_TO_XREFS.get(mondo_id, set())
    hps = set()
    for x in xrefs:
        norm = _normalize_disease_id(x)
        hps.update(_DISEASE_TO_HP.get(norm, set()))
    return sorted(hps)


def lookup_hps_for_gene(gene_symbol):
    '''
    Return a sorted list of HP IDs for a gene from genes_to_phenotype.txt.
    Returns an empty list if the gene has no HP annotations.
    '''
    global _GENE_TO_HP
    if _GENE_TO_HP is None:
        _GENE_TO_HP = _build_gene_to_hp()
    return sorted(_GENE_TO_HP.get(gene_symbol, set()))

def _build_disease_names():
    '''Load phenotype.hpoa and build {disease_id: disease_name}.'''
    path = _config['paths']['hpoa']
    df = pd.read_csv(path, sep='\t', comment='#', dtype=str,
                     usecols=['database_id', 'disease_name'])
    df = df.dropna(subset=['database_id', 'disease_name']).drop_duplicates(subset=['database_id'])
    return dict(zip(df['database_id'], df['disease_name']))

def lookup_name_for_disease(disease_id):
    '''
    Return the human-readable name for a disease ID (OMIM/ORPHA/DECIPHER)
    from phenotype.hpoa. Returns None if the ID isn't found.
    Accepts 'Orphanet:778' and normalizes to 'ORPHA:778' before lookup.
    '''
    global _DISEASE_NAMES
    if _DISEASE_NAMES is None:
        _DISEASE_NAMES = _build_disease_names()
    return _DISEASE_NAMES.get(_normalize_disease_id(disease_id))


def lookup_name_for_mondo(mondo_id):
    '''
    Return a human-readable disease name for a MONDO ID by walking
    its xref list and asking phenotype.hpoa for the name under OMIM
    (preferred) or Orphanet (fallback).

    Returns None if no name is available.
    '''
    global _MONDO_TO_XREFS
    if _MONDO_TO_XREFS is None:
        _MONDO_TO_XREFS = _build_mondo_xrefs()

    xrefs = _MONDO_TO_XREFS.get(mondo_id, set())
    omim = sorted(x for x in xrefs if x.startswith('OMIM:'))
    orph = sorted(x for x in xrefs if x.startswith('Orphanet:'))

    for xref in omim + orph:
        name = lookup_name_for_disease(xref)
        if name:
            return name
    return None
