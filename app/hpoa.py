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

from app.data_util import count_header_lines

with open('app/config.yaml') as f:
    _config = yaml.safe_load(f)

# Lazy singletons
_DISEASE_TO_HP = None       # {OMIM:/ORPHA:/DECIPHER: id -> set of HP:... strings}
_GENE_TO_HP = None          # {gene_symbol -> set of HP:... strings}
_DISEASE_NAMES = None     # {OMIM:/ORPHA:/DECIPHER: id -> disease name string}
_MONDO_TO_XREFS = None      # {MONDO:id -> set of other disease IDs}
_HP_NAMES = None            # {HP:id (and alt_id) -> term name}
_HP_CURRENT = None          # {retired or obsolete HP:id -> the current term's HP:id}


def _load_hp_obo():
    '''
    Read hp.obo once for two things:
      - term names, with each alt_id given the name of the term it was merged into
      - the current term for retired IDs (alt_id) and obsolete terms (replaced_by)
    Falls back to the names in genes_to_phenotype.txt, and no ID updates, when
    hp.obo hasn't been downloaded (python app/data_refresh.py fetches it).
    '''
    global _HP_NAMES, _HP_CURRENT
    from pathlib import Path
    names, current = {}, {}
    path = Path(_config['paths']['hp_obo'])
    if path.exists():
        term = {}

        def _finish(term):
            tid = term.get('id')
            if not tid:
                return
            if term.get('name'):
                names[tid] = term['name']
            for alt in term['alt']:
                current[alt] = tid
            if term.get('obsolete') and term.get('replaced_by'):
                current[tid] = term['replaced_by']

        with open(path, encoding='utf-8') as f:
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('['):
                    _finish(term)
                    term = {'alt': []} if line == '[Term]' else {}
                elif not term:
                    continue
                elif line.startswith('id: HP:'):
                    term['id'] = line[4:]
                elif line.startswith('name: '):
                    term['name'] = line[6:]
                elif line.startswith('alt_id: HP:'):
                    term['alt'].append(line[8:])
                elif line == 'is_obsolete: true':
                    term['obsolete'] = True
                elif line.startswith('replaced_by: HP:'):
                    term['replaced_by'] = line[len('replaced_by: '):]
        _finish(term)
        # A retired ID shows the name of the term it was merged into
        for old, new in current.items():
            if old not in names and new in names:
                names[old] = names[new]
        logging.info(f'hp.obo: {len(names)} HP term names, {len(current)} retired or obsolete IDs')
    else:
        df = pd.read_csv(_config['paths']['genes_to_phenotype'], sep='\t', comment='#', dtype=str,
                         usecols=['hpo_id', 'hpo_name']).dropna().drop_duplicates('hpo_id')
        names = dict(zip(df['hpo_id'], df['hpo_name']))
        logging.info(f'hp.obo not found; {len(names)} HP term names from genes_to_phenotype.txt')
    _HP_NAMES, _HP_CURRENT = names, current


def hp_names(hp_ids):
    '''{HP ID: term name} for the given IDs; IDs without a known name are left out.'''
    if _HP_NAMES is None:
        _load_hp_obo()
    return {h: _HP_NAMES[h] for h in hp_ids if h in _HP_NAMES}


def current_hp_ids(hp_ids):
    '''
    The given HP IDs with retired and obsolete ones replaced by their current
    term, without duplicates, sorted. ClinVar still cites some old IDs, e.g.
    HP:0000624 and HP:0007930, which are both now Epicanthus (HP:0000286).
    '''
    if _HP_CURRENT is None:
        _load_hp_obo()
    out = set()
    for h in hp_ids:
        seen = set()
        # Follow the chain in case a replacement was itself later retired
        while h in _HP_CURRENT and h not in seen:
            seen.add(h)
            h = _HP_CURRENT[h]
        out.add(h)
    return sorted(out)


def _build_mondo_xrefs():
    '''Invert mondo_xref_map (xref -> MONDO) into MONDO -> set of xrefs.'''
    with open(_config['paths']['mondo_xref_map']) as f:
        xref = json.load(f)
    inv = defaultdict(set)
    for xref_id, mondo in xref.items():
        inv[mondo].add(xref_id)
    logging.info(f'MONDO->xref index: {len(inv)} MONDO IDs')
    return dict(inv)

def _build_hpoa_indexes():
    '''
    Load phenotype.hpoa once into a {disease_id: set(HP:...)} map and a
    {disease_id: disease name} map.
    '''
    path = _config['paths']['hpoa']
    logging.info(f'Loading HPO annotations: {path}')
    df = pd.read_csv(
        path, sep='\t', skiprows=count_header_lines(path), dtype=str,
    )

    names_df = df.dropna(subset=['database_id', 'disease_name']).drop_duplicates(subset=['database_id'])
    names = dict(zip(names_df['database_id'], names_df['disease_name']))

    # Skip NOT annotations (the disease explicitly does NOT have this phenotype)
    df = df[df['qualifier'].fillna('') != 'NOT']
    index = defaultdict(set)
    for did, hp in df[['database_id', 'hpo_id']].itertuples(index=False):
        if isinstance(did, str) and isinstance(hp, str) and hp.startswith('HP:'):
            index[did].add(hp)
    logging.info(f'phenotype.hpoa indexed {len(index)} diseases')
    return dict(index), names


def _ensure_hpoa():
    global _DISEASE_TO_HP, _DISEASE_NAMES
    if _DISEASE_TO_HP is None or _DISEASE_NAMES is None:
        _DISEASE_TO_HP, _DISEASE_NAMES = _build_hpoa_indexes()


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
    global _MONDO_TO_XREFS
    if _MONDO_TO_XREFS is None:
        _MONDO_TO_XREFS = _build_mondo_xrefs()
    _ensure_hpoa()

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

def lookup_name_for_disease(disease_id):
    '''
    Return the human-readable name for a disease ID (OMIM/ORPHA/DECIPHER)
    from phenotype.hpoa. Returns None if the ID isn't found.
    Accepts 'Orphanet:778' and normalizes to 'ORPHA:778' before lookup.
    '''
    _ensure_hpoa()
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
