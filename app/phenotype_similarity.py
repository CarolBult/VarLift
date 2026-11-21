'''
MP (Mammalian Phenotype) to HP (Human Phenotype Ontology) equivalency scoring.

Translates mouse MP terms to HP terms via the Monarch mapping-commons SSSOM file,
then computes an asymmetric coverage score: the weighted share of a mouse
allele's translated HP terms that a human reference HP set accounts for, where
each term counts by its semantic similarity to the closest reference term
(app.hpo_similarity; 1 for the identical term).
'''
import logging
import yaml
import pandas as pd

from app.data_util import count_header_lines

with open('app/config.yaml') as f:
    _config = yaml.safe_load(f)

# Weight each SSSOM predicate by how equivalent the mapping is. The weight sets how much a
# translated HP term counts in the score (numerator and denominator alike), so a loose
# translation can neither earn nor lose as much as an exact one:
#   - exactMatch (1.0): the MP and HP terms describe the same phenotype.
#   - narrowMatch / broadMatch (0.75): one term is a more specific form of the other. The
#     phenotype is the same kind but at a different level of detail; neither direction is
#     favoured, since the score cares whether the phenotype is present, not how specific.
#   - closeMatch (0.5): similar enough to use in some contexts but, by the SKOS definition,
#     not interchangeable.
#   - relatedMatch (0.25): associated terms (for example, a structure and its function) that
#     are not the same phenotype; kept at a low weight rather than dropped.
# The steps are an ordinal judgement of mapping strength, not fitted values. Each weight is
# multiplied by the mapping's SSSOM confidence (1.0 when it gives none), and predicates not
# listed here (none in the current file) are left out.
PREDICATE_WEIGHTS = {
    'skos:exactMatch': 1.0,
    'skos:narrowMatch': 0.75,
    'skos:broadMatch': 0.75,
    'skos:closeMatch': 0.5,
    'skos:relatedMatch': 0.25,
}


# The SSSOM file still targets some HP IDs that HPO has since retired. A mapping to a
# retired ID can never match current HPO annotations, so it silently lowers scores.
# Each one was reviewed (2026-11-10, HPO release 2026-09-01) against the term HPO
# now points it to:
#   - replaced: the current term keeps the meaning of the retired one
#   - dropped:  the current term changes the meaning (it loses part of the phenotype or is
#               narrower), so the mapping is left out rather than redirected
# A retired ID not listed here (e.g. after a data refresh) keeps the old ID and is logged
# for review.
RETIRED_HP_REPLACED = {
    'HP:0007222',  # CNS demyelination -> CNS demyelination (HP:0007305)
    'HP:0030782',  # Abnormal circulating interleukin concentration -> same name (HP:0011117)
    'HP:0002275',  # Poor motor coordination -> Incoordination (HP:0002311)
    'HP:0001388',  # Joint laxity -> Joint hypermobility (HP:0001382)
    'HP:0002648',  # Abnormality of calvarial morphology -> Abnormal calvaria morphology (HP:0002683)
    'HP:0008572',  # External ear malformation -> Abnormal pinna morphology (HP:0000377)
}
RETIRED_HP_DROPPED = {
    'HP:0000368',  # Low-set, posteriorly rotated ears -> Posteriorly rotated ears: loses "low-set"
    'HP:0002744',  # Bilateral cleft lip and palate -> Bilateral cleft palate: loses the lip
    'HP:0001180',  # Hand oligodactyly -> Finger aplasia: narrower
    'HP:0008165',  # Decreased helper T cell proportion -> Decreased total CD4+ T cell proportion
    'HP:0010960',  # Bronchopulmonary sequestration -> Pulmonary sequestration (mapped from lung cysts)
}


def _update_retired_hp(df):
    '''Apply the reviewed replacements and drops above to the mapping's HP IDs.'''
    from app.hpoa import current_hp_ids
    current = {h: current_hp_ids([h])[0] for h in df['object_id'].dropna().unique()}
    retired = {h for h, c in current.items() if c != h}

    df = df[~df['object_id'].isin(RETIRED_HP_DROPPED)].copy()
    replace = {h: current[h] for h in retired & RETIRED_HP_REPLACED}
    df['object_id'] = df['object_id'].replace(replace)

    unreviewed = sorted(retired - RETIRED_HP_REPLACED - RETIRED_HP_DROPPED)
    if unreviewed:
        logging.warning(f'MP-HP mapping uses retired HP IDs not yet reviewed (left unchanged): {unreviewed}')
    return df


def _load_mp_hp_map():
    '''Build {MP_ID: [(HP_ID, weight), ...]} and {(MP_ID, HP_ID): predicate} from the SSSOM TSV.'''
    df = pd.read_csv(
        _config['paths']['mp_hp_mapping'],
        sep='\t',
        skiprows=count_header_lines(_config['paths']['mp_hp_mapping']),
        dtype=str,
    )
    df['confidence'] = pd.to_numeric(df['confidence'], errors='coerce').fillna(1.0)
    df['pred_weight'] = df['predicate_id'].map(PREDICATE_WEIGHTS).fillna(0.0)
    df['weight'] = df['pred_weight'] * df['confidence']
    df = df[df['weight'] > 0]
    df = _update_retired_hp(df)

    mp_to_hp = {}
    predicates = {}
    for _, row in df.iterrows():
        pairs = mp_to_hp.setdefault(row['subject_id'], [])
        key = (row['subject_id'], row['object_id'])
        # A replacement can land on an HP the MP already maps to; keep the strongest weight
        for i, (hp, w) in enumerate(pairs):
            if hp == row['object_id']:
                if row['weight'] > w:
                    pairs[i] = (hp, row['weight'])
                    predicates[key] = row['predicate_id']
                break
        else:
            pairs.append((row['object_id'], row['weight']))
            predicates[key] = row['predicate_id']
    return mp_to_hp, predicates


# {MP_ID: [(HP_ID, weight), ...]} and {(MP_ID, HP_ID): SKOS predicate} for display
MP_TO_HP, MP_HP_PREDICATE = _load_mp_hp_map()


# MP terms left out of phenotype scoring:
# - the top-level system headers (children of MP:0000001, e.g. "nervous system phenotype"),
#   which MGI's allele file lists alongside the specific terms. They translate to equally
#   broad HP terms (e.g. "Abnormality of the nervous system") that disease annotations
#   almost never use, so they could only lower a score.
# - "no abnormal phenotype detected", which is not a phenotype.
NO_ABNORMAL_PHENOTYPE = 'MP:0002169'


def _load_unscored_mp():
    from pathlib import Path
    unscored = {'MP:0000001', NO_ABNORMAL_PHENOTYPE}
    path = Path(_config['paths']['mp_obo'])
    if not path.exists():
        logging.warning('mp.obo not found (python app/data_refresh.py fetches it); '
                        'top-level MP terms will be scored')
        return unscored
    term = None
    with open(path, encoding='utf-8') as f:
        for line in f:
            if line.startswith('[Term]'):
                term = None
            elif line.startswith('id: MP:'):
                term = line[4:].strip()
            elif line.startswith('is_a: MP:0000001 ') or line.strip() == 'is_a: MP:0000001':
                unscored.add(term)
    return unscored


UNSCORED_MP = _load_unscored_mp()


def scored_mp_terms(mp_terms):
    '''Split an allele's MP terms into (scored, not scored); see UNSCORED_MP.'''
    mp_terms = set(mp_terms or [])
    return mp_terms - UNSCORED_MP, mp_terms & UNSCORED_MP


def mp_to_weighted_hp(mp_terms):
    '''
    Translate a set of MP IDs to {HP_ID: max_weight}.
    When several MPs map to the same HP, keep the strongest weight.
    '''
    hp_weights = {}
    for mp in mp_terms or []:
        for hp, w in MP_TO_HP.get(mp, []):
            if w > hp_weights.get(hp, 0):
                hp_weights[hp] = w
    return hp_weights


def phenotype_coverage_detailed(mouse_mp_terms, reference_hp_terms):
    '''Score only; see phenotype_mapping_detail.'''
    d = phenotype_mapping_detail(mouse_mp_terms, reference_hp_terms)
    return {k: d[k] for k in ('score', 'mapped_mp_count', 'total_mp_count')}


def phenotype_mapping_detail(mouse_mp_terms, reference_hp_terms, mp_names=None):
    '''
    Asymmetric phenotype score with the term mapping behind it.

    1. Each mouse MP term is translated to HP terms through the SSSOM mapping; each
       translated HP term counts once, with the strongest translation weight w any of
       the allele's MP terms gives it.
    2. Each translated HP term is compared with its closest term in the human
       reference set by Lin semantic similarity s (app.hpo_similarity): 1 when the
       identical term is there, less for a parent, child or more distant relative.
    3. score = sum(w * s) / sum(w): how much of the allele's phenotype the human
       reference accounts for, from 0 to 1. A large reference set is not penalized.

    Returns a dict with score, mapped_mp_count, total_mp_count, matched_weight
    (sum of w * s), total_weight (sum of w), and
      - 'pairs': one per MP -> HP translation {mp, mp_name, hp, relation, weight,
        human_hp, similarity, link, mica, counted}; human_hp is the closest reference
        term, link says how the two are related, and counted is False when another MP
        gives the same HP a stronger weight (the HP is counted once)
      - 'unmapped': [{mp, mp_name}] for MP terms with no HP translation
    '''
    from app.hpo_similarity import closest, relation as link_of

    mp_names = mp_names or {}
    mouse_mp_terms = sorted(set(mouse_mp_terms or []))
    reference = frozenset(reference_hp_terms or [])
    hp_weights = mp_to_weighted_hp(mouse_mp_terms)
    best = {hp: closest(hp, reference) for hp in hp_weights} if reference else {}

    pairs, unmapped, counted_hp = [], [], set()
    for mp in mouse_mp_terms:
        targets = MP_TO_HP.get(mp, [])
        if not targets:
            unmapped.append({'mp': mp, 'mp_name': mp_names.get(mp)})
            continue
        for hp, w in targets:
            # Each HP counts once: the first MP to reach its strongest weight carries it
            counted = w == hp_weights[hp] and hp not in counted_hp
            if counted:
                counted_hp.add(hp)
            sim, human_hp, mica = best.get(hp, (0.0, None, None))
            pairs.append({'mp': mp, 'mp_name': mp_names.get(mp), 'hp': hp,
                          'relation': MP_HP_PREDICATE.get((mp, hp), '').replace('skos:', ''),
                          'weight': round(w, 3), 'human_hp': human_hp, 'similarity': round(sim, 3),
                          'link': link_of(hp, human_hp, mica), 'mica': mica, 'counted': counted,
                          'matched': hp in reference})

    total_weight = sum(hp_weights.values())
    matched_weight = sum(w * best[hp][0] for hp, w in hp_weights.items()) if best else 0.0
    score = (round(matched_weight / total_weight, 3)
             if hp_weights and reference and total_weight > 0 else None)
    return {'score': score, 'mapped_mp_count': len(mouse_mp_terms) - len(unmapped),
            'total_mp_count': len(mouse_mp_terms), 'matched_weight': round(matched_weight, 3),
            'total_weight': round(total_weight, 3), 'pairs': pairs, 'unmapped': unmapped}
