'''
MP (Mammalian Phenotype) to HP (Human Phenotype Ontology) equivalency scoring.

Translates mouse MP terms to HP terms via the Monarch mapping-commons SSSOM file,
then computes an asymmetric coverage score: the weighted fraction of a mouse
allele's translated HP terms that appear in a human reference HP set.
'''
import yaml
import pandas as pd

from app.data_util import count_header_lines

with open('app/config.yaml') as f:
    _config = yaml.safe_load(f)

# Weight each SSSOM predicate by how equivalent the mapping is. The weights are
# an ordinal choice, not fitted values:
# - exactMatch: the MP and HP terms describe the same phenotype.
# - narrowMatch/broadMatch: one term is a more specific form of the other, so the
#   phenotype is the same at a different level of detail; direction doesn't
#   matter for coverage, so both get the same weight.
# - closeMatch: similar enough to substitute in some contexts, but not asserted
#   to be the same phenotype.
# - relatedMatch: associated phenotypes only, so the weakest evidence.
# Each weight is multiplied by the mapping's own SSSOM confidence.
PREDICATE_WEIGHTS = {
    'skos:exactMatch': 1.0,
    'skos:narrowMatch': 0.75,
    'skos:broadMatch': 0.75,
    'skos:closeMatch': 0.5,
    'skos:relatedMatch': 0.25,
}


def _load_mp_hp_map():
    '''Build {MP_ID: [(HP_ID, weight), ...]} from the SSSOM TSV.'''
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

    mp_to_hp = {}
    for _, row in df.iterrows():
        mp_to_hp.setdefault(row['subject_id'], []).append(
            (row['object_id'], row['weight'])
        )
    return mp_to_hp


MP_TO_HP = _load_mp_hp_map()


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
    '''
    Asymmetric phenotype score: what fraction of the mouse allele's translated HPs
    appear in the reference HP set.

    Returns a dict:
      - 'score': float 0.0-1.0 or None. |mouse_HP ∩ reference_HP| / |mouse_HP|
      - 'mapped_mp_count': int, how many of the mouse MP terms are in the SSSOM map
      - 'total_mp_count': int, total MP terms attached to the allele

    Unlike Jaccard (symmetric), this metric does not penalize the reference set
    for being large and comprehensive — which is the right choice when the
    reference is a gene's full phenotype landscape and the mouse is one specific
    allele's observed phenotypes.
    '''
    mouse_mp_terms = list(mouse_mp_terms or [])
    total = len(mouse_mp_terms)
    mapped = sum(1 for mp in mouse_mp_terms if mp in MP_TO_HP)

    mouse_hp_weights = mp_to_weighted_hp(mouse_mp_terms)
    reference_hp = set(reference_hp_terms or [])

    if not mouse_hp_weights or not reference_hp:
        return {'score': None, 'mapped_mp_count': mapped, 'total_mp_count': total}

    intersection_weight = sum(w for hp, w in mouse_hp_weights.items() if hp in reference_hp)
    mouse_total_weight = sum(mouse_hp_weights.values())
    score = round(intersection_weight / mouse_total_weight, 3) if mouse_total_weight > 0 else None
    return {'score': score, 'mapped_mp_count': mapped, 'total_mp_count': total}
