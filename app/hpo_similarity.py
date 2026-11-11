'''
Semantic similarity between HPO terms, for phenotype scoring.

Two HPO terms rarely match exactly: a mouse allele's "Vertebral fusion" and a
disease's "Fused cervical vertebrae" are parent and child, and an identical-term
test scores them as unrelated. As in Monarch's Phenodigm, terms are compared
through the HPO hierarchy and information content (IC):

  IC(t)   = -log(share of HPO-annotated diseases annotated with t or a descendant)
  MICA    = the most informative common ancestor of two terms (highest IC)
  Lin     = 2 * IC(MICA) / (IC(a) + IC(b))

Lin is 1 for identical terms and high for a parent and child. Terms whose most
informative shared ancestor is a top-level category or an organ system (e.g.
"Abnormality of the skeletal system") score 0: sharing an organ system is not
phenotypic similarity.
Phenodigm combines Resnik IC with Jaccard instead; Lin is used here because it
stays on the 0-1 scale of VarLift's coverage score.

Data: hp.obo (is_a hierarchy) and phenotype.hpoa (the annotation corpus for IC).
'''
import logging
import math
from collections import defaultdict
from functools import lru_cache

import yaml

with open('app/config.yaml') as f:
    _config = yaml.safe_load(f)

_PARENTS = None     # {HP:id: [parent HP:ids]}
_IC = None          # {HP:id: information content}


def _load():
    global _PARENTS, _IC
    from app.hpoa import _ensure_hpoa
    import app.hpoa as hpoa

    parents, term, obsolete = defaultdict(list), None, False
    with open(_config['paths']['hp_obo'], encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('['):
                term, obsolete = None, False
            elif line.startswith('id: HP:'):
                term = line[4:]
            elif line == 'is_obsolete: true':
                obsolete = True
                parents.pop(term, None)
            elif line.startswith('is_a: HP:') and term and not obsolete:
                parents[term].append(line[6:16])
    _PARENTS = dict(parents)

    # IC from the disease annotation corpus: each disease counts once for every
    # term it is annotated with and every ancestor of those terms
    _ensure_hpoa()
    counts = defaultdict(int)
    diseases = hpoa._DISEASE_TO_HP
    for hps in diseases.values():
        for t in set().union(*(ancestors(h) for h in hps)) if hps else ():
            counts[t] += 1
    n = len(diseases)
    # A term no disease uses is as specific as a term used once
    floor = -math.log(1 / n)
    _IC = defaultdict(lambda: floor, {t: -math.log(c / n) for t, c in counts.items()})
    logging.info(f'HPO similarity: {len(_PARENTS)} terms, IC from {n} diseases')


@lru_cache(maxsize=None)
def ancestors(term):
    '''The term and all its is_a ancestors.'''
    if _PARENTS is None:
        _load()
    out = {term}
    for p in _PARENTS.get(term, ()):
        out |= ancestors(p)
    return frozenset(out)


@lru_cache(maxsize=1)
def _too_general():
    '''
    Ancestors too general to show similarity: HP:0000001 (All), the top-level categories
    (Phenotypic abnormality, Clinical modifier, Mode of inheritance, ...) and their direct
    children (organ systems such as Abnormality of the skeletal system, Clinical course).
    '''
    if _PARENTS is None:
        _load()
    top = {'HP:0000001'} | {t for t, ps in _PARENTS.items() if 'HP:0000001' in ps}
    return frozenset(top | {t for t, ps in _PARENTS.items() if set(ps) & top})


def ic(term):
    if _IC is None:
        _load()
    return _IC[term]


@lru_cache(maxsize=None)
def lin(a, b):
    '''(Lin similarity, MICA term) for two HPO terms.'''
    if a == b:
        return 1.0, a
    common = ancestors(a) & ancestors(b)
    if not common:
        return 0.0, None
    mica = max(common, key=ic)
    # Sharing only a top-level category or its direct child (Phenotypic abnormality, an organ
    # system such as Abnormality of the skeletal system, Clinical course, ...) is not phenotypic
    # similarity, e.g. Neonatal death and Congenital onset share only Clinical course
    if mica in _too_general():
        return 0.0, None
    denom = ic(a) + ic(b)
    return (2 * ic(mica) / denom if denom > 0 else 0.0), mica


def closest(term, reference):
    '''
    The reference term most similar to term: (similarity, reference term, MICA).
    reference: frozenset of HPO IDs. Returns (0.0, None, None) for an empty reference.
    '''
    return _closest(term, reference if isinstance(reference, frozenset) else frozenset(reference))


@lru_cache(maxsize=200_000)
def _closest(term, reference):
    if term in reference:
        return 1.0, term, term
    best = (0.0, None, None)
    for r in reference:
        s, mica = lin(term, r)
        if s > best[0]:
            best = (s, r, mica)
    return best


def relation(mouse_hp, human_hp, mica):
    '''How a translated mouse term relates to its closest human term, in words.'''
    if human_hp is None:
        return 'none'
    if mouse_hp == human_hp:
        return 'identical'
    if human_hp in ancestors(mouse_hp):
        return 'mouse term more specific'
    if mouse_hp in ancestors(human_hp):
        return 'mouse term more general'
    return 'shared ancestor'
