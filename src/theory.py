#!/usr/bin/env python3

import property as prop
import rcc5, checklist, workspace, simple, exemplar, typify

from util import log
from checklist import *
from workspace import *
from rcc5 import rcc5_symbol

import specimen
import exemplar
import ranks
from specimen import same_specimens, get_exemplar, same_type_ufs, \
  sid_to_opposite_record
from exemplar import find_exemplars
from block import analyze_blocks, get_mono
from block import block_relation, same_block
from block import is_empty_block, get_block, BOTTOM_BLOCK
from cross_mrca import compute_cross_mrcas, get_cross_mrca

# Assumes that name matches are already stored in AB.
# Assumes that exemplars are known.

def theorize(AB, compute_exemplars=True):
  # ... exemplar.exemplars(A_iter, B_iter, m_iter) ...
  # TBD: option to read them from a file
  # Assumes links found already  - find_links(AB)
  if compute_exemplars:
    exemplar.find_exemplars(AB)
  # else: read them from a file

  analyze_blocks(AB)               # does set_block(...)
  compute_cross_mrcas(AB)

#-----------------------------------------------------------------------------
# compare: The implementation of the RCC-5 theory of AB (A+B).
# Returns a Predicate.
# central = non-peripheral

def compare(AB, u, v):
  assert get_workspace(u)
  assert get_workspace(v)
  assert separated(u, v)
  rel1 = get_central(AB, u)          # u <= u1
  rel3r = get_central(AB, v)         # v <= v1
  if not (rel1 and rel3r):
    log("trees not connected: %s, %s" % (blurb(u), blurb(v)))
    return predicate(NOINFO, v, "trees not connected")
  u1 = rel1.record
  v1 = rel3r.record
  assert separated(u1, v1)
  rel3 = reverse_predicate(v, rel3r)  # v1 >= v
  # Compare u1 and v1, which are central, in opposite checklists
  rel2 = compare_centrally(AB, u1, v1)   # u1 ? v1
  assert rel2
  # u -> u1 -> v1 -> v
  assert get_workspace(rel3.record)
  return compose_final(u, rel1, rel2, rel3)

def get_central(AB, u):
  u_central = u
  while u_central:
    if get_block(u_central):
      break
    u_central = local_sup(AB, u_central).record
  assert u_central, ("No central", blurb(u))
  if u_central is u:
    return predicate(EQ, u_central)
  else:
    sup = local_sup(AB, u)
    if sup and sup.record is u_central:
      return sup                # u < u_central = sup
    else:
      return predicate(LT, u_central, note="get_central")

# Compare taxa that are known to have nonempty exemplar sets.
# Returns non-None as long as there are any exemplars.

def compare_centrally(AB, u, v):
  assert separated(u, v)
  b1 = get_block(u)
  b2 = get_block(v)
  assert b1 != BOTTOM_BLOCK
  assert b2 != BOTTOM_BLOCK
  ship = block_relation(b1, b2)
  if ship == COMPARABLE:
    #! In same block.  Use ranks to figure out relationship.
    return compare_within_block(AB, u, v)
  else:
    # Different blocks (exemplar sets)
    return predicate(ship, v, note="exemplar set comparison")

# --------------------
# u and v are in opposite checklists but same block

def compare_within_block(AB, u, v):
  uu = unique_in_block(AB, u)
  vv = unique_in_block(AB, v)
  if uu and vv:
    #log("# both unique in block: %s, %s" % (blurb(u), blurb(v)))
    return predicate(EQ, v, "unique in both blocks")

  # synonym < accepted regardless ... ?
  if not is_accepted_locally(AB, u) and is_accepted_locally(AB, v):
    answer = predicate(LE, v, "synonym <= accepted")
    if False:
      log("# %s %s %s" %
          (blurb(u), rcc5_symbol(answer.relation), blurb(v)))
  elif is_accepted_locally(AB, u) and not is_accepted_locally(AB, v):
    answer = predicate(GE, v, "accepted >= synonym")
    if False:
      log("# %s %s %s" %
          (blurb(u), rcc5_symbol(answer.relation), blurb(v)))
  else:
    u_rank_n = ranks.ranks_dict.get(get_rank(u, None), 0)
    v_rank_n = ranks.ranks_dict.get(get_rank(v, None), 0)

    # Assume a treelike model...
    if u_rank_n > 0 and v_rank_n > 0:
      if u_rank_n < v_rank_n:
        answer = predicate(LT, v, "ranks <")
      elif u_rank_n == v_rank_n:
        answer = predicate(EQ, v, "ranks =")
      else:
        answer = predicate(GT, v, "ranks >")
    else:
      answer = predicate(COMPARABLE, v, "missing rank information")
    if False and get_canonical(u) != get_canonical(v):
      # This actually happens quite a bit
      # Jackpot
      log("# Using ranks to decide %s %s %s" %
          (blurb(u), rcc5_symbol(answer.relation), blurb(v)))
  return answer

# True iff there is no v in u's checklist congruent to u.

def unique_in_block(AB, u):
  b = get_block(u)
  sup = get_superior(u, None)
  if sup and get_block(sup.record) == b:
    # parent is also in block b
    return False
  if True:
    return not get_mono(u, None)
  else:
    # method to use if we're not computing 'monotypes'
    for c in get_inferiors(get_outject(u)):
      w = AB.in_left(c) if isinA(AB, u) else AB.in_right(c)
      if get_block(w) == b:
        return False
  return True

def compose_final(u, rel1, rel2, rel3):
  assert rel1                   # could be < or <= or =
  assert rel2                   # one of < > = >< !
  assert rel3                   # could be > or >= or =
  rel13 = compose_paths(u, rel1, compose_paths(rel1.record, rel2, rel3))
  return rel13

def compare_locally(AB, u, v):
  rel = simple.compare_per_checklist(get_outject(u), get_outject(v)) # in A or B
  # Returns NOINFO for sibling + any synonym.
  if rel.relation == NOINFO:
    if same_type_ufs(get_type_uf(u), get_type_uf(v)):
      return predicate(EQ, v, "homotypic synonym")   # treat as aliases.
      # (perhaps INTERSECT instead ??)
    else:
      return predicate(NEQ, v, "heterotypic synonym")  # treat as distinct.
  else:
    return predicate(rel.relation,
                    v,
                    note=rel.note,
                    span=rel.span)  # in A or B


# Similar to compose but assumes... assumptions
# Remember all the nodes involved are central (no synonyms)

def compose_paths(z, rel1, rel2):
  rel = compose_predicates(rel1, rel2)
  if False and (monitor(z) or monitor(rel2.record)):
    log("# %s ; %s -> %s" % (rcc5_symbol(rel1.relation),
                             rcc5_symbol(rel2.relation),
                             rcc5_symbol(rel.relation)))
  return rel

# Both relationships hold

def parallel_relation(ship1, ship2):
  return ship1 & ship2

# RCC-5 relationship across the two checklists
# x and y are in AB
# Could be made more efficient by skipping unused calculations

def cross_lt(AB, u, v):
  ship = compare(AB, u, v).relation
  return ship == LT or ship == LE or ship == PERI

def cross_le(AB, u, v):
  ship = compare(AB, u, v).relation
  return ship == LT or ship == LE or ship == PERI or ship == EQ

# -----------------------------------------------------------------------------

# Returns [v1, v2, ...]

def get_intersecting_species(AB, u):
  inters = []
  ids = set()
  for v in exemplar_opposite_records(AB, u):
    s = get_species(v)
    if s and not s.id in ids:
      ids = ids | {s.id}
      inters.append(s)
  return inters

# Nearest ancestor that is an accepted species

def get_species(u):             # u is in workspace
  AB = get_workspace(u)
  s = local_accepted(AB, u)
  while not is_species(s):
    # TBD: use ranks_dict
    rel = local_sup(AB, s)        # predicate
    if rel: s = rel.record
    else: return None
  return s

# A c = A b c, but A b b != A b  ... ?

def same_protonym(u, v):
  e = get_exemplar(u)
  if e:
    f = get_exemplar(v)
    if f:
      return same_specimens(e, f)
  return False

# Taxon in other checklist having same epithet, and also same rank if possible.
# Might be a synonym or subspecies; use get_species if a species is needed.

def get_buddy(AB, u):
  uf = get_exemplar(u)
  if not uf:
    if monitor(u): log("# get_buddy: No exemplar: %s" % blurb(u))
    return False
  r = uf.payload()
  (_, u2, v) = r
  if in_same_tree(AB, u, v):
    v = u2
  v = local_accepted(AB, v)
  if monitor(u): log("# get_buddy: 1 %s -> %s" % (blurb(u), blurb(v)))
  if get_rank(u, None) != get_rank(v, None):
    # e.g. u rank is species, v rank is subspecies
    q = local_sup(AB, v).record
    if get_rank(u, None) == get_rank(q, None):
      # Ranks don't match - maybe it's a promotion/demotion?
      v = q
  if monitor(u): log("# get_buddy: 2 %s -> %s" % (blurb(u), blurb(v)))
  vf = get_exemplar(v)
  if vf and same_specimens(uf, vf):
    return v
  return None

# -----------------------------------------------------------------------------

# The records on z's "side" corresponding to the exemplars
# in the block for z.  (z is in AB)

def exemplar_opposite_records(AB, z):
  b = get_block(z)
  #if len(b) > 0 and len(b) < 50: log("block %s" % list(b))
  return (sid_to_opposite_record(AB, id, z) for id in b)
