#!/usr/bin/env python3

import property as prop
import rcc5, checklist, workspace, simple, exemplar, typify

from util import log
from checklist import *
from workspace import *
from rcc5 import rcc5_symbol

import specimen
import exemplar
import estimate
from specimen import same_specimens, get_exemplar, same_typifications
from estimate import find_estimates, get_estimate, get_equivalent
from estimate import is_empty_block, get_block, BOTTOM_BLOCK
from estimate import block_relationship, same_block, opposite_exemplar_records

# Assumes that name matches are already stored in AB.

def theorize(AB, compute_exemplars=True):
  # ... exemplar.exemplars(A_iter, B_iter, m_iter) ...
  # TBD: option to read them from a file
  # Assumes links found already  - find_links(AB)
  if compute_exemplars:
    exemplar.find_exemplars(get_estimate, AB)
  # else: read them from a file

  find_estimates(AB)

#-----------------------------------------------------------------------------
# compare: The implementation of the RCC-5 theory of AB (A+B).
# Returns a Relation.

def compare(AB, u, v):
  if separated(u, v):           # in different checklists?
    return cross_compare(AB, u, v)
  else:
    return compare_locally(AB, u, v)

# u and v are Records in AB, not in A or B
# Returns a Relative to v (?)

# u <= u1 ? w1 >= w
def cross_compare(AB, u, v):
  assert separated(u, v)
  u1 = local_accepted(AB, u)
  v1 = local_accepted(AB, v)
  assert separated(u1, v1)
  rel = compare_accepted(AB, u1, v1)
  # Cf. simple.compare_per_checklist

  if False:
    if u is u1 and v is v1:
      return rel
    elif u is u1:
      assert not v is v1
      syn = get_superior(v)       # v -> v1
      rev_syn = reverse_relation(v, syn) # v1 -> v
      return compose_relations(rel, rev_syn) # u -> v1 -> v
    elif v is v1:
      syn = get_superior(u)       # u -> u1
      return compose_relations(syn, rel)  # u -> u1 -> v1
    elif u1 is v1:
      # u and v are sibling synonyms...?  No, they're in different checklists
      # NO NO NO NO
      return relation(NOINFO, v, "sibling synonyms")
    else:
      # u and v are not siblings.  probably disjoint?
      pass

  if not u is u1:
    syn = get_superior(u)       # u -> u1
    rel = compose_relations(syn, rel)  # u -> u1 -> v1
  if not v is v1:
    # u1 -> v1 -> v
    syn = get_superior(v)       # v -> v1
    rev_syn = reverse_relation(v, syn) # v1 -> v
    rel = compose_relations(rel, rev_syn) # u1 -> v1 -> v
  return rel

# Compare two nodes that are known to be accepted

def compare_accepted(AB, u, v):
  assert get_workspace(u)
  assert get_workspace(v)
  assert separated(u, v)
  rel1 = get_central(AB, u)          # u <= u1
  rel3r = get_central(AB, v)         # v <= v1
  if not (rel1 and rel3r):
    log("trees not connected: %s, %s" % (blurb(u), blurb(v)))
    return relation(NOINFO, v, "trees not connected")
  u1 = rel1.record
  v1 = rel3r.record
  assert separated(u1, v1)
  rel3 = reverse_relation(v, rel3r)  # v1 >= v
  # Compare u1 and v1, which are central, in opposite checklists
  rel2 = compare_centrally(AB, u1, v1)   # u1 ? v1
  assert rel2
  # u -> u1 -> v1 -> v
  assert get_workspace(rel3.record)
  return compose_final(u, rel1, rel2, rel3)

def maybe_graft(start, rel):
  if is_graft(start, rel):
    return relation(rel.relationship | DISJOINT,
                    rel.record,
                    rel.note,
                    rel.span)
  else:
    return rel

def is_graft(start, rel):
  return (rel.span == 1 and
          is_empty_block(get_block(start)) != is_empty_block(rel.record))

# Compare taxa that are known to have nonempty exemplar sets.
# Returns non-None as long as there are any exemplars.

def compare_centrally(AB, u, v):
  assert separated(u, v)
  b1 = get_block(u)
  b2 = get_block(v)
  assert b1 != BOTTOM_BLOCK
  assert b2 != BOTTOM_BLOCK
  if b1 == b2:
    #! In same block.  Use names to figure out relationships.
    return compare_within_block(AB, u, v)
  else:
    # i.e. exemplar sets are different
    ship = block_relationship(b1, b2)
    return optimize_relation(AB, u,
                             relation(ship, v, note="exemplar set comparison"))

# u and v are inequivalent, but they are in the same nonempty block
# (in parallel chains)

def compare_within_block(AB, u, v):
  assert separated(u, v)
  assert not is_empty_block(get_block(u))
  assert same_block(get_block(u), get_block(v))
  # Look for a record match for u or v, and punt to simple case

  # Set up a u/m/v/n diagram and see how well it commutes.

  ship1 = ship2 = None

  # if get_rank(u) == get_rank(v) ...

  # Path 1: u = m ? v
  rel_um = get_estimate(u, None)      # u <= m   in same checklist
  m = rel_um.record
  assert separated(u, m)
  rel_mv = compare_locally(AB, m, v)
  rel_uv = compose_paths(u, rel_um, rel_mv)
  ship1 = rel_uv.relationship   # u <= m ? v

  # Path 2: u ? n = v   (starting with v)
  rel_vn = get_estimate(v, None)     # n = v
  n = rel_vn.record
  assert separated(v, n)
  rel_nu = compare_locally(AB, n, u)
  rel_vu = compose_paths(v, rel_vn, rel_nu)
  rev_rel_vu = reverse_relation(v, rel_vu)
  rev_ship2 = rev_rel_vu.relationship   # u ? v

  # Take intersection to see where they agree
  ship = parallel_relationship(ship1, rev_ship2)
  rel = relation(ship, v, "parallel")
  if ship == NOINFO:   # probably means sibling synonyms ?
    log("# No info: %s %s,%s %s" % (blurb(u), rcc5_symbol(ship1), rcc5_symbol(rev_ship2), blurb(v)))
  if ship == INCONSISTENT:
    log("# Baffled: %s is inconsistent with %s" %
        (rcc5_symbol(ship1), rcc5_symbol(rev_ship2)))
    log("#   u        : %s" % blurb(u))         # u
    log("#     %s m    :  %s" % (rcc5_symbol(rel_um.relationship), blurb(rel_um)))
    log("#         %s v:   %s" % (rcc5_symbol(rel_mv.relationship), blurb(rel_mv)))
    log("#   u   %s   v:  %s" % (rcc5_symbol(rel_uv.relationship), blurb(rel_uv)))
    log("#   v        : %s" % blurb(v)),        # v
    log("#     %s n    :  %s" % (rcc5_symbol(rel_vn.relationship), blurb(rel_vn)))
    log("#         %s u:   %s" % (rcc5_symbol(rel_nu.relationship), blurb(rel_nu)))
    log("#   v   %s   u:  %s" % (rcc5_symbol(rel_vu.relationship), blurb(rel_vu)))
    log('')

  return rel

def compose_final(u, rel1, rel2, rel3):
  assert rel1                   # could be < or <= or =
  assert rel2                   # one of < > = >< !
  assert rel3                   # could be > or >= or =
  rel13 = compose_paths(u, rel1, compose_paths(rel1.record, rel2, rel3))
  return rel13

def compare_locally(AB, u, v):
  rel = simple.compare_per_checklist(get_outject(u), get_outject(v)) # in A or B
  if rel.relationship & DISJOINT and same_typifications(u, v):
    # rel.relationship is NOINFO or DISJOINT
    # They're not disjoint because type is in both
    return relation(INTERSECT, v, "homotypic synonyms")
  return relation(rel.relationship,
                  v,
                  note=rel.note,
                  span=rel.span)  # in A or B


# Similar to compose but assumes... assumptions
# Remember all the nodes involved are central (no synonyms)

def compose_paths(z, rel1, rel2):
  rel = compose_relations(rel1, rel2)
  if False and (monitor(z) or monitor(rel2.record)):
    log("# %s ; %s -> %s" % (rcc5_symbol(rel1.relationship),
                             rcc5_symbol(rel2.relationship),
                             rcc5_symbol(rel.relationship)))
  return rel

# Both relationships hold

def parallel_relationship(ship1, ship2):
  return ship1 & ship2

# RCC-5 relationship across the two checklists
# x and y are in AB
# Could be made more efficient by skipping unused calculations

def cross_lt(AB, u, v):
  ship = cross_compare(AB, u, v).relationship
  return ship == LT or ship == LE or ship == PERI

def cross_le(AB, u, v):
  ship = cross_compare(AB, u, v).relationship
  return ship == LT or ship == LE or ship == PERI or ship == EQ

# -----------------------------------------------------------------------------

# sort of like: exemplar  ???

# Find a 'central' (non-peripheral) ancestor node (contains at least
# one exemplar) in SAME checklist, and u's relation to it

def get_central(AB, u):
  u_central = u
  while is_empty_block(get_block(u_central)):
    u_central = local_sup(AB, u_central).record
  # "optimize"
  if u_central == u:
    return relation(EQ, u_central)
  else:
    sup = local_sup(AB, u)
    if sup and sup.record == u_central:
      return sup
    else:
      return relation(LT, u_central, note="get_central")

def is_central(u):
  return not is_empty_block(u)

# Attempt to set span to 1 (parent/child) if possible.

def optimize_relation(AB, u, rel):
  assert separated(u, rel.record)
  v = rel.record
  sup = find_cross_sup_rel(AB, u, v)
  if sup:
    # Make a copy of within-checklist relationship, retaining note
    assert sup.relationship & rel.relationship != 0
    return relation(sup.relationship, v, note=sup.note, span=sup.span)
  else:
    return rel
  
# Cannot use find_estimate because it relies on optimize_relation!
# Warning: rel can be in either checklist... we really just care
# about ship and note; we already know the target will be v

def find_cross_sup_rel(AB, u, v):
  # Are u/v in a child/parent configuration?
  assert u
  assert v
  assert separated(u, v)
  q = None
  # Option 1. u -> p = q ?= v
  rel = local_sup(AB, u)
  if rel:
    p = rel.record
    p_eq = get_equivalent(AB, p)
    if p_eq:
      q = p_eq.record
  if not q:
    # Option 2. u = z -> q ?= v
    u_eq = get_equivalent(AB, u)
    if u_eq:
      z = u_eq.record
      rel = local_sup(AB, z)
      if rel:
        q = rel.record
        # oops, could use rel directly without copying it!
  # Option 3. u ?= z = q -> v   -- cannot go from v to child q.
  return rel if q == v else None

# -----------------------------------------------------------------------------

# Returns [v1, v2, ...]

def get_intersecting_species(u):
  inters = []
  ids = set()
  AB = get_workspace(u)
  for v in opposite_exemplar_records(AB, u):
    s = get_species(v)
    if s and not s.id in ids:
      ids = ids | {s.id}
      inters.append(s)
  return inters

def get_species(u):
  AB = get_workspace(u)
  s = local_accepted(AB, u)
  while not is_species(s):
    rel = local_sup(AB, s)        # relation
    if rel: s = rel.record
    else: return None
  return s

def is_species(u):              # z local
  if u == False: return False
  x = get_outject(u)
  return get_rank(x, None) == 'species' and is_accepted(x)

# A c = A b c, but A b b != A b  ... ?

def same_protonym(u, v):
  e = get_exemplar(u)
  if e:
    f = get_exemplar(v)
    if f:
      return same_specimens(e, f)
  return False

# Taxon in v with same epithet, and also same rank if possible.
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

# ----------------------------------------
# Extra stuff

# Least taxon z of AB with u= < z and v <= z

def mrca(AB, u, v):
  while True:
    m = AB.in_left(simple.mrca(get_outject(u), get_outject(get_estimate(v))))
    n = AB.in_right(simple.mrca(get_outject(v), get_outject(get_estimate(u))))
    rel = compare(AB, m, n)
    if   rel.relationship == EQ: return n
    elif rel.relationship == LT: return n
    elif rel.relationship == LE: return n
    elif rel.relationship == GT: return m
    elif rel.relationship == GE: return m
    elif rel.relationship == OVERLAP:
      # Does this terminate?  Yes, because simple.mrca always goes
      # rootward after an overlap.
      log("# Reducing mrca(%s, %s) to mrca(%s, %s)" %
          (blurb(x), blurb(y), blurb(m), blurb(n)))
      u = m
      v = n
    else: assert False

# Least taxon q of B such that w < q 

def get_dominator(AB, w):
  if isinA(AB, w):
    u = w
    v = get_estimate(w)
  else:
    u = get_estimate(w)
    v = w
  p = get_superior(u) if u is w else u
  q = get_superior(v) if v is w else v
  while True:
    rel = compare(AB, p, q)
    if rel.relationship == LT: return p
    elif rel.relationship == LE: return p
    elif rel.relationship == EQ: return q
    elif rel.relationship == GT: return q
    elif rel.relationship == GE: return q
    elif rel.relationship == OVERLAP:
      log("# Overlap: %s with %s" % (blurb(p), blurb(q)))
      # iterate
      p = get_superior(p)       # 'Break' the taxon
    else: assert False
