#!/usr/bin/env python3

import math
import property as prop
import rcc5, checklist, workspace, simple, exemplar, some_exemplar

from util import log
from checklist import *
from workspace import *
from rcc5 import rcc5_symbol

import exemplar
import estimate
from estimate import get_estimate, get_equivalent

# Assumes that name matches are already stored in AB.

def theorize(AB, compute_exemplars=True):
  # ... exemplar.exemplars(A_iter, B_iter, m_iter) ...
  # TBD: option to read them from a file
  # Assumes links found already  - find_links(AB)
  if compute_exemplars:
    exemplar.find_exemplars(AB)
  # else: read them from a file

  estimate.find_estimates(AB)
  analyze_blocks(AB)               # does set_block(...)

#-----------------------------------------------------------------------------
# compare: The implementation of the RCC-5 theory of AB (A+B).
# Returns a relation.

def compare(AB, v, w):
  if separated(v, w):
    return cross_compare(AB, v, w)
  else:
    return simple.compare_per_checklist(get_outject(v), get_outject(w))

# v and w are Records in AB, not in A or B
# Returns a Relative to w

# u <= u1 ? w1 >= w
def cross_compare(AB, u, w):
  assert separated(u, w)
  rel3r = get_accepted_relation(w)     # w <= w1
  rel1 = get_accepted_relation(u)      # u <= u1
  u1 = rel1.record
  w1 = rel3r.record
  rel3 = reverse_relation(w, rel3r)    # w1 >= w
  assert get_workspace(rel3.record) # ***********************
  # compare rel1 and rel3, which are accepted, in opposite checklists
  assert separated(u1, w1)
  rel2 = compare_accepted(AB, u1, w1)
  if rel2.relationship == EQ:
    # Species and one synonym, or two synonyms
    rel = simple.compare_siblings(u, u1, w1, w)
    assert get_workspace(rel.record) # ***********************
  else:
    rel = compose_final(u, rel1, rel2, rel3)
    assert get_workspace(rel.record) # ***********************
  return rel

# Compare two nodes that are known to be accepted

def compare_accepted(AB, u, v):
  assert get_workspace(u)
  assert get_workspace(v)
  rel1 = get_central(AB, u)          # u <= u1
  rel3r = get_central(AB, v)         # v <= v1
  if not (rel1 and rel3r):
    return relation(NOINFO, v, "trees not connected")
  u1 = rel1.record
  v1 = rel3r.record
  rel3 = reverse_relation(v, rel3r)  # v1 >= v
  assert separated(u1, v1)
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
  ship = block_relationship(b1, b2)
  if ship == EQ:
    #! In same block.  Use names to figure out relationships.
    return compare_within_block(AB, u, v)
  else:
    # ship is not EQ (i.e. exemplar sets are different)
    return optimize_relation(AB, u,
                             relation(ship, v, note="exemplar set comparison"))

# u and w are inequivalent, but they are in the same nonempty block
# (in parallel chains)

def compare_within_block(AB, u, v):
  assert separated(u, v)
  assert not is_empty_block(get_block(u))
  assert same_block(get_block(u), get_block(v))
  # Look for a record match for u or v, and punt to simple case

  # Set up a u/m/v/n diagram and see how well it commutes.

  ship1 = ship2 = None

  # Path 1: u = m ? v
  rel_um = get_estimate(u, None)      # u <= m   in same checklist
  m = rel_um.record
  assert separated(u, m)
  rel_mv = compare_locally(m, v)
  rel_uv = compose_paths(u, rel_um, rel_mv)
  ship1 = rel_uv.relationship   # u <= m ? v

  # Path 2: u ? n = v   (starting with v)
  rel_vn = get_estimate(v, None)     # n = v
  n = rel_vn.record
  assert separated(v, n)
  rel_nu = compare_locally(n, u)
  rel_vu = compose_paths(v, rel_vn, rel_nu)
  rev_rel_vu = reverse_relation(v, rel_vu)
  rev_ship2 = rev_rel_vu.relationship   # u ? v

  # Take intersection to see where they agree
  ship = parallel_relationship(ship1, rev_ship2)

  if ((monitor(u) or monitor(v)) and
      ((ship & (LT|GT)) == LT|GT or ship == INCONSISTENT)):
    log("# Path %s from %s to %s" % ("inconsistency" if ship == INCONSISTENT else "incompleteness",
                                     blurb(u), blurb(v)))
    log("#   u        : %s" % blurb(u))         # u
    log("#     %s m    :  %s" % (rcc5_symbol(rel_um.relationship), blurb(rel_um)))
    log("#         %s v:   %s" % (rcc5_symbol(rel_mv.relationship), blurb(rel_mv)))
    log("#   u   %s   v:  %s" % (rcc5_symbol(rel_uv.relationship), blurb(rel_uv)))
    log("#   v        : %s" % blurb(v)),        # v
    log("#     %s n    :  %s" % (rcc5_symbol(rel_vn.relationship), blurb(rel_vn)))
    log("#         %s u:   %s" % (rcc5_symbol(rel_nu.relationship), blurb(rel_nu)))
    log("#   v   %s   u:  %s" % (rcc5_symbol(rel_vu.relationship), blurb(rel_vu)))
    log('')

  if ship == ship1:
    rel = rel_uv
  elif ship == rev_ship2:
    rel = rev_rel_vu
  elif ship == NOINFO:
    log("# Baffled: %s" % rcc5_symbol(ship))
    rel = relation(INTERSECT, v, "tightened")
  elif ship == INCONSISTENT:
    if monitor(u) or monitor(v):
      log("# Recovering %s ! %s from inconsistency" % (blurb(u), blurb(v)))
    rel = relation(INTERSECT, v, "recovered")
  return rel

def compose_final(u, rel1, rel2, rel3):
  assert rel1                   # could be < or <= or =
  assert rel2                   # one of < > = >< !
  assert rel3                   # could be > or >= or =
  rel13 = compose_paths(u, rel1, compose_paths(rel1.record, rel2, rel3))
  return rel13

def compare_locally(u, v):
  rel = simple.compare_per_checklist(get_outject(u), get_outject(v)) # in A or B
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

def get_cross_glb(u):
  v = v1 = get_cross_mrca(u)
  while True:
    if simple.gt(get_cross_mrca(v1), u):
      break
    v = v1
    v1 = get_superior(v1).record
  return v                      # Top of chain

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
# Precompute 'blocks' (exemplar sets, implemented in one of various ways).
# A block is represented as a set of exemplar ids.
# Blocks are stored on nodes in AB.
# Assumes exemplars have already been chosen and are available
# via `get_exemplar`.

def analyze_blocks(ws):
  def doit(AB):
    def traverse(x):
      u = AB.in_left(x)
      if monitor(u): log("# theory: computing block for %s" % (blurb(u),))
      # initial e = exemplars from descendants
      e = BOTTOM_BLOCK
      for c in get_inferiors(x):  # inferiors in A/B
        e = combine_blocks(e, traverse(c))

      exem = some_exemplar.get_exemplar(u) # returns none or (id, u, v)
      if exem:
        e = adjoin_exemplar(exem[0], e)
      # **************** TBD
      set_block(u, e)
      if monitor(u):
        show_exemplars(u, blurb(u), ws)
      return e
    traverse(AB.A.top)
  doit(ws)
  doit(swap(ws))

  # Sanity check
  b1 = get_block(ws.in_left(ws.A.top))
  b2 = get_block(ws.in_right(ws.B.top))
  if b1 != b2:
    assert b1 == b2
  if b1 == BOTTOM_BLOCK:
    assert b1 != BOTTOM_BLOCK

# For debugging

def show_exemplars(z, tag, AB):
  def foo(id):
    return blurb(xid_to_record(AB, id, z))
  log("# theory: %s: {%s}" %
      (tag, ", ".join(map(foo, get_block(z)))))

# -----------------------------------------------------------------------------

def get_intersecting_species(u):
  o = []
  ids = set()
  AB = get_workspace(u)
  for v in opposite_exemplar_records(AB, u):
    s = get_species(v)
    if s and not s.id in ids:
      ids = ids | {s.id}
      o.append(s)
  return o

# The records on z's "side" corresponding to the exemplars
# in the block for z.  (z is in AB)

def exemplar_records(AB, z):
  return (xid_to_record(AB, id, z) for id in exemplar_ids(AB, z))

def opposite_exemplar_records(AB, z):
  return (xid_to_opposite_record(AB, id, z) for id in exemplar_ids(AB, z))

# record -> list of exemplar ids

def exemplar_ids(AB, z):
  return list(get_block(z))

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

# Apply this to an exemplar id to obtain an exemplar union/find node,
# and return the associated taxon record that's in same checklist as z.

def xid_to_record(AB, xid, z):
  uf = AB.exemplar_ufs[xid]
  (_, u, v) = uf.payload()
  return u if isinA(AB, z) else v

def xid_to_opposite_record(AB, xid, z):
  uf = AB.exemplar_ufs[xid]
  (_, u, v) = uf.payload()
  return v if isinA(AB, z) else u

# -----------------------------------------------------------------------------
# Implementation of blocks as Python sets of 'exemplars'.
# A 'block' is just a set of exemplars, implemented as ... a python set.
# The term 'block' comes from the mathematical treatment of partitions.

def get_block(x):
  return really_get_block(x, BOTTOM_BLOCK)

(really_get_block, set_block) = prop.get_set(prop.declare_property("block"))

# RCC-5 relationship between two blocks

def block_relationship(e1, e2):   # can assume intersecting
  if e1 == e2: return EQ          # same block
  elif e1.issubset(e2): return LT
  elif e2.issubset(e1): return GT
  elif e1.isdisjoint(e2): return DISJOINT
  else: return OVERLAP

def same_block(e1, e2):
  return e1 == e2

def block_ge(e1, e2):
  return e1 >= e2

def block_lt(e1, e2):
  return e1 < e1

def block_le(e1, e2):
  return block_ge(e2, e1)

def block_size(e):
  return len(e)

def adjoin_exemplar(exemplar_id, e):
  return combine_blocks(e, {exemplar_id})

# Lattice join (union) of two blocks

def combine_blocks(e1, e2):
  if e1 == BOTTOM_BLOCK: return e2
  if e2 == BOTTOM_BLOCK: return e1
  return e1 | e2

BOTTOM_BLOCK = set()
def same_block(e1, e2): return e1 == e2
def is_empty_block(e): return e == BOTTOM_BLOCK
