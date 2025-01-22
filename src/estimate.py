#!/usr/bin/env python3

import property as prop
import checklist, workspace, simple, ranks

from util import log
from checklist import *
from workspace import *
from simple import BOTTOM, compare_per_checklist
from typify import get_link
from specimen import get_exemplar_info, get_exemplar, get_exemplar_id
from specimen import sid_to_record, sid_to_opposite_record

# -----------------------------------------------------------------------------

(get_estimate, set_estimate) = prop.get_set(prop.declare_property("estimate"))

def find_estimates(AB):
  # Cross-mrcas might need to be replaced 
  #  (LUBs get tighter/smaller on 2nd pass).
  #  Easily overwritten.
  compute_cross_mrcas(AB)          # does set_cross_mrca(...)
  analyze_blocks(AB)               # does set_block(...)

  counts = [0, 0]
  def findem(AB):
    def doit(AB):
      for x in checklist.postorder_records(AB.A): # end with top
        u = AB.in_left(x)
        rel = find_estimate(AB, u)
        set_estimate(u, rel)
        if rel and rel.relationship == EQ: # Metering
          counts[0] += 1
        else:
          counts[1] += 1
    doit(AB)
    doit(swap(AB))
  findem(AB)
  findem(swap(AB))              # swap is in checklist
  log("# Estimates: %s (=), %s (<)" % (int(counts[0]/2), counts[1]))

# Can't make sense of the following:
# "Given a model M, let [u] be the least node in the opposite checklist
# that contains u.  Then we're interested in the minimal [u] taken over
# all models.  Does that work?  Is it unique?"

# Low-numbered ranks are closer to root

# For u in checklist 1, find smallest v in checklist 2 such that u <= v
# Careful, u might be top or None.

def find_estimate(AB, u):
  rel1 = get_central(AB, u)
  # rel1 ascends from u up to top, looking for a node with a cross_mrca
  u_central = rel1.record
  v = get_cross_mrca(u_central)
  # n.b. v >= every exemplar in u
  #  but it could still be < u (or > u)

  # Fall through with u_central, v
  # How does u relate to v?

  a = get_block(u_central)
  u_rank_n = ranks.ranks_dict.get(get_rank(u, 0))

  # Find a v ancestor that's >= u

  while True:
    b = get_block(v)
    if a < b:
      rel2 = relation(LT, v, "estimate by block")
      break

    elif a == b:
      rel2 = compare_in_block(AB, u_central, v)
      ship = rel2.relationship
      if ship == EQ or ship == LT or ship == LE:
        break

    else:   # a > b or a >< b or a ! b
      sup = local_sup(AB, v)         # B.Tupaia montana
      if not sup:                    # v is top
        log("# %s hasn't an estimate" % blurb(u))
        return None
      v = sup.record
      ship = LT
      # iterate

  assert separated(u, rel2.record)
  rel = compose_relations(rel1, rel2) # u -> u_central -> v
  # convert LT to LE when synonym

  return rel

# u and v are in opposite checklists

def compare_in_block(AB, u, v):
  if unique_in_block(u) and unique_in_block(v):
    return relation(EQ, v, "unique in both blocks")

  # This method is stupid and unreliable.  Name comparison and/or a
  # search would be better.

  a = get_block(u)
  u_rank_n = ranks.ranks_dict.get(get_rank(u, 0))

  v_rank_n = ranks.ranks_dict.get(get_rank(v, 0))
  if u_rank_n and v_rank_n:
    if u_rank_n < v_rank_n:
      answer = relation(LT, v, "estimate by rank<")
    elif u_rank_n == v_rank_n:
      answer = relation(EQ, v, "estimate by rank")
    else:
      answer = relation(GT, v, "estimate by rank>")
  else:
    answer = relation(NOINFO, v, "missing rank information")
  return answer

# More complicated: compare names

def unique_in_block(u):
  b = get_block(u)
  sup = get_superior(u, None)
  if sup and get_block(sup.record) == b:
    # parent is also in block b
    return False
  for c in get_inferiors(u):
    if get_block(c) == b:
      return False
  return True

# -----------------------------------------------------------------------------
# u assumed central

# Given u, find unique node v in opposite checklist such that v = u.
# x and y are in AB, not A or B.
# Returns a Relation or None.

def get_equivalent(AB, u):
  assert u
  assert get_workspace(u)
  est = get_estimate(u, None)   # In opposite checklist
  if (est and est.relationship == EQ and
      (is_accepted_locally(AB, u) == is_accepted_locally(AB, est.record))):
    return est
  else: return None

# -----------------------------------------------------------------------------

# The records on z's "side" corresponding to the exemplars
# in the block for z.  (z is in AB)

def exemplar_records(AB, z):    # not used?
  return (sid_to_record(AB, id, z) for id in exemplar_ids(AB, z))

def exemplar_opposite_records(AB, z): # see theory.py
  return (sid_to_opposite_record(AB, id, z) for id in exemplar_ids(AB, z))

# record -> list of ids for subtended exemplars

def exemplar_ids(AB, z):
  return list(get_block(z))

# For debugging

def show_exemplars(z, tag, AB):
  def foo(id):
    return blurb(sid_to_record(AB, id, z))
  log("# estimate: %s: {%s}" %
      (tag, ", ".join(map(foo, get_block(z)))))

# -----------------------------------------------------------------------------

# If w (in B) is the 'estimate' of u (in A), and u contains at least one 
# exemplar, then:
#      w is the smallest taxon in B containing all the exemplars 
#      that are in v (might contain more).
#
# i.e.   x ≲ y ≲ x' ≤ x   →   x ≅ y   (same 'block')
# i.e.   x ≲ y ≲ x' > x   →   x < y
#
# Cached in AB nodes under the 'estimate' property.
# Needed for equivalent and cosuperior calculations

def compute_cross_mrcas(AB):
  def do_cross_mrcas(WS):        # WS is AB or swap(AB)
    def traverse(x):            # arg in A, result in B
      u = WS.in_left(x)          # in WS
      exem = get_exemplar_info(u)       # exemplar record (not uf)
      if exem:
        (_, u1, v1) = exem
        assert get_outject(v1), blurb(v1) # fails
        assert get_outject(u1), blurb(u1)
        m = get_outject(v1) if isinA(AB, u) else get_outject(u1)
      else:
        m = BOTTOM                # identity for mrca
      for c in get_inferiors(x):  # c in WS.A
        q = traverse(c)           # in WS.B
        m = simple.mrca(m, q)     # in WS.B
      # Sanity checks
      if m != BOTTOM:
        v = WS.in_right(m)
        assert separated(u, v)
        set_cross_mrca(u, v)
      return m
    traverse(WS.A.top)
  do_cross_mrcas(AB)
  do_cross_mrcas(swap(AB))

(get_cross_mrca, set_cross_mrca) = \
  prop.get_set(prop.declare_property("cross_mrca"))

# -----------------------------------------------------------------------------

# Find a 'central' (non-peripheral) ancestor node (contains at least
# one exemplar) in SAME checklist.  Returns a Relation.

def get_central(AB, u):
  u_central = u
  while is_empty_block(get_block(u_central)):
    u_central = local_sup(AB, u_central).record
  # Make a Relation.
  if u_central is u:
    return relation(EQ, u_central)
  else:
    sup = local_sup(AB, u)
    if sup and sup.record == u_central:
      return sup                # u -> u_central = sup
    else:
      return relation(LT, u_central, note="get_central")

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
      if monitor(u): log("# estimate: computing block for %s" % (blurb(u),))
      # initial e = exemplars from descendants
      e = BOTTOM_BLOCK
      mono = None
      for c in get_inferiors(x):  # inferiors in A/B
        b = traverse(c)
        if not is_empty_block(b):
          e = combine_blocks(e, b)
          mono = c if mono == None else False
      if mono != None and mono != False: set_mono(u, AB.in_left(mono))
      uf = get_exemplar(u) # returns None or... (sid, u, v)?
      if uf:
        e = adjoin_exemplar(get_exemplar_id(uf), e)
        if get_redundant(x, None):
          log("# Redundant record's exemplar suppressed: %s" % blurb(x))
          return BOTTOM_BLOCK
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
  assert b1 == b2
  assert b1 != BOTTOM_BLOCK
  log("# top block size: %s" % len(b1))

def adjoin_exemplar(exemplar_id, e):
  return combine_blocks(e, {exemplar_id})

# -----------------------------------------------------------------------------
# Implementation of blocks as Python sets of 'exemplars'.
# A 'block' is just a set of exemplars, implemented as ... a python set.
# The term 'block' comes from the mathematical treatment of partitions.

def get_block(x):
  return really_get_block(x, BOTTOM_BLOCK)

(really_get_block, set_block) = prop.get_set(prop.declare_property("block"))
(get_mono, set_mono) = prop.get_set(prop.declare_property("mono"))

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

# Lattice join (union) of two blocks

def combine_blocks(e1, e2):
  if e1 == BOTTOM_BLOCK: return e2
  if e2 == BOTTOM_BLOCK: return e1
  return e1 | e2

BOTTOM_BLOCK = set()
def same_block(e1, e2): return e1 == e2
def is_empty_block(e): return e == BOTTOM_BLOCK
