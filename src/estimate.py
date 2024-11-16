#!/usr/bin/env python3

import property as prop
import checklist, workspace, simple

from util import log
from checklist import *
from workspace import *
from simple import BOTTOM, compare_per_checklist
from specimen import get_exemplar_record, get_exemplar, get_exemplar_id
from typify import get_link
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
      for x in checklist.postorder_records(AB.A):
        u = AB.in_left(x)
        rel = find_estimate(AB, u)
        set_estimate(u, rel)
        if rel.relationship == EQ: # Metering
          counts[0] += 1
        else:
          counts[1] += 1
    doit(AB)
    doit(swap(AB))
  findem(AB)
  findem(swap(AB))              # swap is in checklist
  log("# Estimates: %s (=), %s (<)" % (int(counts[0]/2), counts[1]))

# Given a model M, let [u] be the least node in the opposite checklist
# that contains u.  Then we're interested in the minimal [u] taken over
# all models.  Does that work?  Is it unique?

def find_estimate(AB, u):
  ship = EQ
  # If u is in A, find smallest v in B with u <= v (sim. B/A)
  u2 = u
  # Ascend from u up to top
  while True:
    # Skip over peripherals, locating non-peripheral u2
    v = get_cross_mrca(u2, None)
    if v != None:
      break
    sup = local_sup(AB, u2)     # superior acccording to u2's checklist
    if not sup:
      # u2 is top
      return relation(ship, u2, "estimate = top")
    ship = sup.relationship if ship == EQ else LT
    u2 = sup.record

  # Fall through with u2, ship

  a = get_block(u2)
  b = get_block(v)
  assert b >= a
  if b != a:                  # a < b, u2 < v
    return relation(LT, v, "estimate")

  # ship is LT or EQ
  # we'll turn EQ into LE if there are multiple options

  # If u2 is at top of chain, EQ is still an option.  Otherwise has to be LE.
  sup = local_sup(AB, u2)       # Synonym?
  if sup and a == get_block(sup.record):
    if ship == EQ: ship = LE

  # Iterate v through superior chain up to top

  while True:
    # Take a look at the block rootward of v
    sup = local_sup(AB, v)         # B.Tupaia montana
    if not sup:
      return relation(ship, v, "estimate = top")
    v2 = sup.record
    b2 = get_block(v2)
    if b2 > b:
      return relation(ship, v, "estimate-chain-top")
    if ship == EQ: ship = LE 
    v = v2
    b = b2

# -----------------------------------------------------------------------------
# u assumed central

# Given u, find unique node v in opposite checklist such that v = u.
# x and y are in AB, not A or B.
# Returns a Relation or None.

def get_equivalent(AB, u):
  assert u
  assert get_workspace(u)
  est = get_estimate(u, None)
  if est and est.relationship == EQ: return est
  else: return None

# -----------------------------------------------------------------------------

# The records on z's "side" corresponding to the exemplars
# in the block for z.  (z is in AB)

def exemplar_records(AB, z):
  return (sid_to_record(AB, id, z) for id in exemplar_ids(AB, z))

def opposite_exemplar_records(AB, z):
  return (sid_to_opposite_record(AB, id, z) for id in exemplar_ids(AB, z))

# record -> list of exemplar ids

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
      exem = get_exemplar_record(u)       # exemplar record (not uf)
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
      mono = True
      for c in get_inferiors(x):  # inferiors in A/B
        b = traverse(c)
        if not is_empty_block(b):
          e = combine_blocks(e, b)
          mono = c if mono == True else None
      if mono != True and mono != None: set_mono(u, AB.in_left(mono))
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
  if b1 != b2:
    assert b1 == b2
  if b1 == BOTTOM_BLOCK:
    assert b1 != BOTTOM_BLOCK

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
