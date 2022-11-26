#!/usr/bin/env python3

import types, functools
import property as prop, checklist, workspace, simple

import exemplar

from util import log
from simple import BOTTOM
from checklist import *
from workspace import *

def theorize(AB):
  AB.exemplar_records = exemplar.choose_exemplars(AB)
  analyze_blocks(AB)                  # sets of examplars
  compute_cross_mrcas(AB)
  find_reflections(AB)

#-----------------------------------------------------------------------------
# cross_relation: The implementation of the RCC-5 theory of AB (A+B).

# v and w are Records in AB, not in A or B
# Returns a Relative to w

# TBD: set status to "accepted" or "synonym" as appropriate

def cross_relation(AB, v, w):
  #! In same summand
  if isinA(AB, v) == isinA(AB, w):
    rel = simple.simple_relationship(get_outject(v), get_outject(w))
    answer = (rel.relationship, rel.note)
  #! Equivalent (should actually check reflection)
  elif equivalent(v, w):        # special case of get_reflection
    answer = (EQ, "equivalent")
  else:                         # inequivalent somehow
    #! Synonym situations
    psup = local_sup(AB, v)
    qsup = local_sup(AB, w)
    if psup and psup.relationship == SYNONYM:
      if qsup and qsup.relationship == SYNONYM:
        if equivalent(psup.record, qsup.record):
          answer = (OVERLAP, "co-synonyms")
      else:
        if equivalent(psup.record, w):
          answer = (SYNONYM, "is a synonym of")
    else:
      if qsup and qsup.relationship == SYNONYM:
        answer = (SYNONYM, "has synonym")

    #! Deal with peripheral nodes (empty block)
    v_up = increase_until_overlap(AB, v)
    w_up = increase_until_overlap(AB, w)
    if v_up != v and w_up != w:
      answer = (DISJOINT, "peripherals")
    elif v_up != v:
      answer = (PERI, "left peripheral")
    elif w_up != w:
      answer = (IREP, "right peripheral")

    else:
      ship = block_relationship(get_block(v, BOTTOM_BLOCK), get_block(w, BOTTOM_BLOCK))
      if ship == EQ:
        answer = compare_within_block(v, w)
      else:
        # ship is not EQ (i.e. exemplar sets are different)
        answer = (ship, "via exemplar set comparison")

    if answer[0] == EQ:
      # Bug in get_reflection
      log("# Recovered latent equivalence %s = %s, %s" %
          (blurb(v), blurb(w), answer[1]))

  (ship, note) = answer
  if monitor(v) or monitor(w):
    log("# %s %s %s / %s" % (blurb(v), rcc5_symbol(ship), blurb(w), note))

  return relation(ship, w, "articulation", note)

# v and w are inequivalent, but they are in the same nonempty block
# (in parallel chains)

def compare_within_block(v, w):
  assert not is_empty_block(get_block(v, BOTTOM_BLOCK))
  # Look for a record match for v or w, and punt to simple case
  rel1 = rel2 = None
  v_eq = get_equivalent(v, None)  # v, v_eq, v_eq.record in AB (from B)
  if v_eq:
    # v = m ? w
    rel1 = simple.simple_relationship(get_outject(v_eq.record), get_outject(w))
  w_eq = get_equivalent(w, None)   # w, w_eq, w_eq.record in AB
  if w_eq:
    # v ? n = w
    rel2 = simple.simple_relationship(get_outject(v), get_outject(w_eq.record))
  # See whether the v/v_eq/w_eq/w diagram commutes
  if w_eq and v_eq:
    ship = rel1.relationship
    assert rel2.relationship == ship, \
      (blurb(v),
       rcc5_symbol(rel1.relationship), 
       blurb(w_eq.record),
       blurb(v_eq.record),
       rcc5_symbol(rel2.relationship), 
       blurb(w))
    answer = (ship, rel1.note)
  elif v_eq:
    answer = (rel1.relationship, rel1.note)
    # log("# v %s %s %s" % (blurb(v), blurb(rel1), blurb(w)))
  elif w_eq:
    answer = (rel2.relationship, rel2.note)
    # log("# w %s %s %s" % (blurb(v), blurb(rel2), blurb(w)))
  else:
    # Neither v nor w has a suitable equivalent, but they're in same block,
    # so COMPARABLE.   (assume existence of total interleaved chain.)
    answer = (COMPARABLE, "in parallel monotypic chains")
  return answer

# Least node in opposite checklist that's > v
# Returns a relation

def alt_superior(AB, v):
  w = get_reflection(v)     # candidate in AB, w >= v
  if w == BOTTOM:
    return None
  ship = cross_relation(AB, v, w).relationship
  if ship == EQ:
    return get_superior(w, None)
  else:
    assert cross_lt(v, w)
    sup = get_superior(v, None)
    return relation(ship, w, sup.status if sup else "accepted", "alt_superior")

# RCC-5 relationship across the two checklists
# x and y are in AB
# Could be made more efficient by skipping unused calculations

def cross_lt(AB, v, w):
  ship = cross_relation(AB, v, w).relationship
  return ship == LT or ship == LE or ship == PERI

def cross_le(AB, v, w):
  ship = cross_relation(AB, v, w).relationship
  return ship == LT or ship == LE or ship == PERI or ship == EQ

# -----------------------------------------------------------------------------
# Equivalence - a special case of reflection

def equivalent(v, w):
  AB = get_source(v)            # BIG KLUDGE
  # Easier: return v == get_equivalent(w, None)
  # But, let's check symmetry here
  e = really_equivalent(v, w)
  assert e == really_equivalent(w, v), \
    (blurb(v), blurb(w), blurb(get_matched(v)), blurb(get_matched(w)))
  return e

# This gets cached, but shouldn't be

def really_equivalent(v, w):
  ref = get_equivalent(v)
  return ref and ref.record == w

# Given x, returns y such that x = y.  x and y are in AB, not A or B.
# Returns a Relation or None.

def get_equivalent(v, default=None):
  ref = get_reflection(v)
  return ref if ref.relationship == EQ else default

# -----------------------------------------------------------------------------
# Find least w in B such that w >= v (in A).

# TBD: Cache this  (and do not cache equivalents)
# def find_reflections(AB): ...

(get_reflection, set_reflection) = prop.get_set(prop.declare_property("reflection"))

def find_reflections(AB):
  count = 1
  def findem(AB):
    for x in checklist.postorder_records(AB.A):
      v = AB.in_left(x)
      set_reflection(v, find_reflection(AB, v))
  findem(AB)
  findem(swap(AB))

def find_reflection(AB, v):
  if is_toplike(v):
    if isinA(AB, v):
      return relation(EQ, AB.in_right(AB.B.top), "top")
    else:
      return relation(EQ, AB.in_left(AB.A.top), "top")

  # 1. Look for match by name
  b1 = get_block(v, BOTTOM_BLOCK)    # Necessarily nonempty
  nm = get_matched(v)     # returns relation in AB
  if nm and same_block(get_block(nm.record, BOTTOM_BLOCK), b1):
    return nm                 # They're EQ

  # 2. Peripherals don't match
  v1 = increase_until_overlap(AB, v)
  b1 = get_block(v1, BOTTOM_BLOCK)
  assert not is_empty_block(b1)
  w = get_cross_mrca(v1, None)
  if v1 != v:
    return relation(LT, w, "reflection", "peripheral")

  # 3. Detect < when mismatched blocks
  b2 = get_block(w, BOTTOM_BLOCK)
  assert block_le(b1, b2), (len(b1), len(b2), b1, b2)
  if not same_block(b2, b1):
    status = "accepted"   # could also be a synonym? look at v's sup rel
    return relation(LT, w, status, "smaller block")

  # 4. If v and w are both at top of their chains, equate them.
  p_rel = get_superior(v, None)
  q_rel = get_superior(w, None)
  b3 = get_block(p_rel.record, BOTTOM_BLOCK) if p_rel else True
  b4 = get_block(q_rel.record, BOTTOM_BLOCK) if q_rel else True
  if not same_block(b3, b1) and not same_block(b4, b2):
    nm = get_matched(w)
    if nm and same_block(get_block(nm.record, BOTTOM_BLOCK), b2):
      # w matches in v's chain but not vice versa.  Need symmetry
      pass
    else:
      log("# matched chain tops: %s = %s" % (blurb(v), blurb(w)))
      return relation(EQ, w, "equivalent", "at tops of chains")

  # 5. Hmph
  return relation(OVERLAP, w, "reflection", "chains not interleaved")

# Find an ancestor of v (in A tree) that overlaps the B tree, i.e.
# has a nonempty block
# TBD: Simplify this.
# Better name: skip_peripherals ?

def increase_until_overlap(AB, v):
  while True:
    if not is_empty_block(get_block(v, BOTTOM_BLOCK)):
      break
    rel = local_sup(AB, v)
    assert rel                  # Assume some intersection
    v = rel.record
  return v

# -----------------------------------------------------------------------------
# ONLY USED WITHIN THIS FILE

# reflection of v (in AB) = w (in AB) satisfies: 
#   1. if y = reflection(x) then x ~<= y  ?
#   2. if furthermore x' = reflection(y), then
#      (a) if x' <= x then x ≅ y, 
#      (b) if x' > x then x > y.  [why?]
#
# i.e.   x ≲ y ≲ x' ≤ x   →   x ≅ y   (same 'block')
# i.e.   x ≲ y ≲ x' > x   →   x < y
#
# Cached in AB nodes under the 'reflection' property.
# Needed for equivalent and cosuperior calculations

def compute_cross_mrcas(AB):
  def do_cross_mrcas(AB):
    def traverse(x):            # arg in A, result in B
      v = AB.in_left(x)         # in AB
      probe = get_tipward(v, None)       # in AB
      if probe:
        m = get_outject(probe.record)  # in B
      else:
        m = BOTTOM                # identity for mrca
      for c in get_inferiors(x):
        q = traverse(c)       # in B
        m = simple.mrca(m, q)      # in B

      # Sanity checks
      b1 = get_block(v, BOTTOM_BLOCK)
      if m == BOTTOM:
        assert is_empty_block(b1)
      else:
        assert not is_empty_block(b1)
        w = AB.in_right(m)

        b2 = get_block(w, BOTTOM_BLOCK)
        assert not is_empty_block(b2)
        assert block_le(b1, b2), (list(b1), list(b2))

        set_cross_mrca(v, w)
      return m
    traverse(AB.A.top)
  do_cross_mrcas(AB)
  do_cross_mrcas(swap(AB))

(get_cross_mrca, set_cross_mrca) = \
  prop.get_set(prop.declare_property("cross_mrca"))

# -----------------------------------------------------------------------------
# Precompute 'blocks' (exemplar sets, implemented in one of various ways).
# A block is represented as a set of exemplar ids.
# Blocks are stored on nodes in AB.
# Assumes exemplars have already been chosen and are available via `get_exemplar`.

def analyze_blocks(AB):
  def traverse(x, in_left, inB):
    # x is in A (or B if inB)
    if monitor(x): log("computing block for %s" % (blurb(x),))
    # initial e = exemplars from descendants
    e = BOTTOM_BLOCK
    for c in get_inferiors(x):  # inferiors in A/B
      e = combine_blocks(e, traverse(c, in_left, inB))
      if monitor(c): log("got subblock %s -> %s" % (blurb(c), len(e),))
    v = in_left(x)              # in A (or B if in B)
    exemplar = get_exemplar(AB, v, inB)
    if e != BOTTOM_BLOCK:
      set_block(v, e)
    #log("# block(%s) = %s" % (blurb(x), e))
    return e
  traverse(AB.A.top, AB.in_left, False)
  traverse(AB.B.top, AB.in_right, True)

def exemplar_id(ex): return ex[0]

# For debugging

def show_exemplars(z, tag, AB):
  def foo(vw):
    (v, w) = vw
    return blurb(w)
  log("# %s: {%s}" % (tag, ", ".join(map(foo, block_exemplar_records(AB, z)))))

# The records in the B checklist corresponding to the exemplars
# in the block for z.

def exemplar_ids(AB, z):
  return list(get_block(z, BOTTOM_BLOCK))

def exemplar_records(AB, z):
  return (AB.exemplar_records[id] for id in exemplar_ids(AB, z))

(get_block, set_block) = prop.get_set(prop.declare_property("block"))

# -----------------------------------------------------------------------------
# Implementation of blocks as Python sets of 'exemplars'.
# A 'block' is just a set of exemplars, implemented as ... a python set.
# The term 'block' comes from the mathematical treatment of partitions.

# RCC-5 relationship between two blocks

def block_relationship(e1, e2):   # can assume overlap
  if e1 == e2: return EQ          # same block
  elif e1.issubset(e2): return LT
  elif e2.issubset(e1): return GT
  elif e1.isdisjoint(e2): return DISJOINT
  else: return CONFLICT

def same_block(e1, e2):
  return e1 == e2

def block_ge(e1, e2):
  return e1 >= e2

def block_le(e1, e2):
  return block_ge(e2, e1)

def block_size(e):
  return len(e)

def adjoin_exemplar(exemplar_id, e):
  return combine_blocks(e, {exemplar_id})

# Lattice join (union) of two blocks

BOTTOM_BLOCK = set()
def combine_blocks(e1, e2):
  if e1 == BOTTOM_BLOCK: return e2
  if e2 == BOTTOM_BLOCK: return e1
  return e1 | e2
def same_block(e1, e2): return e1 == e2
def is_empty_block(e): return e == BOTTOM_BLOCK

# -----------------------------------------------------------------------------
# General workspace utilities (maybe move to workspace.py ?)

def swap(AB):
  BA = AB.swap()
  BA.A = AB.B
  BA.B = AB.A
  BA.exemplar_records = AB.exemplar_records
  return BA

