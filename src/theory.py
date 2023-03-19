#!/usr/bin/env python3

import types, functools
import property as prop, checklist, workspace, simple

from util import log
from checklist import *
from workspace import *
from simple import BOTTOM, compare_per_checklist, compare_siblings


import exemplar

# Assumes that name matches are already stored in AB.

def theorize(AB):
  # ... exemplar.exemplars(A_iter, B_iter, m_iter) ...
  # TBD: option to read them from a file
  # Assumes links found already
  exemplar.analyze_exemplars(AB)   # does set_exemplar(...)
  analyze_blocks(AB)               # does set_block(...)
  compute_cross_mrcas(AB)
  find_estimates(AB)

#-----------------------------------------------------------------------------
# compare: The implementation of the RCC-5 theory of AB (A+B).

def compare(AB, v, w):
  if separated(v, w):
    return cross_compare(AB, v, w)
  else:
    return compare_per_checklist(get_outject(v), get_outject(w))

# v and w are Records in AB, not in A or B
# Returns a Relative to w

# v <= v1 ? w1 >= w
def cross_compare(AB, v, w):
  assert separated(v, w)
  rel3r = get_accepted_relation(w)     # w <= w1
  rel1 = get_accepted_relation(v)      # v <= v1
  v1 = rel1.record
  w1 = rel3r.record
  rel3 = reverse_relation(rel3r, w)    # w1 >= w
  # compare rel1 and rel3, which are accepted, in opposite checklists
  assert separated(v1, w1)
  rel2 = compare_accepted(AB, v1, w1)
  if rel2.relationship == EQ:
    # Species and one synonym, or two synonyms
    rel = compare_siblings(v, v1, w1, w)
  else:
    rel = compose_final(rel1, rel2, rel3)
  return rel

# Compare two nodes that are known to be accepted

def compare_accepted(AB, u, v):
  rel1 = get_central(AB, u)          # u <= u1
  rel3r = get_central(AB, v)         # v <= v1
  if not (rel1 and rel3r):
    return relation(NOINFO, v, "trees not connected")
  u1 = rel1.record
  v1 = rel3r.record
  rel3 = reverse_relation(rel3r, v)  # v1 >= v
  assert separated(u1, v1)
  # Compare u1 and v1, which are central, in opposite checklists
  rel2 = compare_centrally(AB, u1, v1)   # u1 ? v1
  assert rel2
  # u -> u1 -> v1 -> v
  return compose_final(rel1, rel2, rel3)

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

def compare_centrally(AB, v, w):
  assert separated(v, w)
  b1 = get_block(v, BOTTOM_BLOCK)
  b2 = get_block(w, BOTTOM_BLOCK)
  assert b1 != BOTTOM_BLOCK
  assert b2 != BOTTOM_BLOCK
  ship = block_relationship(b1, b2)
  if ship != EQ:
    # ship is not EQ (i.e. exemplar sets are different)
    return optimize_relation(AB, v, relation(ship, w, note="different exemplar sets"))
  else:
    #! In same block.  Use names to order.
    return compare_within_block(AB, v, w)

# v and w are inequivalent, but they are in the same nonempty block
# (in parallel chains)

def compare_within_block(AB, v, w):
  assert separated(v, w)
  assert not is_empty_block(get_block(v, BOTTOM_BLOCK))
  assert same_block(get_block(v), get_block(w))
  # Look for a record match for v or w, and punt to simple case

  # Set up a v/m/w/n diagram and see how well it commutes.

  ship1 = ship2 = None

  # Path 1: v = m ? w
  rel_vm = get_estimate(v, None)      # v <= m   in same checklist
  assert rel_vm 
  assert rel_vm.record
  m = rel_vm.record
  assert separated(v, m)
  rel_mw = compare_per_checklist(get_outject(m), get_outject(w)) # in A or B
  rel_vw = compose_paths(rel_vm, rel_mw)
  ship1 = rel_vw.relationship   # v <= m ? w

  # Path 2: v ? n = w   (starting with w)
  rel_wn = get_estimate(w, None)     # n = w
  assert rel_wn
  n = rel_wn.record
  assert separated(w, n)
  rel_nv = compare_per_checklist(get_outject(n), get_outject(v))
  rel_wv = compose_paths(rel_wn, rel_nv)
  rship2 = rel_wv.relationship   # v ? w
  ship2 = reverse_relationship(rship2)

  if ship1==None and ship2==None:
    return None                 # can't get there from here
  if ship1!=None and ship2==None:
    return rel_vw
  if ship2!=None and ship1==None:
    return reverse_relation(w, rel_nv)

  # Take intersection to see where they agree
  ship = parallel_relationship(ship1, ship2)
  if ship1 != ship2:
    if ship != ship1:
      log("# (%s, %s); Reducing %s to %s" %
          (blurb(v), blurb(w), rcc5_symbol(ship1), rcc5_symbol(ship)))
    elif ship != ship2:
      log("# (%s, %s): Reducing %s to %s" %
          (blurb(v), blurb(w), rcc5_symbol(ship2), rcc5_symbol(ship)))

  if ship & (LT|GT) == LT|GT:
    log("# Path inconsistency from %s to %s" % (blurb(v), blurb(w)))
    log("# v <= m %s w:" % (rcc5_symbol(ship1),))
    log("# v         : %s" % blurb(v)),        # v
    log("#   <= m    : %s" % blurb(rel_vm)),   #   = m
    log("#        ? w: %s" % blurb(rel_mw))    #       ? w
    log("# v      ? w: %s" % blurb(rel_vw))    #       ? w

    log("# w <= n %s v:" % (rcc5_symbol(ship2),))
    log("# w         : %s" % blurb(w)),        # w
    log("#   <= n    : %s" % blurb(rel_wn)),   #   = n
    log("#        ? v: %s" % blurb(rel_nv))    #       ? v

    return relation(ship, w, "unknown order in same block")

  else:
    # All four notes are relevant, actually
    return relation(ship, w, "same block")

def compose_final(rel1, rel2, rel3):
  assert rel1                   # could be < or <= or =
  assert rel2                   # one of < > = >< !
  assert rel3                   # could be > or >= or =
  rel23 = compose_paths(rel2, rel3)
  rel13 = compose_paths(rel1, rel23)
  return rel13

def compose_paths(rel1, rel2):
  ship1 = rel1.relationship
  ship2 = rel2.relationship
  if ship1 == EQ: return rel2
  if ship2 == EQ: return rel1
  if ship1 == DISJOINT: return rel1
  if ship2 == DISJOINT: return rel2
  if ship1 == OVERLAP: return rel1
  if ship2 == OVERLAP: return rel2
  note = compose_notes(rel1.note, rel2.note)
  # now < or <= then >= or >
  if rel1.span <= 1 and rel2.span <= 1:
    # Neither is EQ, so we have <= and >= .. sibling synos.
    # INTERSECT vs. NEQ depends on whether homotypic or not.
    return relation(INTERSECT, rel2.record, note=note)
  else:
    assert ship1 == LT
    assert ship2 == GT
    # < >
    return relation(DISJOINT, rel2.record, note=note)

def parallel_relationship(ship1, ship2):
  return ship1 & ship2

# RCC-5 relationship across the two checklists
# x and y are in AB
# Could be made more efficient by skipping unused calculations

def cross_lt(AB, v, w):
  ship = cross_compare(AB, v, w).relationship
  return ship == LT or ship == LE or ship == PERI

def cross_le(AB, v, w):
  ship = cross_compare(AB, v, w).relationship
  return ship == LT or ship == LE or ship == PERI or ship == EQ

# -----------------------------------------------------------------------------
# TBD: Cache this  (and do not cache equivalents)
# def find_estimates(AB): ...

(really_get_estimate, set_estimate) = prop.get_set(prop.declare_property("estimate"))

def get_estimate(v, default=-17):
  if default == -17:
    probe = really_get_estimate(v)
  else:
    probe = really_get_estimate(v, default)
  return probe

def find_estimates(AB):
  count = 1
  def findem(AB):
    for x in checklist.postorder_records(AB.A):
      v = AB.in_left(x)
      rel = find_estimates_from(AB, v) # Cache it
  findem(AB)
  findem(swap(AB))              # swap is in checklist

def find_estimates_from(AB, u):
  probe = get_estimate(u, None)
  if probe != None: return probe
  assert u, 1
  u_overlap = increase_until_overlap(AB, u)
  if not u_overlap:
    log("# %s is peripheral", blurb(u))
    return None        # very peripheral
  rel = find_estimates_when_central(AB, u_overlap)
  assert rel
  v = rel.record
  assert separated(u_overlap, v)
  if u_overlap != u:
    v = rel.record
    # u is peripheral and can be grafted at v.
    # u might be a synonym of v.
    rel = relation(LT, v, note="peripheral")
    rel = optimize_relation(AB, u, rel)
  set_estimate(u, rel)
  log("# %s estimate is %s" % (blurb(u), blurb(rel)))
  return rel

# Find estimates (>= in opposite checklist) in the case where u is not peripheral.
# 'Central' = 'has a nonempty block' i.e. nonperipheral

def find_estimates_when_central(AB, u):
  log("# Estimating %s ..." % blurb(u))
  probe = really_get_estimate(u, None)
  if probe != None:
    log("#  ... cached: %s", (blurb(probe)))
    return probe
  v = get_cross_mrca(u, None)   # u <= v
  assert v

  # increase v until its block is >= u's block
  b = get_block(u)
  while not block_le(b, get_block(v)):
    sup = get_superior(v)
    assert sup, (blurb(v), "at top")
    v = sup.record
  assert v

  # u's block <= v's block

  # overshoot is OK.  done.
  if block_lt(u, get_block(v)):
    rel = relation(LT, v, "in smaller block")
    log("#  ... in smaller block: %s" % blurb(rel))
  # u's block == v's block

  # If no choice, go ahead and match.  Might still be unsound...
  elif at_top_of_chain(AB, u) and at_top_of_chain(AB, v):
    rel = relation(EQ, v, note="unique same-block match")
    log("#  ... chain top: %s" % blurb(rel))

  # "ladder" - lineage from u (bottom of block) up to top of block,
  # and from v (bottom of block) up to top of block.
  else:
    # End of u1 loop - go back to the original u
    rel = really_get_estimate(u, None)
    if not rel:
      analyze_ladder(u, v)
      rel = really_get_estimate(u, None)
      assert rel
    if rel and rel.record: assert separated(u, rel.record)

    if rel:
      log("#  ... optimized: %s", blurb(rel))
      return optimize_relation(AB, u, rel)
    log("#  ... not recovered")
  return rel

# Zip up parallel lineages in the two checklists

def analyze_ladder(u, v):
  u1 = u
  w1 = v
  assert w1
  while u1:   # v still in u's block ... ?
    rel = None
    m = get_match_in_block(u1)
    if m:
      while w1:
        n = get_match_in_block(w1)
        if not n:
          continue

        # Yes, u1 matches something, call it m!
        if m == w1:
          # u1 matches m = w1!  Match is symmetric, so w1 also matches u.
          assert n == u
          rel = relation(EQ, w1, note="match in block")
        else:
          rel = relation(NEQ, w1, note="mismatch in block")
        w1 = get_superior_in_block(w1)
        # End w1 loop

    if not rel:
      # No match for u1, no information on upper bound (!?).   Tombstone for caching.
      rel = False   # relation(NOINFO, None, "no match in block")
    set_estimate(u1, rel)
    u1 = get_superior_in_block(u1)

# u assumed central

def get_superior_in_block(u):
  assert u
  AB = get_source(u)
  sup = local_sup(AB, u)
  if not sup: return None
  u1 = sup.record
  if same_block(get_block(u1), get_block(u)):
    return u1
  else:
    return None

def get_match_in_block(u):
  v = get_link(u, None)     # returns relation in AB
  if v and same_block(get_block(v), get_block(u)):
    return v
  return None

# This performs a bit better but there are still lapses.
# Not consistent with what get_estimate does.

def get_match_in_block_hairy(u):
  b = get_block(u, None)
  if b == None: return None     # peripheral
  m = get_link(u)
  if m:
    n = get_link(m)
    if n and get_block(n, None) == b:
      return m
  return None

# Given u, find unique node v in opposite checklist such that v = u.
# x and y are in AB, not A or B.
# Returns a Relation or None.

def get_equivalent(AB, u):
  assert u
  est = get_estimate(u, None)
  if est and est.relationship == EQ: return est
  else: return None


# Attempt to set span to 1 (parent/child) if possible.

def optimize_relation(AB, v, rel):
  assert separated(v, rel.record)
  w = rel.record
  sup = find_cross_sup_rel(AB, v, w)
  if sup:
    # Make a copy of within-checklist relationship, retaining note
    assert sup.relationship & rel.relationship != 0
    return relation(sup.relationship, w, note=sup.note, span=sup.span)
  else:
    return rel
  
# Cannot use find_estimate because it relies on optimize_relation!
# Warning: rel can be in either checklist... we really just care
# about ship and note; we already know the target will be w

def find_cross_sup_rel(AB, v, w):
  # Are v/w in a child/parent configuration?
  assert v
  assert w
  assert separated(v, w)
  q = None
  # Option 1. v -> p = q ?= w
  rel = local_sup(AB, v)
  if rel:
    p = rel.record
    p_eq = get_equivalent(AB, p)
    if p_eq:
      q = p_eq.record
  if not q:
    # Option 2. v = z -> q ?= w
    v_eq = get_equivalent(AB, v)
    if v_eq:
      z = v_eq.record
      rel = local_sup(AB, z)
      if rel:
        q = rel.record
        # oops, could use rel directly without copying it!
  # Option 3. v ?= z = q -> w   -- cannot go from w to child q.
  return rel if q == w else None

def debug_block(AB, v, w, v0):
      log("# v0=v %s" % (v0 == v))
      log("# %s top %s" % (blurb(v), at_top_of_chain(AB, v)))
      log("# %s top %s" % (blurb(w), at_top_of_chain(AB, w)))

def at_top_of_chain(AB, v):
  sup = local_sup(AB, v)
  if not sup: return True
  b = get_block(v)
  # each peripheral has its own chain
  if is_empty_block(b): return True
  return not same_block(get_block(sup.record), b)

# Find an ancestor of v (in A tree) that overlaps the B tree, i.e.
# has a nonempty block
# TBD: Simplify this.
# Better name: skip_peripherals ?

def increase_until_overlap(AB, u):
  while True:
    if not is_empty_block(get_block(u, BOTTOM_BLOCK)):
      break
    rel = local_sup(AB, u)
    assert rel
    u = rel.record
  return u

# Find 'central' node (at least one exemplar) in SAME checklist,
# and v's relation to it

def get_central(AB, v):
  vc = v
  state = False
  while True:
    if get_cross_mrca(vc, None): # i.e. nonempty block
      break
    sup = local_sup(AB, vc)
    if not sup:
      return None               # No connection
    vc = sup.record
    if not state:
      state = sup               # 'optimize'
    else:
      state = True              # Lost chance to 'optimize'
  if state == True:
    return relation(LT, vc, note="get_central")
  elif state:
    return state
  else:
    return relation(EQ, vc)

# -----------------------------------------------------------------------------
# ONLY USED WITHIN THIS FILE

# estimate of v (in AB) = w (in AB) satisfies: 
#   1. if y = estimate(x) then x ~<= y  ?
#   2. if furthermore x' = estimate(y), then
#      (a) if x' <= x then x ≅ y, 
#      (b) if x' > x then x > y.  [why?]
#
# i.e.   x ≲ y ≲ x' ≤ x   →   x ≅ y   (same 'block')
# i.e.   x ≲ y ≲ x' > x   →   x < y
#
# Cached in AB nodes under the 'estimate' property.
# Needed for equivalent and cosuperior calculations

def compute_cross_mrcas(AB):
  def do_cross_mrcas(AB):
    def traverse(x):            # arg in A, result in B
      u = AB.in_left(x)         # in AB
      exem = exemplar.get_exemplar(u)       # in AB
      if exem:
        (id, u1, v1) = exem
        v = v1 if separated(u, v1) else u1
        m = get_outject(v)
      else:
        m = BOTTOM                # identity for mrca
      for c in get_inferiors(x):  # c in A
        q = traverse(c)       # in B
        m0 = m
        m = simple.mrca(m, q)      # in B
      # Sanity checks
      b1 = get_block(u, BOTTOM_BLOCK)
      if m == BOTTOM:
        assert is_empty_block(b1)
      else:
        assert not is_empty_block(b1)
        v = AB.in_right(m)
        assert separated(u, v)

        b2 = get_block(v, BOTTOM_BLOCK)
        # AssertionError: ([122, 123, 124], [124])
        assert block_le(b1, b2), (list(b1), list(b2))

        #if is_toplike(v):
        #  log("# toplike cross_mrca of %s = %s" % (blurb(u), blurb(v)))
        set_cross_mrca(u, v)
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
# Assumes exemplars have already been chosen and are available
# via `get_exemplar`.

def analyze_blocks(AB):
  def traverse(x, in_left, inB):
    # x is in A (or B if inB)
    if monitor(x): log("# theory: computing block for %s" % (blurb(x),))
    # initial e = exemplars from descendants
    e = BOTTOM_BLOCK
    for c in get_inferiors(x):  # inferiors in A/B
      e = combine_blocks(e, traverse(c, in_left, inB))
      if monitor(c):
        log("# theory: got subblock %s -> %s" % (blurb(c), len(e),))

    u = in_left(x)              # in A (or B if in B)
    exem = exemplar.get_exemplar(u)
    if exem:
      (id, _, _) = exem
      e = adjoin_exemplar(id, e)
    # **************** TBD
    if e != BOTTOM_BLOCK:
      set_block(u, e)
    #log("# theory: block(%s) = %s" % (blurb(x), e))
    return e
  traverse(AB.A.top, AB.in_left, False)
  traverse(AB.B.top, AB.in_right, True)

  # Sanity check
  b1 = get_block(AB.in_left(AB.A.top), BOTTOM_BLOCK)
  b2 = get_block(AB.in_right(AB.B.top), BOTTOM_BLOCK)
  assert b1 != BOTTOM_BLOCK
  assert b1 == b2

# For debugging

def show_exemplars(z, tag, AB):
  def foo(id):
    return blurb(exemplar_record(AB, id, z))
  log("# theory: %s: {%s}" % (tag, ", ".join(map(foo, get_block(z, BOTTOM_BLOCK)))))

# Apply this to something that 'belongs' to a single exemplar
# id -> record in same checklist as z

def exemplar_record(AB, id, z):
  (_, vs, ws) = AB.exemplar_records[id]
  return vs[0] if isinA(AB, z) else ws[0]

# id -> record not in same checklist as z

def opposite_exemplar_record(AB, id, z):
  (_, vs, ws) = AB.exemplar_records[id]
  return vs[0] if isinB(AB, z) else ws[0]

# record -> list of exemplar ids

def exemplar_ids(AB, z):
  return list(get_block(z, BOTTOM_BLOCK))

# The records in z's checklist corresponding to the exemplars
# in the block for z.

def exemplar_records(AB, z):
  return (exemplar_record(AB, id, z) for id in exemplar_ids(AB, z))

def opposite_exemplar_records(AB, z):
  return (opposite_exemplar_record(AB, id, z) for id in exemplar_ids(AB, z))

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
