#!/usr/bin/env python3

import types, functools
import property as prop, checklist, workspace, simple

from util import log
from simple import BOTTOM, compare_in_checklist, compare_siblings
from checklist import *
from workspace import *

import exemplar

# Assumes that name matches are already stored in AB.

def theorize(AB):
  # ... exemplar.exemplars(A_iter, B_iter, m_iter) ...
  # TBD: option to read them from a file
  exemplar_records = {}
  for (xid, u, v) in exemplar.choose_exemplars(AB):
    if u and v:
      set_exemplar(u, xid)
      set_exemplar(v, xid)
      exemplar_records[xid] = (u, v)
    elif u:
      set_exemplar(u, xid)
      if not xid in exemplar_records:
        exemplar_records[xid] = (u, v)
    elif v:
      set_exemplar(v, xid)
      if not xid in exemplar_records:
        exemplar_records[xid] = (u, v)
    
  AB.exemplar_records = exemplar_records
  analyze_blocks(AB)                  # sets of exemplars
  compute_cross_mrcas(AB)
  find_estimates(AB)

#-----------------------------------------------------------------------------
# compare: The implementation of the RCC-5 theory of AB (A+B).

def compare(AB, v, w):
  if separated(v, w):
    return cross_compare(AB, v, w)
  else:
    return compare_in_checklist(get_outject(v), get_outject(w))

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

def compare_accepted(AB, v, w):
  rel1 = get_central(AB, v)          # v <= v1
  rel3r = get_central(AB, w)         # w <= w1
  if not (rel1 and rel3r):
    return relation(NOINFO, w, "trees not connected")
  v1 = rel1.record
  w1 = rel3r.record
  rel3 = reverse_relation(rel3r, w)  # w1 >= w
  assert separated(v1, w1)
  # Compare v1 and w1, which are central, in opposite checklists
  rel2 = compare_centrally(AB, v1, w1)   # v1 ? w1
  # v -> v1 -> w1 -> w
  return compose_final(rel1, rel2, rel3)

def compose_final(rel1, rel2, rel3):
  assert rel1
  assert rel2
  assert rel3
  return compose_relations(rel1, compose_relations(rel2, rel3))

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

def compare_centrally(AB, v, w):
  assert separated(v, w)
  ship = block_relationship(get_block(v, BOTTOM_BLOCK),
                            get_block(w, BOTTOM_BLOCK))
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

  # Path 1: v = m ? w
  rel_vm = get_estimate(v, None)      # v <= m   in same checklist
  assert rel_vm
  m = rel_vm.record
  assert separated(v, m)
  rel_mw = compare_in_checklist(get_outject(m), get_outject(w)) # in A or B
  rel_vw = compose_relations(rel_vm, rel_mw)
  ship1 = rel_vw.relationship   # v <= m ? w

  # Path 2: v ? n = w   (starting with w)
  rel_wn = get_estimate(w, None)     # n = w
  assert rel_wn
  n = rel_wn.record
  assert separated(w, n)
  rel_nv = compare_in_checklist(get_outject(n), get_outject(v))
  rel_wv = compose_relations(rel_wn, rel_nv)
  rship2 = rel_wv.relationship   # v ? w

  ship2 = reverse_relationship(rship2)

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
  if is_toplike(v):
    return None                 # genuine, there is no such w
  else:
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
  findem(swap(AB))

def find_estimates_from(AB, v):
  probe = get_estimate(v, None)
  if probe != None: return probe
  assert v, 1
  vo = increase_until_overlap(AB, v)
  if not vo: return None        # very peripheral
  rel = find_estimates_when_central(AB, vo)
  assert rel
  w = rel.record
  assert separated(vo, w)
  if vo != v:
    w = rel.record
    # v is peripheral and can be grafted at w.
    # v might be a synonym of w.
    rel = relation(LT, w, note="peripheral")
    rel = optimize_relation(AB, v, rel)
  set_estimate(v, rel)
  return rel

# Find estimates (>= in opposite checklist) in the case where v is not peripheral.
# 'Central' = 'has a nonempty block' i.e. nonperipheral

def find_estimates_when_central(AB, v):
  probe = get_estimate(v, None)
  if probe != None: return probe
  w = get_cross_mrca(v, None)   # v <= w

  # increase w until its block is >= v's block
  b = get_block(v)
  while not block_le(v, get_block(w)):
    sup = get_superior(w)
    assert sup, (blurb(w), "at top")
    w = sup.record

  # v's block <= w's block

  # overshoot is OK.  done.
  if block_lt(v, get_block(w)):
    rel = relation(LT, w, "in smaller block")
  # v's block == w's block

  # If no choice, go ahead and match.  Might still be unsound...
  elif at_top_of_chain(AB, v) and at_top_of_chain(AB, w):
    rel = relation(EQ, w, note="unique same-block match")

  # "ladder" - lineage from v (bottom of block) up to top of block,
  # and from w (bottom of block) up to top of block.
  else:
    v1 = v
    w1 = w
    while v1:   # w still in v's block ... ?
      rel = None
      m_rel = get_match_in_block(v1, None)
      if m_rel:
        while w1:
          n_rel = get_match_in_block(w1, None)
          if not n_rel:
            continue

          # Yes, v1 matches something, call it m!
          m = m_rel.record
          if m == w1:
            # v1 matches m = w1!  Match is symmetric, so w1 also matches v.
            assert n_rel and n_rel.record == v
            rel = relation(EQ, w1, note="match in block")
          else:
            rel = relation(NEQ, None, note="mismatch in block")
          w1 = get_superior_in_block(w1)
          # End w1 loop

      if not rel:
        # No match for v1, no information on upper bound (!?).   Tombstone for caching.
        rel = relation(NOINFO, None, "no match in block")
      set_estimate(v1, rel)
      v1 = get_superior_in_block(v1)

    # End of v1 loop - go back to the original v
    rel = get_estimate(v, None)
    assert rel
    if rel.record: assert separated(v, rel.record)

  return optimize_relation(AB, v, rel)

# v assumed central

def get_superior_in_block(v):
  sup = local_sup(v1, None)
  if not sup: return None
  v1 = sup.record
  if same_block(get_block(v1), get_block(v)):
    return v1
  else:
    return None

def get_match_in_block(v):
  rel = get_match(v, None)     # returns relation in AB
  if rel:
    w = rel1.record
    if w and same_block(get_block(w), b):
      return w
  return None

# This performs a bit better but there are still lapses.
# Not consistent with what get_estimate does.

def get_match_in_block_hairy(v):
  b = get_block(v, None)
  if b == None: return None     # peripheral
  ms0 = get_matches(v)
  ms = [m for m in ms0 if get_block(m, None) == b]
  if len(ms) == 1:
    m = ms[0]
    ns0 = get_matches(m)
    ns = [n for n in ns0 if get_block(n, None) == b]
    #log("#   %s matches: %s" % (blurb(m), "; ".join(map(blurb, ns))))
    if len(ns) == 1:
      n = ns[0]
      if n == v:
        if len(ms) < len(ms0):
          log("# Used block to disambiguate %s between %s" %
              (blurb(v), ";".join(map(blurb, ms))))
        if len(ns) < len(ns0):
          log("# Used block to disambiguate %s between %s (return)" %
              (blurb(m), ";".join(map(blurb, n))))
        return m
  return None

# Given v, find unique node w in opposite checklist such that w = v.
# x and y are in AB, not A or B.
# Returns a Relation or None.

def get_equivalent(AB, v):
  assert v
  est = get_estimate(v, None)
  if est and est.relationship == EQ: return est.record
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

def increase_until_overlap(AB, v):
  while True:
    if not is_empty_block(get_block(v, BOTTOM_BLOCK)):
      break
    rel = local_sup(AB, v)
    if not rel: return None
    assert rel                  # Assume some intersection
    v = rel.record
    #assert not(is_toplike(v))
  return v

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
      v = AB.in_left(x)         # in AB
      probe = get_exemplar(v, None)       # in AB
      if probe:
        w_rel = get_matched(v)
        assert w_rel
        w = w_rel.record
        m = get_outject(w)
      else:
        m = BOTTOM                # identity for mrca
      for c in get_inferiors(x):  # c in A
        q = traverse(c)       # in B
        m0 = m
        m = simple.mrca(m, q)      # in B
        #if m and is_top(m) and m != m0:
        #  log("# toplike mrca of %s, %s = %s" % (blurb(m0), blurb(q), blurb(m)))

      # Sanity checks
      b1 = get_block(v, BOTTOM_BLOCK)
      if m == BOTTOM:
        assert is_empty_block(b1)
      else:
        assert not is_empty_block(b1)
        w = AB.in_right(m)
        assert separated(v, w)

        b2 = get_block(w, BOTTOM_BLOCK)
        # AssertionError: ([122, 123, 124], [124])
        assert block_le(b1, b2), (list(b1), list(b2))

        #if is_toplike(w):
        #  log("# toplike cross_mrca of %s = %s" % (blurb(v), blurb(w)))
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
# Assumes exemplars have already been chosen and are available
# via `get_exemplar`.  NOT.

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

    if False and len(e) > 1000:
      log("# analyze_blocks: %s has %s exemplars" % (blurb(x), len(e)))

    v = in_left(x)              # in A (or B if in B)
    xid = get_exemplar(v, None)
    if xid: e = adjoin_exemplar(xid, e)
    # **************** TBD
    if e != BOTTOM_BLOCK:
      set_block(v, e)
    #log("# theory: block(%s) = %s" % (blurb(x), e))
    return e
  traverse(AB.A.top, AB.in_left, False)
  traverse(AB.B.top, AB.in_right, True)

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

# -----------------------------------------------------------------------------
# General workspace utilities (maybe move to workspace.py ?)

def swap(AB):
  BA = AB.swap()
  BA.A = AB.B
  BA.B = AB.A
  BA.exemplar_records = AB.exemplar_records
  return BA

