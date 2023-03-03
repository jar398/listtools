#!/usr/bin/env python3

import types, functools
import property as prop, checklist, workspace, simple

import exemplar

from util import log
from simple import BOTTOM
from checklist import *
from workspace import *

def theorize(AB):
  # ... exemplar.exemplars(A_iter, B_iter, m_iter) ...
  # TBD: option to read them from a file
  AB.exemplar_records = exemplar.choose_exemplars(AB) #list of (id, vs, ws)
  for ex in AB.exemplar_records:
    (_, vs, ws) = ex
    for v in vs: set_exemplar_record(v, ex)
    for w in ws: set_exemplar_record(w, ex)
  analyze_blocks(AB)                  # sets of examplars
  compute_cross_mrcas(AB)
  find_reflections(AB)

(get_exemplar_record, set_exemplar_record) = \
   prop.get_set(prop.declare_property("exemplar_record"))

#-----------------------------------------------------------------------------
# cross_relation: The implementation of the RCC-5 theory of AB (A+B).

# v and w are Records in AB, not in A or B
# Returns a Relative to w

# TBD: set status to "accepted" or "synonym" as appropriate

def cross_relation(AB, v, w):
  #! In same summand?
  if isinA(AB, v) == isinA(AB, w):
    rel = simple.simple_relationship(get_outject(v), get_outject(w))
    answer = (rel.relationship, rel.note)
  else:

    #! Co-synonyms?
    rel = equivalent(local_accepted(AB, v), local_accepted(AB, w))
    if rel:
      answer = simple.sibling_relationship(get_outject(v), get_outject(w))
    else:

      #! Is either one peripheral?
      # Compose v <= v_up ship w_up >= w.  Cases:
      #  v < v_up  ? w_up > w    !
      #  v < v_up <= w_up = w    <
      #  v = v_up >= w_up > w    >
      #  v < v_up  ? w_up = w    !
      #  v = v_up  ? w_up > w    !
      #  v = v_up  ? w_up = w    ?
      v_up = increase_until_overlap(AB, v)
      w_up = increase_until_overlap(AB, w)
      if v_up != v and w_up != w:
        answer = (DISJOINT, "different peripherals")
      elif v_up != v or w_up != w:
        (ship, _) = blocked_relationship(AB, v_up, w_up)
        if v_up != v and (ship & ~LE == 0):
          answer = (PERI, "left peripheral")
        elif w_up != w and (ship & ~GE == 0):
          answer = (IREP, "right peripheral")
        else:
          answer = (DISJOINT, "peripheral elsewhere")
      else:

        #! Do they have distinct exemplar sets?
        answer = blocked_relationship(AB, v, w)

  (ship, note) = answer
  if monitor(v) or monitor(w):
    log("# theory: [%s %s %s] because %s" % (blurb(v), rcc5_symbol(ship), blurb(w), note))

  return relation(ship, w, note=note)

# Compare taxa that both have nonempty exemplar sets.

def blocked_relationship(AB, v, w):
  ship = block_relationship(get_block(v, BOTTOM_BLOCK),
                            get_block(w, BOTTOM_BLOCK))
  if ship != EQ:
    # ship is not EQ (i.e. exemplar sets are different)
    answer = (ship, "different exemplar sets")
  else:
    #! In same block.  Use names to order.
    answer = compare_within_block(AB, v, w)
  return answer

# v and w are inequivalent, but they are in the same nonempty block
# (in parallel chains)

def compare_within_block(AB, v, w):
  assert not is_empty_block(get_block(v, BOTTOM_BLOCK))
  # Look for a record match for v or w, and punt to simple case

  # Set up a v/m/w/n diagram and see how well it commutes.

  # Path 1: v = m ? w
  rel_vm = get_equivalent(v)      # v = m
  if rel_vm:
    m = rel_vm.record
    rel_mw = simple.simple_relationship(get_outject(m), get_outject(w)) # in A or B
    ship1 = rel_mw.relationship   # v ? w

  # Path 2: v ? n = w   (starting with w)
  rel_wn = get_equivalent(w)     # n = w
  if rel_wn:
    n = rel_wn.record
    rel_vn = simple.simple_relationship(get_outject(v), get_outject(n))
    ship2 = rel_vn.relationship   # v ? w

  if rel_vm and rel_wn:
    # Take intersection to see where they agree
    ship = ship1 & ship2
    # Take intersection of relationships (both are true)??
    assert ship != 0, \
      (blurb(v),              # v
       blurb(rel_vm),         #   = m
       blurb(rel_mw),         #       ? w
       'and',                 # v
       blurb(rel_vn),         #   ? n
       cross_relation(n, w))  #       = w
    if ship1 != ship2:
      if ship != ship1:
        log("# (%s, %s); Reducing %s to %s" %
            (blurb(v), blurb(w), rcc5_symbol(ship1), rcc5_symbol(ship)))
      elif ship != ship2:
        log("# (%s, %s): Reducing %s to %s" %
            (blurb(v), blurb(w), rcc5_symbol(ship2), rcc5_symbol(ship)))

    # All four notes are relevant, actually
    answer = (ship, rel_mw.note)

  elif rel_vm:
    # have m, but no equivalent to w
    answer = (ship1, rel_vm.note)
    # log("# theory: v %s %s %s" % (blurb(v), blurb(rel_m), blurb(w)))
  elif rel_wn:
    answer = (ship2, rel_wn.note)
    # log("# theory: w %s %s %s" % (blurb(v), blurb(rel2), blurb(w)))
  else:
    # Neither v nor w has a suitable equivalent, but they're in same block,
    # so COMPARABLE.   (assume existence of total interleaved chain.)
    # Happens a lot e.g. Spermophilopsis leptodactyla
    #log("# theory: Comparable: [%s %s %s]" % (blurb(v), rcc5_symbol(COMPARABLE), blurb(w)))
    answer = (COMPARABLE, "'monotypic chain'")
  return answer

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
# Equivalence - a special case of reflection - returns relation

def equivalent(v, w):
  rel = get_equivalent(v)
  if rel and rel.record != w: rel = None
  if True:
    # Let's check symmetry
    rel2 = get_equivalent(w)
    if rel2 and rel2.record != v: rel2 = None
    assert (not rel) == (not rel2), \
      (blurb(v), blurb(rel),
       blurb(w), blurb(rel2))
  return rel

# Given x, returns y such that x = y.  x and y are in AB, not A or B.
# Returns a Relation or None.

def get_equivalent(v, default=None):
  rel = get_reflection(v)
  return rel if rel.relationship == EQ else default

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

# We conclude that v = w iff we would conclude that w = v

def find_reflection(AB, v):
  if is_toplike(v):
    if isinA(AB, v):
      return relation(EQ, AB.in_right(AB.B.top), note="top")
    else:
      return relation(EQ, AB.in_left(AB.A.top), note="top")

  # 1. Peripherals don't match, period
  vo = increase_until_overlap(AB, v)
  w = get_cross_mrca(vo, None)
  if vo != v:
    # v doesn't overlap the checklist v isn't in
    assert not is_empty_block(get_block(vo, BOTTOM_BLOCK))
    return relation(LT, w, note="peripheral")

  # 2. Detect < when mismatched blocks
  b = get_block(v)    # Necessarily nonempty
  b2 = get_block(w)
  assert block_le(b, b2), (len(b), len(b2), b, b2)
  if not same_block(b2, b):
    return relation(LT, w, note="smaller block")

  # 3. Blocks match.  Look for match within block by name.
  nm = get_matched(v)     # returns relation in AB
  if nm and same_block(get_block(nm.record, BOTTOM_BLOCK), b):
    return nm                 # They're EQ

  # 4. If v and w are both at top of their chains, and neither has
  # a match inside the block, equate them.
  p_rel = local_sup(AB, v)
  b3 = get_block(p_rel.record) if p_rel else True
  if not same_block(b3, b):     # v is at top of its chain
    while True:
      q_rel = local_sup(AB, w)
      b4 = get_block(q_rel.record) if q_rel else True
      if not same_block(b4, b): # w is at top of its chain
        break
      w = q_rel.record
    # v and w are proper children (not monotypes) of their parents
    nm = get_matched(w)
    if nm and same_block(get_block(nm.record, BOTTOM_BLOCK), b):
      # w matches in v's chain but not vice versa.  Asymmetric
      # log("# theory: unmatched goes to matched: %s ? %s" % (blurb(v), blurb(w)))
      return relation(OVERLAP, w, note="w at top of chain has a match")
    else:
      # Both unmatched and both at top of their chains
      assert get_block(w) == b
      return relation(EQ, w, note="tops of chains, names unmatched")

  # 5. Hmph
  return relation(OVERLAP, w, note="v not at top of chain")

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
      probe = get_exemplar_record(v, None)       # in AB
      if probe:
        w_rel = get_matched(v)
        assert w_rel
        w = w_rel.record
        m = get_outject(w)
      else:
        m = BOTTOM                # identity for mrca
      for c in get_inferiors(x):  # c in A
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
        # AssertionError: ([122, 123, 124], [124])
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
# Assumes exemplars have already been chosen and are available via `get_exemplar_record`.

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
      print("# analyze_blocks: %s has %s exemplars" % (blurb(x), len(e)),
            file=sys.stderr)

    v = in_left(x)              # in A (or B if in B)
    exem = get_exemplar_record(v, None) # (id, v, w)
    if exem: e = adjoin_exemplar(exem[0], e)
    if e != BOTTOM_BLOCK:
      set_block(v, e)
    #log("# theory: block(%s) = %s" % (blurb(x), e))
    return e
  traverse(AB.A.top, AB.in_left, False)
  traverse(AB.B.top, AB.in_right, True)

def exemplar_id(ex): return ex[0]

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

