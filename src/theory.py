#!/usr/bin/env python3

import math
import property as prop
import checklist, workspace, simple, exemplar

from util import log
from checklist import *
from workspace import *
from simple import BOTTOM, compare_per_checklist, compare_siblings


import exemplar

# Assumes that name matches are already stored in AB.

def theorize(AB):
  # ... exemplar.exemplars(A_iter, B_iter, m_iter) ...
  # TBD: option to read them from a file
  # Assumes links found already  - find_links(AB)
  exemplar.analyze_exemplars(AB)   # does set_exemplar(...)
  analyze_blocks(AB)               # does set_block(...)
  compute_cross_mrcas(AB)          # does set_cross_mrca(...)
  find_estimates(AB)               # does set_estimate(...)

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
  rel3 = reverse_relation(w, rel3r)    # w1 >= w
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
  rel3 = reverse_relation(v, rel3r)  # v1 >= v
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

def compare_centrally(AB, u, w):
  assert separated(u, w)
  b1 = get_block(u)
  b2 = get_block(w)
  assert b1 != BOTTOM_BLOCK
  assert b2 != BOTTOM_BLOCK
  ship = block_relationship(b1, b2)
  if ship != EQ:
    # ship is not EQ (i.e. exemplar sets are different)
    return optimize_relation(AB, u, relation(ship, w, note="different exemplar sets"))
  else:
    #! In same block.  Use names to order.
    return compare_within_block(AB, u, w)

# u and w are inequivalent, but they are in the same nonempty block
# (in parallel chains)

def compare_within_block(AB, u, w):
  assert separated(u, w)
  assert not is_empty_block(get_block(u))
  assert same_block(get_block(u), get_block(w))
  # Look for a record match for u or w, and punt to simple case

  # Set up a u/m/w/n diagram and see how well it commutes.

  ship1 = ship2 = None

  # Path 1: u = m ? w
  rel_um = get_estimate(u, None)      # u <= m   in same checklist
  assert rel_um, blurb(u)
  m = rel_um.record
  assert separated(u, m)
  rel_mw = compare_per_checklist(get_outject(m), get_outject(w)) # in A or B
  rel_uw = compose_paths(rel_um, rel_mw)
  ship1 = rel_uw.relationship   # u <= m ? w

  # Path 2: u ? n = w   (starting with w)
  rel_wn = get_estimate(w, None)     # n = w
  assert rel_wn
  n = rel_wn.record
  assert separated(w, n)
  rel_nu = compare_per_checklist(get_outject(n), get_outject(u))
  rel_wu = compose_paths(rel_wn, rel_nu)
  rship2 = rel_wu.relationship   # u ? w
  ship2 = reverse_relationship(rship2)

  if ship1==None and ship2==None:
    return None                 # can't get there from here
  if ship1!=None and ship2==None:
    return rel_uw
  if ship2!=None and ship1==None:
    return reverse_relation(w, rel_nu)

  # Take intersection to see where they agree
  ship = parallel_relationship(ship1, ship2)
  if ship1 != ship2:
    if ship != ship1:
      if monitor(u):
        log("# (%s, %s); Reducing %s to %s" %
            (blurb(u), blurb(w), rcc5_symbol(ship1), rcc5_symbol(ship)))
    elif ship != ship2:
      if monitor(u):
        log("# (%s, %s): Reducing %s to %s." %
            (blurb(u), blurb(w), rcc5_symbol(ship2), rcc5_symbol(ship)))

  if ship & (LT|GT) == LT|GT:
    log("# Path inconsistency from %s to %s" % (blurb(u), blurb(w)))
    log("# u <= m %s w:" % (rcc5_symbol(ship1),))
    log("# u         : %s" % blurb(u)),        # u
    log("#   <= m    : %s" % blurb(rel_um)),   #   = m
    log("#        ? w: %s" % blurb(rel_mw))    #       ? w
    log("# u      ? w: %s" % blurb(rel_uw))    #       ? w

    log("# w <= n %s u:" % (rcc5_symbol(ship2),))
    log("# w         : %s" % blurb(w)),        # w
    log("#   <= n    : %s" % blurb(rel_wn)),   #   = n
    log("#        ? u: %s" % blurb(rel_nu))    #       ? u

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
    ship = ship1 | ship2
    if ship == LT|GT:
      return relation(DISJOINT, rel2.record, note=note)
    else:
      return relation(ship, rel2.record, note=note)

def parallel_relationship(ship1, ship2):
  return ship1 & ship2

# RCC-5 relationship across the two checklists
# x and y are in AB
# Could be made more efficient by skipping unused calculations

def cross_lt(AB, u, w):
  ship = cross_compare(AB, u, w).relationship
  return ship == LT or ship == LE or ship == PERI

def cross_le(AB, u, w):
  ship = cross_compare(AB, u, w).relationship
  return ship == LT or ship == LE or ship == PERI or ship == EQ

# -----------------------------------------------------------------------------

# 'nearness' - reciprocal of distance, approximately - it's thresholded,
# so only matters whether it's small or large

# For mammals, tip to root is expected to be about 13... -> nearness 3
# For 2M species, tip to root is expected to be about 20... -> nearness 2

def compute_nearness(u, v):
  return -compute_distance(u, v)

def compute_distance(u, v):
  return (compute_half_distance(u, v) +
          compute_half_distance(v, u))/2.

def compute_half_distance(u, v):
  # u < u1 <= (v1 < m > v)
  u1 = get_center(u)
  v1 = get_cross_mrca(u1)
  m = simple.mrca(v1, v)
  dist = (distance_on_lineage(u, u1) +
          distance_on_lineage(v1, m) +
          distance_on_lineage(v, m))
  return dist

def distance_on_lineage(u1, u2):
  if u1 == u2:
    return 0
  return (distance_to_parent(u) +
          distance_on_lineage(get_superior(u).record, u2))

def distance_to_parent(u):
  sup = get_superior(u)
  lg(len(get_children(sup.record)))

def lg(x):
  return math.log(x)/log2
log2 = math.log(2)

# sort of like: exemplar

def get_cross_glb(u):
  v = v1 = get_cross_mrca(u)
  while True:
    if simple.gt(get_cross_mrca(v1), u):
      break
    v = v1
    v1 = get_superior(v1).record
  return v                      # Top of chain

# -----------------------------------------------------------------------------

(get_estimate, set_estimate) = prop.get_set(prop.declare_property("estimate"))

def find_estimates(AB):
  count = 1
  def findem(AB):
    def doit(AB):
      if True:
        for x in checklist.postorder_records(AB.A):
          u = AB.in_left(x)
          if False:
            find_estimates_from(AB, u) # Cache it
          else:
            v = find_estimate(AB, u)
            set_estimate(u, v)
    doit(AB)
    doit(swap(AB))
  findem(AB)
  findem(swap(AB))              # swap is in checklist

def find_estimate(AB, u):
  notes = []
  u1 = u
  peri = under = scan = False
  while True:
    v = get_cross_mrca(u1, None)
    if v != None:
      break
    sup = local_sup(AB, u1)
    assert sup, (blurb(u), blurb(u1))
    u1 = sup.record
    peri = True
  # v = xmrca(u), non-None
  while descends_from(get_cross_mrca(v), u):
    sup = local_sup(AB, v)
    if not sup:
      break
    v = sup.record
    under = True
  # v is >= xmrca(u)
  while True:
    sup = local_sup(AB, v)
    if not sup:
      break
    if not get_cross_mrca(sup.record) is u:
      break
    v = sup.record
    scan = True
  note = '+'.join((name
                   for (flag, name)
                   in ((peri, "peri"), (under, "under"), (scan, "scan"))
                   if flag))
  # v is >= u
  return relation(LE, v, '+'.join(notes))

def descends_from(u1, u2):
  return simple.descends_from(get_outject(u1), get_outject(u2))

# -----------------------------------------------------------------------------

# Set estimates per find_estimates_when_central, adjusting for
# peripherals as needed.

def find_estimates_from(AB, u):
  assert u, blurb(u)
  probe = get_estimate(u, None)    # see whether cached already
  if probe != None: return probe
  # Nearest opposite non-peripheral
  u_central = get_central(AB, u).record
  #log("# Central of %s is %s" % (blurb(u), blurb(u_central)))
  assert is_central(u_central)
  assert get_cross_mrca(u_central, None), \
    (blurb(u_central), get_block(u_central))
  rel = find_estimates_when_central(AB, u_central)
  assert rel
  v = rel.record
  assert separated(u_central, v)
  if u_central != u:
    # u is peripheral and can be grafted at u_central -> v.
    # E.g. u might be a synonym of v.
    rel = optimize_relation(AB, u, relation(LT, v, note="peripheral"))
    #log("# Peripheral %s estimate %s." % (blurb(u), blurb(rel)))
    set_estimate(u, rel)

# Find estimates (>= in opposite checklist) in the case where u is not peripheral.
# 'Central' = 'has a nonempty block' i.e. nonperipheral.
# The estimate for u *always* gets set to something non-None.
# Possibly some of u's ancestors too.

def find_estimates_when_central(AB, u):
  if monitor(u):
    log("# Estimating central %s ..." % blurb(u))
  probe = get_estimate(u, None)
  if probe != None:
    if monitor(u):
      log("#  ... cached: %s" % (blurb(probe)))
    return probe

  v = get_cross_mrca(u, None)   # u <= v
  assert v, (blurb(u), get_block(u))

  # v's block should be >= u's block by construction
  assert block_le(get_block(u), get_block(v))

  # overshoot is OK.  done.
  if not same_block(get_block(u), get_block(v)):
    rel = relation(LT, v, "block is smaller")
    rel = optimize_relation(AB, u, rel)
    if monitor(u):
      log("#  ...estimate: %s %s" % (blurb(u), blurb(rel)))
    set_estimate(u, rel)

  # If no choice, go ahead and match.  Might still be unsound...
  elif at_top_of_chain(AB, u) and at_top_of_chain(AB, v):
    rel = relation(EQ, v, note="unique same-block match")
    if monitor(u):
      log("#  ... chain top: %s" % blurb(rel))
    set_estimate(u, rel)

  # "ladder" - lineage from u (bottom of block) up to top of block,
  # and from v (bottom of block) up to top of block.
  else:
    if monitor(u): log("#  ... ladder ...")

    if get_estimate(u, None) == None:
      # End of u1 loop - go back to the original u
      analyze_ladder(u, v)        # sets estimates
      # Check to make sure it got cached
      rel = get_estimate(u, None)
      assert rel != None, (blurb(u), "not cached")
      if monitor(u):
        log("#  ... ladder: %s %s" % (blurb(u), blurb(rel)))
    else:
      if monitor(u):
        log("#  ... ladder cached")

  assert rel, blurb(u)
  assert separated(u, rel.record)
  return rel


# Zip up parallel lineages in the two checklists.

def analyze_ladder(u, v):

  u_up = get_superior_in_block(u)
  v_up = get_superior_in_block(v)

  if not u_up and not v_up:
    # Top of chain, last resort
    set_estimate(u, relation(EQ, v, note="tops of chains"))
    set_estimate(v, relation(EQ, u, note="tops of chains"))
    if monitor(u):
      log("#    ... top: %s = %s" % (blurb(u), blurb(v)))
    (u2, v2) = (u, v)

  else:
    m = get_link_in_block(u)  # a relation
    n = get_link_in_block(v)

    if v_up and u_up:
      if m is v and n is u:
        (u2, v2) = analyze_ladder(u_up, v_up)
        if monitor(u):
          log("#    ... match: %s = %s" % (blurb(u), blurb(v)))
        set_estimate(u, relation(EQ, v, note="match in block"))
        set_estimate(v, relation(EQ, u, note="match in block"))
        set_link(m, u); set_link(v, n) # ???????
        (u2, v2) = (u, v)
      else:
        (u2, v2) = analyze_ladder(u_up, v_up)
        if monitor(u):
          log("#    ... mismatch: %s <> %s" % (blurb(u), blurb(v)))
        set_estimate(u, relation(LT, v2, note="left mismatch"))
        set_estimate(v, relation(LT, u2, note="right mismatch"))

    elif not u_up:
      (u2, v2) = analyze_ladder(u, v_up)
      assert u is u2
      set_estimate(v, relation(LT, u, note="skip left at top"))
      if monitor(u):
        log("#    ... skip left: %s -> %s" % (blurb(u), blurb(v2)))

    else:
      assert not v_up
      (u2, v2) = analyze_ladder(u_up, v)
      assert v is v2
      set_estimate(u, relation(LT, v, note="skip right at top")) # clobbered
      if monitor(u):
        log("#    ... skip right: %s <- %s" % (blurb(u2), blurb(v)))

  # Theorem: u has an estimate.

  assert get_estimate(u, None), blurb(u)
  assert get_estimate(v, None), blurb(v)
  return (u2, v2)

# foo
def same_exemplar(u, v):
  if u is v: return True
  e1 = exemplar.get_exemplar(u)
  return e1 and e1 is exemplar.get_exemplar(v)

# -----------------------------------------------------------------------------
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

# wait.

def get_link_in_block(u):
  v = get_mutual_link(u, None)     # returns relation to node in AB
  # what if v is False (u is ambiguous)?
  if v and same_block(get_block(v), get_block(u)):
    return v
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

def optimize_relation(AB, u, rel):
  assert separated(u, rel.record)
  w = rel.record
  sup = find_cross_sup_rel(AB, u, w)
  if sup:
    # Make a copy of within-checklist relationship, retaining note
    assert sup.relationship & rel.relationship != 0
    return relation(sup.relationship, w, note=sup.note, span=sup.span)
  else:
    return rel
  
# Cannot use find_estimate because it relies on optimize_relation!
# Warning: rel can be in either checklist... we really just care
# about ship and note; we already know the target will be w

def find_cross_sup_rel(AB, u, w):
  # Are u/w in a child/parent configuration?
  assert u
  assert w
  assert separated(u, w)
  q = None
  # Option 1. u -> p = q ?= w
  rel = local_sup(AB, u)
  if rel:
    p = rel.record
    p_eq = get_equivalent(AB, p)
    if p_eq:
      q = p_eq.record
  if not q:
    # Option 2. u = z -> q ?= w
    u_eq = get_equivalent(AB, u)
    if u_eq:
      z = u_eq.record
      rel = local_sup(AB, z)
      if rel:
        q = rel.record
        # oops, could use rel directly without copying it!
  # Option 3. u ?= z = q -> w   -- cannot go from w to child q.
  return rel if q == w else None

def debug_block(AB, u, w, v0):
      log("# v0=u %s" % (v0 == u))
      log("# %s top %s" % (blurb(u), at_top_of_chain(AB, u)))
      log("# %s top %s" % (blurb(w), at_top_of_chain(AB, w)))

def at_top_of_chain(AB, u):
  sup = local_sup(AB, u)
  if not sup: return True
  b = get_block(u)
  # each peripheral has its own chain
  if is_empty_block(b): return True
  return not same_block(get_block(sup.record), b)

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

# -----------------------------------------------------------------------------
# ONLY USED WITHIN THIS FILE

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
  def do_cross_mrcas(AB):
    def traverse(x):            # arg in A, result in B
      u = AB.in_left(x)         # in AB
      exem = exemplar.get_exemplar(u)       # in AB
      if exem:
        (_, u1, v1) = exem
        v = v1 if separated(u, v1) else u1
        m = get_outject(v)
      else:
        m = BOTTOM                # identity for mrca
      for c in get_inferiors(x):  # c in A
        q = traverse(c)       # in B
        m0 = m
        m = simple.mrca(m, q)      # in B
      # Sanity checks
      b1 = get_block(u)
      assert (m == BOTTOM) == is_empty_block(b1)
      if m != BOTTOM:
        v = AB.in_right(m)

        assert separated(u, v)
        b2 = get_block(v)
        # AssertionError: ([122, 123, 124], [124])
        assert block_le(b1, b2), (list(b1), list(b2))

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
  b1 = get_block(AB.in_left(AB.A.top))
  b2 = get_block(AB.in_right(AB.B.top))
  assert b1 != BOTTOM_BLOCK
  assert b1 == b2

# For debugging

def show_exemplars(z, tag, AB):
  def foo(id):
    return blurb(exemplar.xid_to_record(AB, id, z))
  log("# theory: %s: {%s}" % (tag, ", ".join(map(foo, get_block(z)))))

# record -> list of exemplar ids

def exemplar_ids(AB, z):
  return list(get_block(z))

# The records in z's checklist corresponding to the exemplars
# in the block for z.

def exemplar_records(AB, z):
  return (exemplar.xid_to_record(AB, id, z) for id in exemplar_ids(AB, z))

def opposite_exemplar_records(AB, z):
  return (exemplar.xid_to_opposite_record(AB, id, z) for id in exemplar_ids(AB, z))

def get_block(x):
  return really_get_block(x, BOTTOM_BLOCK)

(really_get_block, set_block) = prop.get_set(prop.declare_property("block"))

# -----------------------------------------------------------------------------
# Implementation of blocks as Python sets of 'exemplars'.
# A 'block' is just a set of exemplars, implemented as ... a python set.
# The term 'block' comes from the mathematical treatment of partitions.

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
