#!/usr/bin/env python3

import property as prop
import checklist, workspace, simple

from util import log
from checklist import *
from workspace import *
from simple import BOTTOM, compare_per_checklist, compare_siblings
from typify import equate_typifications, get_typification_record, \
  get_link

# -----------------------------------------------------------------------------

(get_estimate, set_estimate) = prop.get_set(prop.declare_property("estimate"))

def find_estimates(AB):
  # Cross-mrcas might need to be replaced 
  #  (LUBs get tighter/smaller on 2nd pass).
  #  Easily overwritten.
  compute_cross_mrcas(AB)          # does set_cross_mrca(...)

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
  while True:
    # Skip over peripherals, locating non-peripheral u2
    v = get_cross_mrca(u2, None)
    if v != None:
      break
    sup = local_sup(AB, u2)
    if not sup:
      # u2 is top
      return relation(ship, u2, "estimate = top")
    assert sup, blurb(u2)       # not at root
    u2 = sup.record
    ship = LT

  # v is lower bound on u2's estimate (and therefore u's estimate)

  # Let a "chain" be a maximal sequence of nodes all having the same
  # cross_mrca (or exemplar set).

  # If  u is part of a chain u1 < ... u ... < u4
  # and cross_mrca(u) = v1,
  # and we have a chain      v1 < ...   ... < v4, then...

  # with v1 = cross_mrca(u), u1 = cross_mrca(v1),    ***** FALLACY.
  # then if   u <= u1   then the answer is v1 ... v4
  #      if   u > u1    then the answer is v1. ??

  v1 = v
  # could have v1 >= u2 (if v1 has more or same exemplars as u2)
  # or v1 < u2 (if same exemplars but u2 has additional non-exemplars)

  u1 = get_cross_mrca(v1)
  # necessarily u1 >= v1     *** DON'T THINK SO.

  if simple.simple_gt(get_outject(u1), get_outject(u2)):
    return relation(LT, v, "goes up ladder")
  # necessarily u1 <= u2... but it can't be <, can it?

  v0 = get_cross_mrca(u1)
  if not v0 is v1:
    # could have u1 ~ v1 <= u1 <= v0... very unusual
    log("# Lose: u2 %s, v1 %s <= u1 %s <= v0 %s" %
        (blurb(u2), blurb(v1), blurb(u1), blurb(v0)))
    return relation(ship, v1, "lose")
  if u1 is u2:
    return relation(ship, v1, "bottom of ladder")
  if monitor(u1) or monitor(u2) or monitor(v1):
    set_scene(AB, u, u1, v1, u2)

  # Scan upwards from v1 looking for a way back to the u chain, either
  # by a linked name or by hitting the top of the chain.

  while True:

    # u top: B  --> ?
    #        A <--> A :v top

    v_sup = local_sup(AB, v)
    assert v_sup
    # We're at top of chain if v's parent is not in u1/v1 block
    v = v_sup.record              # Try next higher v

    if not get_cross_mrca(v) is u1:
      note="v off top of chain"
      ship = LT                 # ?
      break

    # Search for links by name
    v3 = v

    # See if chains are linked here
    u3 = get_link(v3, None)
    if (u3 and
        get_cross_mrca(u3, None) == v1 and
        get_cross_mrca(v3) == u1):
      # Candidate's (v3's) partner (u3) is in same block
      if u2 is u3:
        # Candidate's (v3's) partner (u3) belongs to our node u2
        if monitor(u3):
          log("# estimate: Jackpot %s, %s" % (blurb(u2), blurb(v3)))
        equate_typifications(u3, v3)
        note="linked"
        break
      if simple.simple_lt(u2, u3):
        # Answer is v3 (= v)
        if monitor(u3):
          log("# estimate: Modestly %s, %s" % (blurb(u2), blurb(v3)))
        ship = LT
        note="goes up to linked"
        break

  return relation(ship, v, note)

# is u correct here?

def set_scene(AB, u, u1, v1, u2):
  us = []
  u = u1
  while get_cross_mrca(u) is v1:
    us.append(blurb(u))
    u = local_sup(AB, u).record
  vs = []
  v = v1
  while get_cross_mrca(v) is u1:
    vs.append(blurb(v))
    v = local_sup(AB, v).record
  u2p = local_sup(AB, u2).record
  if u is u2:
    log("# u = u2 %s | %s" % (blurb(u), blurb(u2p)))
  else:
    log("# u %s --> u2 %s | %s" % (blurb(u), blurb(u2), blurb(u2p)))
  log("#   v1...: %s | %s" % ("; ".join(vs), blurb(v)))
  log("#   u1...: %s | %s" % ("; ".join(us), blurb(u)))

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
      exem = get_typification_record(u)       # exemplar record (not uf)
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
