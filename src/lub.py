#!/usr/bin/env python3

import math
import property as prop
import checklist, workspace, simple, exemplar

from util import log
from checklist import *
from workspace import *
from simple import BOTTOM, compare_per_checklist, compare_siblings
from linkage import really_find_links
from exemplar import equate_exemplars

import exemplar


def iteration(AB):
  # Links grow monotonically.  Ambiguities are recorded as False but that's OK.
  really_find_links(AB, compute_distance)

  # Exemplars grow monotonically, no need to overwrite.
  exemplar.analyze_exemplars(AB)   # does set_exemplar(...)

  # Cross-mrcas need to be replaced (LUBs get tighter/smaller).
  #  Easily overwritten.
  compute_cross_mrcas(AB)          # does set_cross_mrca(...)

  # Estimates depend highly on cross-mrcas.  Easily overwritten.
  find_estimates(AB)               # does set_estimate(...)

# -----------------------------------------------------------------------------

(get_estimate, set_estimate) = prop.get_set(prop.declare_property("estimate"))

def find_estimates(AB):
  counts = [0, 0]
  def findem(AB):
    def doit(AB):
      if True:
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
  log("-- Estimates: %s (=), %s (<)" % (int(counts[0]/2), counts[1]))

# Not sure about this

def find_estimate(AB, u):
  ship = EQ
  while True:                   # Skip over peripherals
    v = get_cross_mrca(u, None)
    if v != None:
      break
    sup = local_sup(AB, u)
    assert sup, (blurb(u1), blurb(u))
    u = sup.record
    ship = LT
  u2 = u
  v1 = v
  u1 = get_cross_mrca(v1)
  if simple.simple_gt(get_outject(u1), get_outject(u2)):
    return relation(LT, v, "goes up ladder")
  if not get_cross_mrca(u1) is v1:
    log("# Lose: u2 %s, v1 %s, u1 %s, v0 %s" %
        (blurb(u2), blurb(v1), blurb(u1), blurb(get_cross_mrca(u1))))
    assert False
  if u1 is u2:
    return relation(ship, v1, "bottom of ladder")

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
    u3 = get_mutual_link(v3, None)
    if (u3 and
        get_cross_mrca(u3, None) == v1 and
        get_cross_mrca(v3) == u1):
      # Candidate's (v3's) partner (u3) is in same block
      if u2 is u3:
        # Candidate's (v3's) partner (u3) belongs to our node u2
        if monitor(u3):
          log("# lub: Jackpot %s, %s" % (blurb(u2), blurb(v3)))
        equate_exemplars(u3, v3)
        note="linked"
        break
      if descends_from(u2, u3):
        # Answer is v3 (= v)
        if monitor(u3):
          log("# lub: Modestly %s, %s" % (blurb(u2), blurb(v3)))
        ship = LT
        note="goes up to linked"
        break

  return relation(ship, v, note)

def descends_from(u1, u2):
  return simple.descends_from(get_outject(u1), get_outject(u2))

# -----------------------------------------------------------------------------
# u assumed central

# Given u, find unique node v in opposite checklist such that v = u.
# x and y are in AB, not A or B.
# Returns a Relation or None.

def get_equivalent(AB, u):
  assert u
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
      if m != BOTTOM:
        v = AB.in_right(m)
        assert separated(u, v)
        set_cross_mrca(u, v)
      return m
    traverse(AB.A.top)
  do_cross_mrcas(AB)
  do_cross_mrcas(swap(AB))

(get_cross_mrca, set_cross_mrca) = \
  prop.get_set(prop.declare_property("cross_mrca"))

# -----------------------------------------------------------------------------

# distance is thresholded, so it only really matters whether it's small
# or large

# For mammals, tip to root is expected to be about 13... 
# For 2M species, tip to root is expected to be about 20... 

def compute_distance(u, v):
  assert separated(u, v)
  if get_estimate(u, None) and get_estimate(v, None):
    return int((compute_half_distance(u, v) +
                compute_half_distance(v, u))/2)
  else:
    return None

def compute_half_distance(u, v):
  # u < u1 <= (v1 < m > v)
  assert separated(u, v)
  v1 = get_estimate(u).record
  y = get_outject(v)
  y1 = get_outject(v1)
  m = simple.mrca(y1, y)
  dist = (distance_on_lineage(y1, m) +
          distance_on_lineage(y, m))
  return dist

def distance_on_lineage(u1, u2):
  if u1 == u2:
    return 0
  return (distance_to_parent(u1) +
          distance_on_lineage(get_superior(u1).record, u2))

def distance_to_parent(u):
  sup = get_superior(u)
  return lg(max(1,len(get_children(sup.record, ()))))

def lg(x):
  return math.log(x)/log2
log2 = math.log(2)

