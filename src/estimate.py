#!/usr/bin/env python3

import property as prop
import checklist, workspace, simple, ranks

from util import log
from checklist import *
from workspace import *
from simple import BOTTOM, compare_per_checklist
from specimen import get_exemplar_info, get_exemplar, get_exemplar_id
from specimen import maybe_get_type_uf
from rcc5 import rcc5_symbol
from cross_mrca import compute_cross_mrcas, get_cross_mrca

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

# Blocks have not been computed at this point (but they could be)

# Can't make sense of the following:
# "Given a model M, let [u] be the least node in the opposite checklist
# that contains u.  Then we're interested in the minimal [u] taken over
# all models.  Does that work?  Is it unique?"

# Low-numbered ranks are closer to root

# For u in checklist 1, find smallest v in checklist 2 such that
# concepts C(u) <= C(v).
# Assume for now it's sufficient to find exemplars X(u) <= X(v).

# XS(u) = exemplar set of u

# Careful, u might be top or None.

def find_estimate(AB, u):       # u is in AB
  rel1 = get_central(AB, u)     # ancestor of u that has >0 exemplars
  # rel1 ascends from u up to top, looking for a node with a cross_mrca
  # u_central >= u
  u_central = rel1.record
  v = get_cross_mrca(u_central) # in AB, from B
  # n.b. XS(u) <= XS(v)
  #  but it could still be that v < u
  
  # Fall through with u_central, v
  # How does u relate to v?

  # loop: v starts at cross_mrca(u_central).  If it's too small it 
  # goes up (rootward) until it's >= to u...
  while True:
    # if u <= v: return v

    w = get_cross_mrca(v)   # w from A.  X(w) >= X(v) >= X(u)

    # w and u are in the same checklist, so are comparable.
    # Three cases:
    #   w < u    -- keep searching
    #   w = u    -- result should be EQ
    #   w > u    -- result should be LT

    if simple.simple_le(get_outject(u), get_outject(w)):
      # C(u) >= C(w)
      if u is w:                # u = v = w
        # TBD: might actually have C(u) > C(v).

        # Easiest improvement: look for same-exemplar-set ancestors of
        # v, and ancestors and descendants of u, and switch to vaguer
        # relationship or higher estimate if not unique.

        rel2 = relation(EQ, v, "reciprocal cross_mrca")
      else:                     # u < v <= w
        # X(u) < X(w), ergo C(u) < C(w)
        rel2 = relation(LT, v, "bigger cross_mrca")
      break

    #   w < u    -- keep searching
    sup = local_sup(AB, v)         # B.Tupaia montana
    if not sup:                    # v is top
      log("# %s has no estimate" % blurb(u))
      return None
    v = sup.record
    # iterate

  assert separated(u, rel2.record)
  return compose_relations(rel1, rel2) # u -> u_central -> v
  # TBD: convert LT to LE when synonym ?

def is_monotype(u):
  b = get_block(u)
  for c in get_inferiors(u):
    if get_block(c) == b:
      return True
  return False


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

# ----------------------------------------
# Extra stuff

# Least taxon q of B such that w < q 

def get_dominator(AB, w):  # not used as of 1/2026
  if isinA(AB, w):
    u = w
    v = get_estimate(w)
  else:
    u = get_estimate(w)
    v = w
  p = get_superior(u) if u is w else u
  q = get_superior(v) if v is w else v
  while True:
    rel = compare(AB, p, q)
    if rel.relationship == LT: return p
    elif rel.relationship == LE: return p
    elif rel.relationship == EQ: return q
    elif rel.relationship == GT: return q
    elif rel.relationship == GE: return q
    elif rel.relationship == OVERLAP:
      log("# Overlap: %s with %s" % (blurb(p), blurb(q)))
      # iterate
      p = get_superior(p)       # 'Break' the taxon
    else: assert False

# Least taxon z of AB with u= < z and v <= z

def mrca(AB, u, v):
  while True:
    m = AB.in_left(simple.mrca(get_outject(u), get_outject(get_estimate(v))))
    n = AB.in_right(simple.mrca(get_outject(v), get_outject(get_estimate(u))))
    rel = compare(AB, m, n)
    if   rel.relationship == EQ: return n
    elif rel.relationship == LT: return n
    elif rel.relationship == LE: return n
    elif rel.relationship == GT: return m
    elif rel.relationship == GE: return m
    elif rel.relationship == OVERLAP:
      # Does this terminate?  Yes, because simple.mrca always goes
      # rootward after an overlap.
      log("# Reducing mrca(%s, %s) to mrca(%s, %s)" %
          (blurb(x), blurb(y), blurb(m), blurb(n)))
      u = m
      v = n
    else: assert False

# -----------------------------------------------------------------------------

# Find a 'central' (non-peripheral) ancestor node (contains at least
# one exemplar) in SAME checklist.  Returns a Relation.

# Deprecated = see theory.py

def get_central(AB, u):
  u_central = u
  while u_central:
    if get_cross_mrca(u_central, None): break
    u_central = local_sup(AB, u_central).record
  # Make a Relation.
  if not u_central:
    log("# No central: %s" % blurb(u))
  if u_central is u:
    return relation(EQ, u_central)
  else:
    sup = local_sup(AB, u)
    if sup and sup.record is u_central:
      return sup                # u -> u_central = sup
    else:
      return relation(LT, u_central, note="get_central")
