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

# Can't make sense of the following:
# "Given a model M, let [u] be the least node in the opposite checklist
# that contains u.  Then we're interested in the minimal [u] taken over
# all models.  Does that work?  Is it unique?"

# Low-numbered ranks are closer to root

# For u in checklist 1, find smallest v in checklist 2 such that u <= v
# Careful, u might be top or None.

def find_estimate(AB, u):       # u is in AB
  rel1 = get_central(AB, u)
  # rel1 ascends from u up to top, looking for a node with a cross_mrca
  u_central = rel1.record
  v = get_cross_mrca(u_central) # in AB
  # n.b. v >= every exemplar in u
  #  but it could still be < u (or > u)

  # Fall through with u_central, v
  # How does u relate to v?

  u_rank_n = ranks.ranks_dict.get(get_rank(u, 0))

  # v starts at cross_mrca u and goes up (rootward)
  while True:
    w = get_cross_mrca(v)   #get_block(v)    # or get_cross_mrca(v)

    if simple.simple_le(get_outject(u), get_outject(w)):
      if u is w:                # u = v = w
        rel2 = relation(EQ, v, "reciprocal cross_mrca")
      else:                     # u < v
        rel2 = relation(LT, v, "bigger cross_mrca")
      break

    sup = local_sup(AB, v)         # B.Tupaia montana
    if not sup:                    # v is top
      log("# %s has no estimate" % blurb(u))
      return None
    v = sup.record
    # iterate

  assert separated(u, rel2.record)
  return compose_relations(rel1, rel2) # u -> u_central -> v
  # TBD: convert LT to LE when synonym ?

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
  count = [0]
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
        count[0] += 1
      return m
    traverse(WS.A.top)
  do_cross_mrcas(AB)
  do_cross_mrcas(swap(AB))
  log("# %s cross_mrcas set" % count[0])

(get_cross_mrca, set_cross_mrca) = \
  prop.get_set(prop.declare_property("cross_mrca"))

# -----------------------------------------------------------------------------

# Find a 'central' (non-peripheral) ancestor node (contains at least
# one exemplar) in SAME checklist.  Returns a Relation.

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

# --------------------

# Convenience.  Phase this out?  Or rename it?

def get_link(u, default=-19):
  uf = maybe_get_type_uf(u, None)
  if uf:
    (_, u2, v) = uf.payload()
    return v if (v and separated(u, v)) else u2
  return None
