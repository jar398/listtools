#!/usr/bin/env python3

import types
import property as prop, checklist, workspace

from util import log
from checklist import *
from workspace import *
from match_records import match_records

def get_left_superior(AB, z):
  #if not z: return None
  x = get_left_persona(AB, z)
  if not x: return None
  ship = AB.case(x,
                 lambda x:get_superior(x, None),
                 lambda y:None)
  if ship:
    return Related(ship.relation, AB.in_left(ship.record), ship.status, ship.note)
  return None

def get_right_superior(AB, z):
  #if not z: return None
  y = get_right_persona(AB, z)
  if not y: return None
  ship = AB.case(y,
                 lambda x:None,
                 lambda y:get_superior(y, None))
  if ship:
    return Related(ship.relation, AB.in_right(ship.record), ship.status, ship.note)
  return None

def get_left_persona(AB, z):
  assert isinstance(z, prop.Record)
  assert get_outject(z, None)
  p = AB.case(z,
              lambda x:z,
              lambda y:equal_partner(z))
  if p: assert get_outject(p, None)
  return p

def get_right_persona(AB, z):
  return AB.case(z,
                 lambda x:equal_partner(z),
                 lambda y:z)       # B case

def equal_partner(z):
  ship = get_equated(z, None)
  if ship:
    assert get_outject(ship.record, None)
    return ship.record
  else: return None

def isinB(AB, z):
  return AB.case(z, lambda x: False, lambda x: True)

def equated(x, y):              # Are x and y equated?
  if x == y: return True
  z = get_equated(x, None)
  return z and z.record == y

def mtrmset_lt(AB, x, y):
  rel = mtrmset_relation(AB, x, y)
  return rel == LT    # or LE for synonyms ...?

# -----------------------------------------------------------------------------
# Precompute MTRM sets.  Assumes compute_mtrmsets has already run.

def mtrmset_relation(AB, x, y):
  assert isinstance(x, prop.Record)
  assert isinstance(y, prop.Record)
  if in_same_tree(AB, x, y):
    rel = simple_relation(AB.case(x, lambda x:x, lambda y:y),
                          AB.case(y, lambda x:x, lambda y:y))
    return rel
  else:
    if is_toplike(x):
      if is_toplike(y): return EQ
      else: return GT
    elif is_toplike(y): return LT
    #log("# Compare %s to %s, neither is top" % (blurb(x), blurb(y)))
    return relate_mtrmsets(get_mtrmset(x), get_mtrmset(y))

def compute_mtrmsets(AB):
  def traverse(x, in_left):
    e = null_mtrmset
    for c in get_inferiors(x):
      e = join_mtrmsets(e, traverse(c, in_left))
    if is_empty(e):
      z = in_left(x)
      w = get_equated(z, None)  # Does z represent an MTRM ?
      if w:
        e = eta(z, w.record)
        set_mtrmset(x, e)
    else:
      set_mtrmset(x, e)
    #log("# mtrmset(%s) = %s" % (blurb(x), e))
    return e
  traverse(AB.B.top, AB.in_right)
  traverse(AB.A.top, AB.in_left)

# Mtrmsets are sets in this instance, and aren't mtrmsets

null_mtrmset = set()
def join_mtrmsets(e1, e2): return e1 | e2
def eta(x, y): return {min(x.id, y.id)}
def is_empty(e): return len(e) == 0

def relate_mtrmsets(e1, e2):  # ASSUMES OVERLAP
  if e1 == e2: return COMPARABLE
  elif e1.issubset(e2): return LT
  elif e2.issubset(e1): return GT
  elif e1.isdisjoint(e2): return DISJOINT
  else: return CONFLICT

(get_mtrmset, set_mtrmset) = prop.get_set(prop.get_property("mtrmset"))

# -----------------------------------------------------------------------------
# Find mutual tipward matches

def find_MTRMs(AB):
  find_TRMs(AB.A, AB.in_left, set_trm)
  counter = [1]
  def set_if_mutual(z, m):      # z in B
    set_trm(z, m)           # Nonmutual, might be of interest
    z2 = m.record           # z2 in A
    m2 = get_trm(z2, None)      # tipward match to/in A
    if m2 and m2.relation == EQ and m2.record == z:
      #if monitor(z): log("# MTRM: %s :=: %s" % (blurb(z), blurb(z2)))
      propose_equation(AB, z2, z,
                       "MTRM;%s;%s" % (m2.note, m.note))
  find_TRMs(AB.B, AB.in_right, set_if_mutual)    # of AB.flip()

(get_trm, set_trm) = prop.get_set(prop.get_property("TRM"))

def find_TRMs(A, in_left, set_trm):
  ensure_inferiors_indexed(A)
  def traverse(x):
    seen = False
    for c in get_inferiors(x):
      seen = traverse(c) or seen
    if not seen and not is_top(x):
      z = in_left(x)
      m = get_match(z)
      if m and m.relation == EQ:
        #if monitor(z): log("# TRM: %s = %s" % (blurb(z), blurb(m.record)))
        set_trm(z, m)
        seen = z
    return seen
  traverse(A.top)

def loser():
  if False: yield "lose"

# -----------------------------------------------------------------------------
# Alignment building...

# Propose that x = y (= ry.record)

def propose_equation(AB, x, y, why_equiv):
  # Polarize
  assert isinB(AB, x) != isinB(AB, y)
  (x, y) = AB.case(x, lambda xx: (x, y), lambda yy: (y, x))

  # Record reason for the equation
  set_equated(y, Related(EQ, x, "equivalent", why_equiv))
  set_equated(x, Related(EQ, y, "equivalent", why_equiv))

  # Deprecate the non-priority record of the two
  set_superior(x, Related(SYNONYM, y, "equivalent", why_equiv))

# Propose that rs.record should be the parent (superior) of z

def propose_superior(AB, z, rs, status, note):
  assert rs
  assert isinstance(z, prop.Record), blurb(z)
  assert isinstance(rs, Related)
  assert rs.relation == ACCEPTED or rs.relation == SYNONYM
  set_superior(z, Related(rs.relation,
                          rs.record,
                          status or rs.status,    # accepted, etc
                          note))  # explanation of why < or <=

# -----------------------------------------------------------------------------
# Same-tree relations

#    x1 ? y1     'larger'
#   /       \
#  x         y   'smaller'

def simple_relation(x, y):             # Within a single tree
  (x1, y1) = find_peers(x, y)    # Decrease levels as needed
  if x1 == y1:     # x <= x1 = y1 >= y
    if x == y:
      return EQ
    elif x1 == x:     # x = x1 = y1 > y, x > y
      return GT
    elif y1 == y:
      return LT
    else:
      assert False
  else:
    # !!! FIX FOR SYNONYMS
    return DISJOINT

def find_peers(x, y):
  i = get_level(x)
  while i < get_level(y):
    y = get_superior(y).record
  j = get_level(y)
  while get_level(x) > j:
    x = get_superior(x).record
  return (x, y)

(get_level, set_level) = prop.get_set(prop.get_property("level", inherit=False))

def simple_le(x, y):
  # Is x <= y?  Scan from x upwards, see if we find y
  x1 = x
  stop = get_level(y)

  while get_level(x1) > stop:
    x1 = get_superior(x1).record    # y1 > previously, level(y1) < previously

  return x1 == y

def simple_lt(x, y):
  return simple_le(x, y) and x != y

def in_same_tree(AB, x, y):
  return (AB.case(x, lambda x:1, lambda x:2) ==
          AB.case(y, lambda x:1, lambda x:2))

def ensure_levels(S):
  def cache(x, n):
    set_level(x, n)
    for c in get_inferiors(x):
      cache(c, n+1)
  cache(S.top, 1)
