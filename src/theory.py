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
    return relation(ship.relationship, AB.in_left(ship.record), ship.status, ship.note)
  return None

def get_right_superior(AB, z):
  #if not z: return None
  y = get_right_persona(AB, z)
  if not y: return None
  ship = AB.case(y,
                 lambda x:None,
                 lambda y:get_superior(y, None))
  if ship:
    return relation(ship.relationship, AB.in_right(ship.record), ship.status, ship.note)
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

# -----------------------------------------------------------------------------
# Decide about relations between two trees

# Returns a pair (overlap, outlier)
#   overlap is in y and x {< = > ><} y (but not !)
#   outlier is in y and x {< >< !} y (but not {> =}) x ⋧ y
#   either can be None

def analyze(AB, x, y):
  (over, out) = (None, None)
  for d in get_inferiors(y):
    (over2, out2) = analyze(x, d)
    over = over or over2; out = out or out2
    if over and out: return (over, out)
  # TBD: if over and out are both present, and one of the two is a synonym ...
  if over or out:      # Exclude peripherals
    return (over, out)
  m = AB.record_match(y)
  j = AB.join(m, y)
  return (j, None) if m else (None, j)

def outlier(AB, x, y):
  (over, _) = analyze(x, y)
  return out

# Cache this.
# Satisfies: if y = shadow(x) then x ~<= y, and if x' = shadow(y),
# then if x' <= x then x ~= y, or if x' > x then x > y.
#
# i.e.   x ≲ y ≲ x' ≤ x   →   x ≅ y   (same 'block')
# i.e.   x ≲ y ≲ x' > x   →   x < y

def get_estimate(AB, x, other):
  m = functools.reduce(mrca,
                       (get_shadow(c) for c in get_inferiors(x)),
                       None)
  if m == None:
    return record_match(x, other) or None
  return m

# RCC5 decision procedure, where x and y are in different sources.

def cross_ge(AB, x, y):
  return not outlier(AB, x, y)

def cross_eq(AB, x, y):
  # see if they're peers ??
  return cross_ge(x, y) and cross_ge(y, x)

def gt(AB, x, y):
  return cross_ge(x, y) and not cross_ge(y, x)

# -----------------------------------------------------------------------------
# Precompute MTRM sets.  Assumes compute_mtrmsets has already run.

def mtrmset_lt(AB, x, y):
  rel = mtrmset_relationship(AB, x, y)
  return rel == LT    # or LE for synonyms ...?

def mtrmset_relationship(AB, x, y):
  assert isinstance(x, prop.Record)
  assert isinstance(y, prop.Record)
  if in_same_tree(AB, x, y):
    rel = simple_relationship(AB.case(x, lambda x:x, lambda y:y),
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
# Alignment building...

# Propose that x = y (= ry.record)

def propose_equation(AB, x, y, why_equiv):
  # Polarize
  assert isinB(AB, x) != isinB(AB, y)
  (x, y) = AB.case(x, lambda xx: (x, y), lambda yy: (y, x))

  # Record reason for the equation
  # This is a kludge - I'm encoding the A/B polarity in the relationship
  set_equated(y, relation(EQ, x, "equivalent", why_equiv))
  set_equated(x, relation(SYNONYM, y, "equivalent", why_equiv))

  # Deprecate the non-priority record of the two
  # This doesn't work; y might already be a synonym
  # set_superior(x, relation(SYNONYM, y, "equivalent", why_equiv))

# Propose that rs.record should be the parent (superior) of z

def propose_superior(AB, z, rs, status, note):
  assert rs
  assert isinstance(z, prop.Record), blurb(z)
  assert isinstance(rs, Relative)
  assert rs.relationship == ACCEPTED or rs.relationship == SYNONYM
  set_superior(z, relation(rs.relationship,
                           rs.record,
                           status or rs.status,    # accepted, etc
                           note))  # explanation of why < or <=

# -----------------------------------------------------------------------------
# Same-tree relations

#    x1 ? y1     'larger'
#   /       \
#  x         y   'smaller'

def simple_relationship(x, y):             # Within a single tree
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

# MRCA within the same tree

def mrca(x, y):
  if x == BOTTOM: return y
  if y == BOTTOM: return x
  (x, y) = find_peers(x, y)
  while not (x is y):
    x = get_superior(x)
    y = get_superior(y)
  return x

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
