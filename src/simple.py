from util import VERSION
import property as prop
from checklist import *
from rcc5 import EQ, GT, LT, DISJOINT, NOINFO

# Implementation of "Darwin Core semantics" from the mss

# Same-tree relationships

#    x1 ? y1     'larger', 'up'
#   /       \
#  x         y   'smaller', 'down'

NEWER = (VERSION > 1)

# Returns a Predicate

def compare_per_checklist(x, y):
  if NEWER:
    return compare_per_checklist_newer(x, y)
  else:
    return compare_per_checklist_older(x, y)

# Newer version, more species and congruences, fewer overlaps
# Easier to be seen as correct, but as of 1/23/2026 it's not exactly
# compatible with the mss
# Leads to 6% overall performance hit

def compare_per_checklist_newer(x, y):
  i = min(get_level(x), get_level(y))

  # x3 starts at x and runs up to its final value at level i
  x3 = x
  r13 = predicate(EQ, x)
  while get_level(x3) > i:
    sup = get_superior(x3, None)
    r13 = compose_predicates(r13, sup)  # extend r13
    x3 = sup.record

  # x4 starts at y and runs up to its final value at level i
  x4 = y
  r64 = predicate(EQ, x4)
  while get_level(x4) > i:
    sup = get_superior(x4, None)
    r64 = compose_predicates(r64, sup)  # extend r64
    x4 = sup.record

  r46 = reverse_predicate(x4, r64)
  
  # x3 and x4 are same 'level' (~ 'rank')
  if r13.record is r64.record:
    # r13 is EQ and/or r64 is EQ, so x <=> y
    ship = EQ                   # gets composed with < and/or >
  elif r13.relation != SYNONYM and r64.relation != SYNONYM:
    # accepted / accepted non-siblings
    # {x3 and descendants} disjoint from {x4 and descendants}
    ship = DISJOINT
  else:
    # accepted / synonym siblings or synonym / synonym siblings
    ship = NOINFO
  r34 = predicate(ship, r13.record, "bridge")
  return compose_predicates(compose_predicates(r13, r34), r46)

# Older version, numbers are as in mss as of 1/23/2026, ugly code, 
# dubious optimization, dubious correctness

def compare_per_checklist_older(x, y):
  # older version, matching 1/23/26 version of mss
  if x is y:
    return predicate(EQ, y) # x = peer = y
  (x_peer, y_peer) = find_peers(x, y)
  if x_peer is y_peer:     # x <= x_peer = y_peer >= y
    # They intersect, so < = >
    if x_peer is x:          # x = x_peer >= y, so x >= y (or > y, so x > y)
      y_sup = get_superior(y)
      if y_sup and x is y_sup.record:
        return reverse_predicate(y, y_sup)
      else:
        return predicate(GT, y, note="x = xsup = ysup > y")
    elif y_peer is y:          # x <= y_peer = y, so x <= y (or x <, so x < y)
      x_sup = get_superior(x)
      if x_sup and y is x_sup.record:
        return x_sup
      else:
        return predicate(LT, y, note="x < xsup = ysup = y")
    else:
      assert False, "shouldnt happen"
  else:
    # They do not intersect, or are sibling synonyms
    xsup = get_superior(x)
    ysup = get_superior(y)
    if (xsup and ysup and xsup.record is ysup.record
        and (not is_accepted(x) or not is_accepted(y))):
      return predicate(NOINFO, y, note="synonym ? sibling")
    return predicate(DISJOINT, y, note="x < xsup = ysup > y")

def find_peers(x, y):
  assert get_level(x, None) and get_level(x, None), "no levels"
  i = min(get_level(x), get_level(y))
  while get_level(x) > i:
    x = get_parent(x)
  while get_level(y) > i:
    y = get_parent(y)
  return (x, y)

# MRCA within the same tree

def mrca(x, y):
  if x == y: return x
  (x, y) = find_peers(x, y)
  while not (x is y):
    x = get_parent(x)
    y = get_parent(y)
  return x

# premature optimization of compare(x, y).relation == LE

def simple_le(x, y):
  # Is x <= y?  Scan from x upwards, see if we find y
  if x == y: return True
  assert not is_accepted(y)
  stop = get_level(y)
  while get_level(x) > stop:
    x = get_parent(x)    # y1 > previously, level(y1) < previously
  return x is y

def simple_lt(x, y):
  if x == y: return False
  return simple_le(x, y)

def simple_gt(x, y): return simple_lt(y, x)
def simple_ge(x, y): return simple_le(y, x)

# Necessary or possible intersection... ?
# So far we've been treating RCC-5 as disjunctive

def simple_intersect(x, y):
  pred = compare_per_checklist(x, y)
  return ((pred.relation & INTERSECT) != 0 and # could
          (pred.relation & DISJOINT) == 0)     # must

# ----- Levels maintenance

(really_get_level, set_level) = prop.get_set(prop.declare_property("level", inherit=False))

# 2nd result is True if x is a synonym, False otherwise
#
def get_level(x, default=None):
  i = really_get_level(x, None)
  if i == None:
    set_level(x, "cycle")
    p = get_parent(x)
    if not p:
      i = 1                 # shouldn't happen.  top is level 1
    else:
      assert p != x
      lev = get_level(p, None)
      assert lev != "cycle", \
        ("**** cycle detected while traversing ancestor chain: %s->%s" %
            (blurb(x), blurb(p)))
      i = lev + 1
    set_level(x, i)
    #if monitor(x):
    #  log("# level(%s) = %s" % (blurb(x), i))
  return i

def get_parent(x):
  rp = get_superior(x, None)
  if rp:                        # and is_accepted(x) ????
    return rp.record
  else: return None

def ensure_levels(S):
  def cache(x, n):
    set_level(x, n)
    for c in get_inferiors(x):
      # This isn't right -- get_level works better
      cache(c, n+1)
  cache(S.top, 1)

def descends_from(x, y):
  return get_level(y) < get_level(x)

