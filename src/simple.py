import math
import property as prop
from checklist import *


# Same-tree relationships

#    x1 ? y1     'larger', 'up'
#   /       \
#  x         y   'smaller', 'down'

# Returns a Relation explaining justification

def compare_per_checklist(x, y):             # Within a single checklist
  x1 = get_accepted(x)
  y1 = get_accepted(y)
  rel = compare_accepted_in_checklist(x1, y1)
  # Cf. theory.cross_compare
  if not x is x1:
    syn = get_superior(x)               # x <= x1
    rel = compose_relations(syn, rel)   # x <= x1 ? y1
  if not y is y1:
    syn = reverse_relation(y1, get_superior(y))  # y1 >= y
    rel = compose_relations(rel, syn)   # x <= x1 ? y1 >= y
  return rel

# x and y are accepted, no synonyming around.

def compare_accepted_in_checklist(x, y):
  (x_peer, y_peer) = find_peers(x, y)    # Decrease levels as needed
  if x_peer is y_peer:     # x <= x_peer = y_peer >= y
    peer = x_peer
    if x is y:
      return relation(EQ, y)   # x = peer = y
    elif peer is x:
      sup = get_superior(y)
      if sup and x is sup.record:
        return reverse_relation(y, sup)
      else:
        return relation(GT, y, note="x = xp = yp > y")
    elif peer is y:              # x < peer = y, so x < y
      # if x_peer is parent of y, then use sup relation
      sup = get_superior(x)
      if sup and y is sup.record:
        return sup
      else:
        return relation(LT, y, note="x < xp = yp = y")
    else:
      assert False
  else:
    return relation(DISJOINT, y, note="x < xp = yp > y")

# (ship, note) =     return relation(ship, y, note=note)


# Assumes synonyms already 'nipped' off

def find_peers(x, y):
  if get_level(x, None) == None or get_level(x, None) == None:
    clog("# No level for one of these:", x, y)
    return NOINFO
  i = min(get_level(x), get_level(y))
  while get_level(x) > i:
    x = get_parent(x)
  while get_level(y) > i:
    y = get_parent(y)
  return (x, y)

# MRCA within the same tree

def mrca(x, y):
  if x == y: return x
  if x == BOTTOM: return y
  if y == BOTTOM: return x

  x = get_accepted(x)
  y = get_accepted(y)
  (x, y) = find_peers(x, y)
  while not (x is y):
    x = get_parent(x)
    y = get_parent(y)
  return x

BOTTOM = None

def simple_le(x, y):
  # Is x <= y?  Scan from x upwards, see if we find y
  if x == y: return True
  y = get_accepted(y)
  x1 = x
  stop = get_level(y)

  while get_level(x1) > stop:
    x1 = get_parent(x1)    # y1 > previously, level(y1) < previously

  return x1 is y

def simple_lt(x, y):
  return simple_le(x, y) and not x is y

def simple_gt(x, y): return simple_lt(y, x)
def simple_ge(x, y): return simple_le(y, x)

(really_get_level, set_level) = prop.get_set(prop.declare_property("level", inherit=False))

# ----- Levels maintenance

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

# -----------------------------------------------------------------------------
# Distance computation

def distance_in_checklist(x1, x2):
  m = mrca(x1, x2)
  assert m
  return (distance_on_lineage(x1, m) +
          distance_on_lineage(x2, m))

def distance_on_lineage(x, m):
  assert simple_le(x, m)
  if x is m:
    return 0
  y = get_superior(x).record
  d = distance_on_lineage(y, m)
  if y is m: d -= 1
  return (distance_to_parent(x) + d)

def distance_to_parent(u):
  sup = get_superior(u)
  return lg(max(1,len(get_children(sup.record, ()))))

def lg(x):
  return math.log(x)//log2
log2 = math.log(2)

