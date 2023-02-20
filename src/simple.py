import property as prop
from checklist import *


# Same-tree relationships

#    x1 ? y1     'larger', 'up'
#   /       \
#  x         y   'smaller', 'down'

# Returns a Relation explaining justification

def simple_relationship(x, y):             # Within a single tree
  if x == y:
    return relation(EQ, y, note="identical")
  x1 = get_accepted(x)
  y1 = get_accepted(y)
  if x1 == y1:                    # same accepted
    (ship, note) = sibling_relationship(x, x1, y1, y)
    return relation(ship, y, note=note)

  (x2, y2) = find_peers(x1, y1)    # Decrease levels as needed
  if x2 == y2:     # x <= x2 = y2 >= y1
    if x2 == x1:     # x = x2 = y2 > y, so x > y
      return relation(GT, y, note="x = x2 = y2 > y")
    elif y2 == y1:
      return relation(LT, y, note="x < x2 = y2 = y")
    else:
      assert False  # can't happen
  else:
    return relation(DISJOINT, y, note="x <= x2 != y2 >= y")

def sibling_relationship(x, x1, y1, y): # x1 equivalent to y1
  if x1 != x and y1 != y:
    r1 = get_superior(x).relationship
    r2 = get_superior(y).relationship
    if r1 == SYNONYM or r2 == SYNONYM or r1 != r2:
      return (NOINFO, "sibling synonyms")
    elif r1 == LT:
      return (NEQ, "sibling heterotypic synonyms") # ?
    else:
      return (EQ, "sibling homotypic synonyms")
  elif x1 != x:
    r1 = get_superior(x).relationship
    return (r1, "synonym of") # could say heterotypic, homotypic, none
  elif y1 != y:
    r2 = get_superior(y).relationship
    return (reverse_relationship(r2), "has synonym")
  else:
    return (EQ, "identical" if x == y else "equivalent")

# Assumes already 'nipped'

def find_peers(x, y):
  if get_level(x, None) == None or get_level(x, None) == None:
    clog("# No level for one of these:", x, y)
    return NOINFO
  i = get_level(x)
  while i < get_level(y):
    y = get_parent(y)
  j = get_level(y)
  while get_level(x) > j:
    x = get_parent(x)
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

  return x1 == y

def simple_lt(x, y):
  return simple_le(x, y) and x != y

def in_same_tree(AB, x, y):
  return (AB.case(x, lambda x:1, lambda x:2) ==
          AB.case(y, lambda x:1, lambda x:2))

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
