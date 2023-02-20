import property as prop
from checklist import *


# Same-tree relationships

#    x1 ? y1     'larger', 'up'
#   /       \
#  x         y   'smaller', 'down'

# Returns a Relation explaining justification

def simple_relationship(x0, y0):             # Within a single tree
  (x, xsyn) = nip_synonym(x0)
  (y, ysyn) = nip_synonym(y0)
  if x == y:
    if xsyn and ysyn:
      (ship, note) = typical(x0, y0, NEQ, NOINFO, "sibling %ss")
    elif xsyn:
      (ship, note) = typical(x0, y0, LT, SYNONYM, "%s of")
    elif ysyn:
      (ship, note) = typical(x0, y0, GT, MYNONYS, "has %s")
    else:
      (ship, note) = (EQ, "identical")
    return relation(ship, y, None, note)
    # return(ship, y, None, note)
  (x1, y1) = find_peers(x, y)    # Decrease levels as needed
  if x1 == y1:     # x <= x1 = y1 >= y
    if x1 == x:     # x = x1 = y1 > y, x > y
      return relation(GT, y, None, "x = x1 = y1 > y")
    elif y1 == y:
      return relation(LT, y, None, "x < x1 = y1 = y")
    else:
      assert False  # can't happen
  else:
    return relation(DISJOINT, y, None, "x <= x1 != y1 >= y")

# If x and y homotypic then EQ, else if heterotypic LT, else NOINFO
# cf. theory.py
def typical(x, y, het_ship, no_ship, note):
  t1 = get_tipe(x, None)
  t2 = get_tipe(y, None)
  if t1 and t2:
    if t1 == t2:
      return (EQ, note % "homotypic synonym")
    else:
      return (het_ship, note % "heterotypic synonym")
  else:
    return (no_ship, note % "synonym")

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

  (x, _) = nip_synonym(x)
  (y, _) = nip_synonym(y)
  (x, y) = find_peers(x, y)
  while not (x is y):
    x = get_parent(x)
    y = get_parent(y)
  return x

BOTTOM = None

def simple_le(x, y):
  # Is x <= y?  Scan from x upwards, see if we find y
  if x == y: return True
  (y, _) = nip_synonym(y)
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
def nip_synonym(x):
  xsyn = False
  while True:
    sup = get_superior(x, None)
    if not sup or is_accepted(x):
      break
    if sup.relationship != EQ:  # ?
      xsyn = True
    x = sup.record
  assert get_level(x, None) != None, blurb(x)
  return (x, xsyn)

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
