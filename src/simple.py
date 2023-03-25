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
  if x1 == y1:                    # same accepted
    return compare_siblings(x, x1, y1, y)
  return compare_accepted_in_checklist(x, y)

def compare_accepted_in_checklist(x, y):
  (x_peer, y_peer) = find_peers(x, y)    # Decrease levels as needed
  if x_peer == y_peer:     # x <= x_peer = y_peer >= y
    if x_peer == x:     # x = x_peer = y_peer > y, so x > y
      sup = get_superior(y)
      if sup and x == sup.record:
        return reverse_relation(y, sup)
      else:
        return relation(GT, y, note="x = xp = yp > y")
    elif y_peer == y:
      # if x_peer is parent of y, then use sup relation
      sup = get_superior(x)
      if sup and y == sup.record:
        return sup
      else:
        return relation(LT, y, note="x < xp = yp = y")
    else:
      return relation(EQ, y)   # x = x_peer = y_peer = y
  else:
    if False:
      log("# peers: %s <= %s != %s >= %s" %
          (blurb(x), blurb(x_peer), blurb(y_peer), blurb(y)))
    return relation(DISJOINT, y, note="x < xp = yp > y")

# Compare x to y under the assumption that accepted(x) = accepted(y).
# x and y might be in different checklists.
# Requires review.

def compare_siblings(x, x1, y1, y):
  if x != x1 and y != y1:       # x <= x1 = y1 >= y
    if known_different_exemplars(x, y):
      return relation(NEQ, y, "sibling heterotypic synonyms", span=2)
    elif known_same_exemplar(x, y):                # Need fuzzy protonym compare
      return relation (EQ, y, "sibling homotypic synonyms", span=2)
    else:
      return relation(NOINFO, y, "sibling synonyms", span=2)
  elif x != x1:
    # LT -> heterotypic, EQ -> homotypic (?), SYNONYM -> unknown
    return relation(synonym_relationship(x, x1), y, "synonym of", span=1)
  elif y != y1:
    return relation(reverse_relationship(synonym_relationship(y, y1)),
                    y, "has synonym", span=1)
  else:
    # Provide match info ??
    return relation(EQ, y, MISSING if x == y else "equivalent", span=0)

# (ship, note) =     return relation(ship, y, note=note)


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

def descends_from(x, y):
  return get_level(x) < get_level(y)
