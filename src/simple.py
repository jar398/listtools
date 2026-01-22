import property as prop
from checklist import *
from rcc5 import EQ, GT, LT, DISJOINT, NOINFO


# Same-tree relationships

#    x1 ? y1     'larger', 'up'
#   /       \
#  x         y   'smaller', 'down'

# Returns a Relation explaining justification

def compare_per_checklist(x, y):             # Within a single checklist
  if x is y:
    return relation(EQ, y) # x = peer = y
  (x_peer, y_peer) = find_peers(x, y)    # Decrease levels as needed
  if x_peer is y_peer:     # x <= x_peer = y_peer >= y
    # They intersect, so < = >
    if x_peer is x:          # x = x_peer >= y, so x >= y (or > y, so x > y)
      y_sup = get_superior(y)
      if y_sup and x is y_sup.record:
        return reverse_relation(y, y_sup)
      else:
        return relation(GT, y, note="x = xsup = ysup > y")
    elif y_peer is y:          # x <= y_peer = y, so x <= y (or x <, so x < y)
      x_sup = get_superior(x)
      if x_sup and y is x_sup.record:
        return x_sup
      else:
        return relation(LT, y, note="x < xsup = ysup = y")
    else:
      assert False
  else:
    # They do not intersect, or are sibling synonyms
    xsup = get_superior(x)
    ysup = get_superior(y)
    if (xsup and ysup and xsup.record is ysup.record
        and (not is_accepted(x) or not is_accepted(y))):
      return relation(NOINFO, y, note="synonym ? sibling")
    return relation(DISJOINT, y, note="x < xsup = ysup > y")

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

# =============================================================================

# Following is a currently unused rewrite of compare_per_checklist
# that I hope is seen, more easily than the above, as an
# implementation of the Darwin Core semantics rules in the paper.
# I hope it's not any slower...

def factor_path(x, y):          # x1 to x6
  i = min(get_level(x), get_level(y))

  # x3 starts at x1 and runs up to its final value at bridge
  x3 = x
  r13 = relation(EQ, x)
  while get_level(x3) > i:
    sup = get_superior(x3, None)
    r13 = compose_relations(r13, sup)  # extend r13
    x3 = sup.record

  # x4 starts at x6 and runs up to its final value at bridge
  x4 = y
  r64 = relation(EQ, x4)
  while get_level(x4) > i:
    sup = get_superior(x4, None)
    r64 = compose_relations(r64, sup)  # extend r64
    x4 = sup.record

  r46 = reverse_relation(x4, r46)
  # Cases
  # lineage, 0 synonyms - EQ
  # lineage, 1 synonym  - <= (r34 can be EQ)
  # lineage, 2 synonyms - cannot happen.
  # separated, 0 synonyms - DISJOINT
  # separated, 1 synonym  - NOINFO
  # separated, 2 synonyms - NOINFO

  # x3 and x4 are same 'level' (~ 'rank')
  if r13.record is r64.record:
    # r13 is EQ and/or r64 is EQ, so x <=> y
    ship = EQ                   # gets composed with < and/or >
  elif r13.relationship != SYNONYM and r64.relation != SYNONYM:
    # accepted / accepted non-siblings
    # {x3 and descendants} disjoint from {x4 and descendants}
    ship = DISJOINT
  else:
    # accepted / synonym siblings or synonym / synonym siblings
    ship = NOINFO
  r34 = relation(ship, r13.record, "bridge")
  return (r13, r34, r46)

def compare_per_checklist_draft(x, y):
  (r13, r34, r46) = factor_path(x, y)
  r14 = compose_relations(r13, r34)
  return compose_relations(r14, r46)

