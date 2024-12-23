import math

from util import log, MISSING

from checklist import get_superior, get_rank, get_canonical, blorb
#from parse import get_family

# Are x and y near enough that y might be a genus change of x or
# vice versa?  Called only when they are otherwise homotypic (same 
# epithet etc.)

# TBD: Use estimates if available!

def near_enough(x, y):          # In same or different checklists
  f = get_family(x)
  if not f:
    # log("# No family: %s" % blorb(x))
    return MISSING
  g = get_family(y)
  if not g:
    # log("# No family: %s" % blorb(y))
    return MISSING
  return f == g

# Not a match:
#  Maxomys pagensis (Miller, 1903) -> Macaca pagensis (G. S. Miller, 1903)
#  (families Muridae, Cercopithecidae)

def get_family(x):              # returns canonical name
  while get_rank(x, None) != 'family':
    sup = get_superior(x, None)
    if sup == None:
      # log("# No sup: %s %s" % (get_primary_key(x), blorb(x)))
      return None
    x = sup.record
  return get_canonical(x, None)



"""

NEAR_THRESHOLD = 30
FAR_THRESHOLD = 60

# distance is thresholded, so it only really matters whether it's small
# or large

# Obsolete (used only by linkage.py)

# For mammals, tip to root is expected to be about 13... 
# For 2M species, tip to root is expected to be about 20... 
# get_estimate returns smallest node in opposite checklist that is
# known to be <= the given node.

def compute_distance(u, v, get_estimate):
  assert separated(u, v)
  half1 = compute_half_distance(u, v, get_estimate)
  half2 = compute_half_distance(v, u, get_estimate)
  return int((half1 + half2)//2)

def compute_half_distance(u, v, get_estimate):
  assert separated(u, v)
  v1 = get_estimate(u)
  # If m is MRCA of u and v, then u <= v1 <= m >= v
  dist = distance_in_checklist(get_outject(v), get_outject(v1.record))
  return dist

"""


# Moved here from simple.py

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

