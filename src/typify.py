
import math
import util
import property as prop
import simple

from util import log, UnionFindable
from rcc5 import EQ, NEQ
from checklist import get_parts, monitor, get_superior, get_children, \
  is_accepted, blurb
from workspace import separated, get_outject, get_workspace
from workspace import isinA, isinB

def find_typifications(AB, subproblems, get_pre_estimate):
  # This sets the 'typification_uf' property of ... some ... records.
  log("# Finding some typifications")

  n = 1
  for (key, (us, vs)) in subproblems.items():
    if n % 1000 == 0:
      log("# Subproblem %s %s %s %s %s" % (n, len(us), len(vs), blurb(us[0]), blurb(vs[0])))
    n += 1
    try_in_same_checklist(us)
    try_in_same_checklist(vs)
    for i in range(0, len(us)):
      u = us[i]
      # If it's a synonym, see if it matches any accepted in same checklist?
      for j in range(i, len(vs)):
        v = vs[j]
        # classify as A1, A2, A3 (HETEROTYPIC, REVIEW, HOMOTYPIC)
        # **** COMPUTE DISTANCE if 2nd pass ****
        dist = compute_distance(u, v, get_pre_estimate)
        score = compute_score(u, v, dist)
        if homotypic_score(score):
          equate_typifications(u, v)

  equate_typifications(AB.in_left(AB.A.top),
                       AB.in_right(AB.B.top))

def try_in_same_checklist(us):
  for i in range(0, len(us)):
    u1 = us[i]
    for j in range(i+1, len(us)):
      u2 = us[j]
      dist = distance_in_checklist(get_outject(u1), get_outject(u2)) * 2
      score = compute_score(u1, u2, dist)
      if homotypic_score(score):
        equate_typifications(u1, u2)

# u and v are in workspace but may or may not be from same checklist

def equate_typifications(u, v):     # opposite checklists. u might be species
  if u != v:
    equate_typification_ufs(get_typification_uf(u), get_typification_uf(v))
  return u

def equate_typification_ufs(uf, vf):
  uf = uf.find()
  vf = vf.find()
  (i1, u1, v1) = uf.payload()
  (i2, u2, v2) = vf.payload()
  assert i1 == None or i2 == None or i1 == i2
  assert u1 or v1
  assert u2 or v2
  # What if ambiguity, i.e. pick_better_record returns False?
  ef = uf.absorb(vf)          # ef happens to be uf
  r = ef.payload()
  r[1] = pick_better_record(u1, u2)
  r[2] = pick_better_record(v1, v2)
  return ef

# Only workspace nodes have uf records

def get_typification_uf(u):
  #log("# Thinking about typification for %s" % blurb(u))
  probe = really_get_typification_uf(u, None)
  if probe: return probe
  AB = get_workspace(u)
  exem = [None, u, None] if isinA(AB, u) else [None, None, u]
  uf = UnionFindable(exem)
  assert exem[1] or exem[2]
  set_typification_uf(u, uf)
  return uf

# Union-find typification nodes start out with no xid, then get xid
# when exemplars are identified

(really_get_typification_uf, set_typification_uf) = \
  prop.get_set(prop.declare_property("typification_uf"))


# Returns typification record (xid, u, v) or None.

def get_exemplar(z):
  uf = really_get_typification_uf(z, None)
  if uf:
    uf = uf.find()
    r = uf.payload()
    (xid, u, v) = r
    if u and v:
      if xid == None:
        # Create exemplar id (for set operations) on demand
        ws = get_workspace(u)
        xid = len(ws.exemplar_ufs)
        r[0] = xid
        ws.exemplar_ufs[xid] = uf
        #log("# Exemplar %s: (%s) <-> (%s)" % (xid, blurb(u), blurb(v)))
      return r
  return None
  
def get_typification_record(z):
  uf = really_get_typification_uf(z, None)
  if uf:
    r = uf.payload()
    (xid, u, v) = r
    if u and v:
      return r
  return None

# This chooses the 'best' record to represent a given exemplar

def pick_better_record(v1, v2):
  if v2 is None: return v1
  if v1 is None: return v2                   # One side, not nec. reciprocal
  if v1 is v2: return v1
  if v1 == False: return v1   # Propagate ambiguity block
  if v2 == False: return v2
  y1 = get_outject(v1)
  y2 = get_outject(v2)
  if is_accepted(y2) and not is_accepted(y1):
    # silently improve the choice
    return v2
  if not is_accepted(y2) and is_accepted(y1):
    # silently keep choice
    return v1
  m = simple.mrca(y1, y2)
  assert not y1 is y2
  if m is y1: return v2
  if m is y2: return v1
  # See whether v2 is an improvement over v1
  parts1 = get_parts(v1)
  parts2 = get_parts(v2)
  if len(parts2.middle) > 0 and len(parts1.middle) == 0:
    # Longer name (subspecies overrides species).
    return v2          # Go tipward (to v2)
  if len(parts2.middle) == 0 and len(parts1.middle) > 0:
    # Shorter name (species overrides subspecies).
    return v1                    # Keep tipward (v1)
  if is_accepted(y1) and is_accepted(y2):
    log('')
    log("# %s" % (parts1,))
    log("# %s" % (parts2,))
    log("# middles '%s' '%s'" % (parts1.middle,  parts2.middle))
    log("# %s, %s <= %s" % (blurb(y1), blurb(y2), blurb(m)))
    log("# typify: Ambiguous: %s & %s" %
        (blurb(v1), blurb(v2)))
    return False  # Ambiguous
  return v1       # arbitrary synonym choice; don't want ambiguous



HOMOTYPIC   = EQ           # Hmm, not sure about this
HETEROTYPIC = NEQ
REVIEW      = HOMOTYPIC | HETEROTYPIC

# For combining negative and positive scores:
NEUTRAL = 0


# Compare potentially contypic taxa.
# distance is distance (perhaps estimated) between u and v in combined model.
# If distance (in tree) is close, give a pass on genus mismatch.

def compute_score(u, v, distance=None):
  score = compute_parts_score(get_parts(u),
                              get_parts(v),
                              distance)
  if (monitor(u) or monitor(v)) and score >= NEUTRAL:
    log("# Score (%s, %s) = %s" % (blurb(u), blurb(v), score))
    #log("#  %s" % (get_parts(u),))
    #log("#  %s" % (get_parts(v),))
  return score

EPITHET_MASK = 32
VICINITY_MASK = 16
YEAR_MASK = 8
TOKEN_MASK = 4
GENUS_MASK = 2
MIDDLE_MASK = 1

# Score potentially contypic names from 0 to 100

def compute_parts_score(p, q, distance=None):
  hits = misses = 0

  if p.epithet != None and q.epithet != None:
    if p.epithet == q.epithet: hits |= EPITHET_MASK
    else: misses |= EPITHET_MASK

  # Proximity within the hierarchy essential if no genus match
  if distance != None:
    # 15 = too far apart to link
    # 9 = neutral, most distant linkable
    # 0 = equated, exact hit
    if distance <= 9: hits |= VICINITY_MASK
    elif distance > 15: misses |= VICINITY_MASK

  if p.year != None and q.year != None:
    if p.year == q.year: hits |= YEAR_MASK
    else: misses |= YEAR_MASK

  if p.token != None and q.token != None:
    if p.token == q.token: hits |= TOKEN_MASK
    else: misses |= TOKEN_MASK

  if p.genus != None and q.genus != None:
    #log("# 1 comparing %s, %s (%s, %s)" % (p.genus, q.genus, p.protonymp, q.protonymp))
    if p.genus == q.genus:
      #log("# 2 comparing %s, %s" % (p.genus, q.genus))
      hits |= GENUS_MASK
    elif p.protonymp == True and q.protonymp == True:
      # If both are protonyms, the genera have to match
      #log("# 4 comparing %s, %s -> miss" % (p.genus, q.genus))
      misses |= GENUS_MASK
    else:
      # Allow motion between genera as long as it's flagged with ()
      #  ... or if protonym status is unknown?
      #log("# 3 comparing %s, %s -> ?" % (p.genus, q.genus))
      pass

  if p.middle != None and q.middle != None:
    pmid = p.epithet if p.middle == '' else p.middle
    qmid = q.epithet if q.middle == '' else q.middle
    if pmid == qmid: hits |= MIDDLE_MASK
    else: misses |= MIDDLE_MASK

  if misses > 0: return NEUTRAL - misses
  else: return NEUTRAL + hits

def explain(score):
  def explode(things):
    return ', '.join((z
                      for (i, z)
                      in zip(range(0, 6),
                             (
                              "middle",   # MIDDLE_MASK
                              "genus",    # GENUS_MASK = 2 << 1
                              "epithet",  # EPITHET_MASK
                              "token",    # TOKEN_MASK
                              "year",     # YEAR_MASK
                              "vicinity", # VICINITY_MASK
                             ))
                      if (things & (1<<i)) != 0))
  if score == NEUTRAL:
    word = "neutral"
    bits = 0
  if score > NEUTRAL:
    word = "homotypic" if homotypic_score(score) else "similar, review"
    bits = score - NEUTRAL
  else:
    word = "heterotypic"
    bits = NEUTRAL - score
  return "%s %s(%s)" % (word, explode(bits))
    
# u and species assumed to be in the same checklist, yes?
# Assumes both are descended from the same species (or genus?).
# Yes, if they're near one another and the epithet stems match.
# 5 is an estimate of typical max distance between a species and 
#  any descendant... hmm...  there has to be a better way

def homotypic(u, v, dist):
  return homotypic_score(compute_score(u, v))

def homotypic_score(score):
  return score > 0 and (score & mask1 == mask1 or
                        score & mask2 == mask2)

mask1 = (EPITHET_MASK | GENUS_MASK | YEAR_MASK)
mask2 = (EPITHET_MASK | VICINITY_MASK)

def score_to_ship(score):
  if homotypic_score(score): return HOMOTYPIC    # EQ
  elif score <  NEUTRAL:     return HETEROTYPIC  # NEQ
  else: return REVIEW                         # NOINFO

# -----------------------------------------------------------------------------

# distance is thresholded, so it only really matters whether it's small
# or large

# For mammals, tip to root is expected to be about 13... 
# For 2M species, tip to root is expected to be about 20... 

def compute_distance(u, v, get_pre_estimate):
  assert separated(u, v)
  up = get_pre_estimate(u, None)
  vp = get_pre_estimate(v, None)
  if up and vp:
    return int((compute_half_distance(u, v, up) +
                compute_half_distance(v, u, vp))/2)
  else:
    return None

def compute_half_distance(u, v, up):
  # u < u1 <= (v1 < m > v)
  assert separated(u, v)
  v1 = up.record
  dist = distance_in_checklist(get_outject(v), get_outject(v1))
  return dist

def distance_in_checklist(x1, x2):
  m = simple.mrca(x1, x2)
  assert m
  return (distance_on_lineage(x1, m) +
          distance_on_lineage(x2, m))

def distance_on_lineage(x, m):
  assert simple.simple_le(x, m)
  if x == m:
    return 0
  return (distance_to_parent(x) +
          distance_on_lineage(get_superior(x).record, m))

def distance_to_parent(u):
  sup = get_superior(u)
  return lg(max(1,len(get_children(sup.record, ()))))

def lg(x):
  return math.log(x)/log2
log2 = math.log(2)

# Convenience.  Phase this out?  Or rename it?

def get_link(u, default=-19):
  uf = really_get_typification_uf(u, None)
  if uf:
    (xid, u2, v) = uf.payload()
    return v if (v and separated(u, v)) else u2
  return None


# -----------------------------------------------------------------------------
# Exemplars (separate module maybe?)

# Apply this to an exemplar id to obtain an exemplar union/find node,
# and return the associated taxon record that's in same checklist as z.

def xid_to_record(AB, xid, z):
  uf = AB.exemplar_ufs[xid]
  (_, u, v) = uf.payload()
  return u if isinA(AB, z) else v

def xid_to_opposite_record(AB, xid, z):
  uf = AB.exemplar_ufs[xid]
  (_, u, v) = uf.payload()
  return v if isinA(AB, z) else u
