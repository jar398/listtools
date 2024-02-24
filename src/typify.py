import math
import util
import property as prop
import simple

from util import log, MISSING, UnionFindable
from rcc5 import EQ, NEQ, DISJOINT
from checklist import get_parts, monitor, get_superior, get_children, \
  is_accepted, blurb, get_scientific
from workspace import separated, get_outject, get_workspace
from workspace import isinA, isinB

# This runs twice.  'second' means we're on the second pass.

def find_typifications(AB, subproblems, get_pre_estimate, second):
  # This sets the 'typification_uf' property of ... some ... records.

  n = 1
  for (key, (us, vs)) in subproblems.items():
    if n % 1000 == 0:
      log("# Subproblem %s %s %s %s %s" % (n, len(us), len(vs), blurb(us[0]), blurb(vs[0])))
    n += 1
    try_in_same_checklist(us, second)
    try_in_same_checklist(vs, second)
    success = False
    partial_comparison = (0,0)  # So far
    for i in range(0, len(us)):
      # Sorted by unimportance
      u = us[i]
      # If it's a synonym, see if it matches any accepted in same checklist?
      winners = []
      ambiguous = False
      for j in range(i, len(vs)):
        v = vs[j]
        # classify as A1, A2, A3 (HETEROTYPIC, REVIEW, HOMOTYPIC)
        # **** COMPUTE DISTANCE if 2nd pass ****
        dist = compute_distance(u, v, get_pre_estimate)
        comparison = compare_records(u, v, dist)
        if comparison > partial_comparison: partial_comparison = comparison
        if homotypic_comparison(comparison):
          success = True
          winners.append(v)
      losers = find_ambiguity(winners)
      if losers and not second:
        # silently treat ambiguities as mismatches in first pass.
        # this will fail if e.g. every species has its type species.
        # may show up as typification_compatible(v1, v2)
        log("* Ambiguous pass 1 match %s -> %s, %s" %
            (blurb(u), blurb(losers[0]), blurb(losers[1])))
      else:
        if losers:
          log("* Ambiguous pass 2 match %s -> %s, %s" %
              (blurb(u), blurb(losers[0]), blurb(losers[1])))
        for v in winners:
          equate_typifications(u, v)
    if not success and False:
      log("* No match via %s in %s x %s; best %s" %
          (key, len(us), len(vs), explain(partial_comparison)))

  equate_typifications(AB.in_left(AB.A.top),
                       AB.in_right(AB.B.top))

def find_ambiguity(winners):
  v0 = None
  for v in winners:
    if (v0 != None and
        not typification_compatible(v, v0)):
      return (v0, v)
  return None

# This implements transitivity of equality, I think?
# u and v are in the same checklist

def typification_compatible(u, v):
  uf = really_get_typification_uf(u, None)
  if uf:
    vf = really_get_typification_uf(v, None)
    if vf:
      return uf.payload() is vf.payload()
  return True

# duplicate.  weird.  used? yes, in theory.py
def known_same_typification(u, v):
  uf = really_get_typification_uf(u, None)
  if uf:
    vf = really_get_typification_uf(v, None)
    if vf:
      return uf.payload() is vf.payload()
  return False

def unimportance(u):
  p = get_parts(u)
  if p.epithet == MISSING: imp = 4      # Foo
  if p.middle == p.epithet: imp = 1     # Foo bar bar
  elif p.middle == MISSING: imp = 2     # Foo bar
  else: imp = 3                         # Foo bar baz
  if not is_accepted(get_outject(u)):
    imp += 10
  return imp

def try_in_same_checklist(us, second):
  for i in range(0, len(us)):
    u1 = us[i]
    for j in range(i, len(us)):
      # if i == j: continue  ???
      u2 = us[j]
      x1 = get_outject(u1)
      x2 = get_outject(u2)
      dist = distance_in_checklist(x1, x2) * 2
      comparison = compare_records(u1, u2, dist)
      if homotypic_comparison(comparison):
        rel = simple.compare_per_checklist(x1, x2)
        if rel.relationship != DISJOINT:
          equate_typifications(u1, u2)
        elif second and is_accepted(x1) and is_accepted(x2):
          # Occurs in COL mussels 2022 vs. 2023...
          # looks like it's always Genus sp subsp <-> Genus subsp,
          # and those can always be equated, i.e. they aren't disjoint.
          # Will this confuse downstream analysis ???
          log("## Determining that %s ! %s are homotypic" % (blurb(u1), blurb(u2)))
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
  if i2 == None: i = i1
  elif i1 == None: i = i2
  else: i = min(i1, i2)         # Not sure if this will work
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


# z is in AB.
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
  if v1 == False: return v1   # Propagate ambiguity block ... ???
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
    if heterotypic_comparison(compare_records(v1, v2)):
      log('')
      log("# %s" % (parts1,))
      log("# %s" % (parts2,))
      log("# middles '%s' '%s'" % (parts1.middle,  parts2.middle))
      log("# %s, %s <= %s" % (blurb(y1), blurb(y2), blurb(m)))
      log("# typify: Ambiguous: %s & %s" %
          (blurb(v1), blurb(v2)))
      return False  # Ambiguous
    elif get_scientific(v1) < get_scientific(v2):
      return v1
    else:
      return v2
  return v1       # arbitrary synonym choice; don't want ambiguous


HOMOTYPIC   = EQ           # Hmm, not sure about this
HETEROTYPIC = NEQ
REVIEW      = HOMOTYPIC | HETEROTYPIC

# For combining negative and positive comparisons:
NEUTRAL = 0


# Compare potentially homotypic taxa.
# distance is distance (perhaps estimated) between u and v in combined model.
# If distance (in tree) is close, give a pass on genus mismatch.

def compare_records(u, v, distance=None):
  comparison = compare_parts(get_parts(u),
                             get_parts(v),
                             distance)
  if (monitor(u) or monitor(v)):
    log("# Comparison (%s, %s) = %s" % (blurb(u), blurb(v), comparison))
    #log("#  %s" % (get_parts(u),))
    #log("#  %s" % (get_parts(v),))
  return comparison

EPITHET_MASK = 32
VICINITY_MASK = 16
YEAR_MASK = 8
TOKEN_MASK = 4
GENUS_MASK = 2
MIDDLE_MASK = 1

NEAR_THRESHOLD = 7
FAR_THRESHOLD = 13              # this is the important one

# Compare potentially homotypic names.  Returns (m1, m2) where m1 and
# m2 are integer masks, m1 for differences and m2 for similarities.

def compare_parts(p, q, distance=None):
  hits = misses = 0

  if p.epithet != None and q.epithet != None:
    if p.epithet == q.epithet: hits |= EPITHET_MASK
    else: misses |= EPITHET_MASK

  # Proximity within the hierarchy essential if no genus match
  if distance != None:
    # 15 = too far apart to link
    # 9 = neutral, most distant linkable
    # 0 = equated, exact hit
    if distance <= NEAR_THRESHOLD: hits |= VICINITY_MASK
    elif distance > FAR_THRESHOLD: misses |= VICINITY_MASK

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
    else: pass    # Unmatched middle does not imply mismatch

  return (misses, hits)

def explain(comparison):
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
  if homotypic_comparison(comparison):
    word = "homotypic"
  if heterotypic_comparison(comparison):
    word = "heterotypic"
  else:
    word = "equivocal"

  (misses, hits) = comparison
  return ("%s. hits(%s) misses(%s)" %
          (word, explode(hits), explode(misses)))
    
# u and species assumed to be in the same checklist, yes?
# Assumes both are descended from the same species (or genus?).
# Yes, if they're near one another and the epithet stems match.
# 5 is an estimate of typical max distance between a species and 
#  any descendant... hmm...  there has to be a better way

def homotypic(u, v, dist):
  return homotypic_comparison(compare_records(u, v))

def homotypic_comparison(comparison):
  (misses, hits) = comparison
  return (misses == 0 and
          (hits & mask1 == mask1 or
           hits & mask2 == mask2))

mask1 = (EPITHET_MASK | GENUS_MASK)   # Formerly: YEAR_MASK, now ambiguity check
mask2 = (EPITHET_MASK | VICINITY_MASK)

def heterotypic_comparison(comparison):
  (misses, hits) = comparison
  return misses > 0

def comparison_to_ship(comparison):
  if homotypic_comparison(comparison): return HOMOTYPIC    # EQ
  elif heterotypic_comparison(comparison): return HETEROTYPIC  # NEQ
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

def xid_epithet(AB, xid):
  uf = AB.exemplar_ufs[xid]
  (_, u, v) = uf.payload()
  parts = get_parts(u)
  return parts.epithet or parts.genus 
