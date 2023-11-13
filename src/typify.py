
import math
import util
import property as prop
import simple

from rcc5 import EQ, NEQ
from checklist import get_parts, monitor, get_superior, get_children, is_accepted
from workspace import separated, get_outject

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
  if not m is y1: return v1
  if not m is y2: return v2
  # See whether v2 is an improvement over v1
  # Put this logic in pick_better_record ??
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
    log("# linkage: Ambiguous: %s -> %s & %s" %
        (blurb(u), blurb(v1), blurb(v2)))
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

def homotypic(u, species):
  return homotypic_score(compute_score(u, species))

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
  y  = get_outject(v)
  y1 = get_outject(v1)
  m = simple.mrca(y1, y)
  dist = (distance_on_lineage(y1, m) +
          distance_on_lineage(y, m))
  return dist

def distance_on_lineage(u1, u2):
  if u1 == u2:
    return 0
  return (distance_to_parent(u1) +
          distance_on_lineage(get_superior(u1).record, u2))

def distance_to_parent(u):
  sup = get_superior(u)
  return lg(max(1,len(get_children(sup.record, ()))))

def lg(x):
  return math.log(x)/log2
log2 = math.log(2)

