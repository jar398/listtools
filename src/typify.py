from parse import PROBE

import util
import property as prop
import simple

from util import log, MISSING, UnionFindable
from rcc5 import EQ, NEQ, OVERLAP, DISJOINT
from checklist import get_parts, monitor, get_superior, get_children, \
  is_accepted, blurb, get_scientific, get_primary_key, get_rank, \
  get_canonical, get_source_tag
from workspace import separated, get_outject, get_workspace, local_sup, get_source
from workspace import isinA, isinB

# This runs twice.  'last' means we're on the last pass.

def endo_typifications(AB, subproblems):
  for (key, (us, vs)) in subproblems.items():
    try_in_same_checklist(AB, us, False)
    try_in_same_checklist(AB, vs, False)

def find_typifications(AB, subproblems, get_estimate, last):
  # This sets the 'typification_uf' property of ... some ... records.

  n = 1
  for (key, (us, vs)) in subproblems.items():
    if n % 1000 == 0 or PROBE in key:
      log("# Subproblem %s %s %s %s %s" % (n, len(us), len(vs), blurb(us[0]), blurb(vs[0])))
    n += 1
    all_movers = []
    for i in range(0, len(us)):
      # Sorted by unimportance
      u = us[i]
      # If it's a synonym, see if it matches any accepted in same checklist?
      winners = []
      targets = []
      for j in range(0, len(vs)):
        v = vs[j]
        if get_estimate:
          dist = compute_distance(u, v, get_estimate)
        else:
          dist = None
        comparison = compare_records(u, v, dist)
        if monitor(u) or monitor(v):
          log("# Comparison: %s %s %s %s" %
              (blurb(u), blurb(v), explain(comparison), dist))
        # classify as (HETEROTYPIC, REVIEW, HOMOTYPIC)
        classified = classify_comparison(comparison)
        if classified == HOMOTYPIC:
          winners.append(v)
          break  # Experimental
        elif classified == HETEROTYPIC:
          pass
        elif classified == MOTION:
          if last:
            # **** COMPUTE DISTANCE if 2nd pass ****
            if get_rank(v) == 'species':
              if same_vicinity(u, v):
                targets.append(v)
              else:
                if monitor(u) or monitor(v):
                  log("# Too far apart: %s -> %s" %
                      (blorb(u), blorb(v)))
        else:                   # REVIEW
          pass
          # log("# Review: '%s' '%s'" % (blorb(u), blorb(v)))
      # end j loop

      if len(winners) > 1:
        log("# %s winners for %s" % (len(winners), blorb(u)))
      for v in winners:
        equate_typifications(u, v)
      if len(targets) > 0 and get_rank(u) == 'species':
        all_movers.append((u, targets))
    # end i loop

    if last:
      unique = None
      us_by_id = {u.id : u for u in us}
      vs_by_id = {v.id : v for v in us}
      u_by_vicinity = {}        # maps record id to record
      v_by_vicinity = {}
      # Find matches that are two-way unique
      for (u, targets) in all_movers:
        for v in targets:       # list
          if v.id in u_by_vicinity:
            u_by_vicinity[v.id] = False
          u_by_vicinity[v.id] = u
          if u.id in v_by_vicinity:
            v_by_vicinity[u.id] = False
          v_by_vicinity[u.id] = v
      for (u_id, v) in v_by_vicinity.items():
        if v == False:
          log("# Ambiguous motion %s -> ..." %
              (blorb(u),))
        else:
          u = us_by_id[u_id]
          u_rev = u_by_vicinity[v.id]
          if u_rev == False:
            log("# Ambiguous motion %s, ... -> %s" %
                (blorb(u), blorb(v)))
          else:
            if monitor(u) or monitor(v):
              log("# Unique match %s <-> %s" %
                  (blorb(u), blorb(v)))
            equate_typifications(u, v)


      # Nope.  This is not a match:
      #   Tylonycteris malayana Chasen, 1940 -> Euroscaptor malayanus (Chasen, 1940)
      # end subproblem loop

  equate_typifications(AB.in_left(AB.A.top),
                       AB.in_right(AB.B.top))

def same_vicinity(u0, v0):
  f = get_family(u0)
  if not f:
    log("# no family: %s" % blurb(u0))
    return False
  g = get_family(v0)
  if not g:
    log("# no family: %s" % blurb(v0))
    return False
  return f == g

# Not a match:
#  Maxomys pagensis (Miller, 1903) -> Macaca pagensis (G. S. Miller, 1903)
#  (families Muridae, Cercopithecidae)

def get_family(u):              # returns canonical name
  while get_rank(u) != 'family':
    sup = local_sup(get_source(u), u)
    if sup == None:
      log("# No sup: %s %s" % (get_primary_key(u), blurb(u)))
      return None
    u = sup.record
  return get_canonical(u)

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
  p = get_parts(get_outject(u))
  if p.epithet == MISSING: imp = 4      # Foo
  elif p.middle == MISSING: imp = 2     # Foo bar
  elif p.middle == p.epithet: imp = 1     # Foo bar bar
  else: imp = 3                         # Foo bar baz
  if not is_accepted(get_outject(u)):
    imp += 10
  return (imp, get_primary_key(u))

def try_in_same_checklist(AB, us, last):
  which = 'A' if isinA(AB, us[0]) else 'B'
  for i in range(0, len(us)):
    u1 = us[i]
    x1 = get_outject(u1)
    others = []
    for j in range(i+1, len(us)):
      u2 = us[j]
      x2 = get_outject(u2)
      dist = simple.distance_in_checklist(x1, x2)
      comparison = compare_records(u1, u2, dist)
      classified = classify_comparison(comparison)
      if classified == HOMOTYPIC:
        rel = simple.compare_per_checklist(x1, x2)
        if rel.relationship != DISJOINT:
          equate_typifications(u1, u2)
        else:
          # Synonyms not classified properly?
          others.append(x2)
      elif classified == HETEROTYPIC:
        pass
      elif classified == MOTION:
        # treat as heterotypic within checklist...
        pass                    # ?
      elif monitor(u1) or monitor(u2):
        # Possible missed synonymy ??  Collate by token/year ?
        log("# To review in %s: '%s' '%s'" %
            (which, blorb(x1), blorb(x2)))
    if len(others) > 0:
      for o in others:
        b1 = blorb(x1)
        b2 = blorb(o)
        if b1 == b2:
          log("** Duplicate record in %s: '%s'" % (which, b1))
        else:
          log("** Duplicate record in %s: '%s', '%s'" % (which, b1, b2))
    else:
      pass

# u and v are in workspace but may or may not be from same checklist

def equate_typifications(u, v):     # opposite checklists. u might be species
  if u is not v:
    equate_typification_ufs(get_typification_uf(u), get_typification_uf(v))
    if monitor(u) or monitor(v):
      log("# Equating type '%s' = type '%s'", (u, v))
  return u

def equate_typification_ufs(uf, vf):
  if same_exemplar(uf, vf): return

  (i1, u1, v1) = uf.payload()
  (i2, u2, v2) = vf.payload()

  ef = uf.absorb(vf)          # ef happens to be uf
  r = ef.payload()

  # Not sure this is necessary, but should at least be harmless
  if i2 == None: r[0] = i1
  elif i1 == None: r[0] = i2
  elif i1 != i2: r[0] = min(i1, i2)

  # What if ambiguity, i.e. pick_better_record returns False?
  r[1] = pick_better_record(u1, u2)
  r[2] = pick_better_record(v1, v2)
  return ef

# Only workspace nodes have uf records

def get_typification_uf(u):
  #log("# Thinking about typification for %s" % blurb(u))
  probe = really_get_typification_uf(u, None)
  if probe: return probe
  AB = get_workspace(u)
  r = [None, u, None] if isinA(AB, u) else [None, None, u]
  uf = UnionFindable(r)
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
    (_, u, v) = uf.payload()
    if u and v:
      return uf
  return None
  
def get_exemplar_id(uf):
  r = uf.payload()
  (xid, u, v) = r
  assert u and v
  if xid == None:
    # Create exemplar id (for set operations) on demand
    ws = get_workspace(u)
    xid = len(ws.exemplar_ufs)
    r[0] = xid
    ws.exemplar_ufs[xid] = uf
    #log("# Exemplar %s: (%s) <-> (%s)" % (xid, blurb(u), blurb(v)))
  return xid

def get_typification_record(z):
  uf = really_get_typification_uf(z, None)
  if uf:
    r = uf.payload()
    (_, u, v) = r
    if u and v:
      return r
  return None

def same_exemplar(uf1, uf2):
  return uf1.find() is uf2.find()

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
  parts1 = get_parts(get_outject(v1))
  parts2 = get_parts(get_outject(v2))

  g1 = badness(parts1)
  g2 = badness(parts2)
  if g1 < g2:
    return v2
  elif g1 > g2:
    return v1

  # Fall through - ties all the way
  if get_scientific(v1) < get_scientific(v2):
    return v1
  else:
    return v2
    
# A b b before A c b before A b
# I think that A c b shouldn't happen
def badness(parts):
  if parts.middle == None or parts.middle == '':
    return 3
  if parts.middle == parts.epithet:
    return 1
  else:
    return 2


HOMOTYPIC   = 10
MOTION      = 6
REVIEW      = 5
HETEROTYPIC = 0

# Compare potentially homotypic taxa.
# distance is distance (perhaps estimated) between u and v in combined model.
# If distance (in tree) is close, give a pass on genus mismatch.

def compare_records(u, v, distance=None):
  comparison = compare_parts(get_parts(get_outject(u)),
                             get_parts(get_outject(v)),
                             distance)
  if (monitor(u) or monitor(v)):
    log("# Comparison (%s, %s) = %s" % (blurb(u), blurb(v), comparison))
    #log("#  %s" % (get_parts(get_outject(u)),))
    #log("#  %s" % (get_parts(get_outject(v)),))
  return comparison

EPITHET_MASK = 32
VICINITY_MASK = 16
YEAR_MASK = 8
TOKEN_MASK = 4
GENUS_MASK = 2
MIDDLE_MASK = 1

NEAR_THRESHOLD = 30
FAR_THRESHOLD = 60

# Compare potentially homotypic names.  Returns (m1, m2) where m1 and
# m2 are integer masks, m1 for differences and m2 for similarities.

def compare_parts(p, q, distance=None):
  hits = misses = 0

  if p.epithet != None and q.epithet != None:
    if p.epithet == q.epithet: hits |= EPITHET_MASK
    else: misses |= EPITHET_MASK

  # Proximity within the hierarchy essential if no genus match
  if distance != None:
    if distance <= NEAR_THRESHOLD:
      hits |= VICINITY_MASK
    elif distance > FAR_THRESHOLD:
      misses |= VICINITY_MASK

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
      if False and (PROBE in p.scientific or PROBE in q.scientific):
        log("# 4 comparing %s, %s -> miss" % (p.genus, q.genus))
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

  return (misses, hits)

def classify_comparison(comp):
  (misses, hits) = comp

  # No epithet -> mismatch (or very likely mismatch)
  if (hits & EPITHET_MASK) == 0:
    return HETEROTYPIC          # ??

  # Genus -> match
  if misses == 0 and (hits & GENUS_MASK) != 0:
    return HOMOTYPIC

  if misses != 0:
    if misses == GENUS_MASK:
      return REVIEW          # ??
    if misses == GENUS_MASK | YEAR_MASK:
      return REVIEW          # ??
    if misses == GENUS_MASK | TOKEN_MASK:
      return REVIEW          # ??
    else:
      return HETEROTYPIC          # ??

  return MOTION

def explain(comparison):
  def explode(things):
    return ', '.join((z
                      for (i, z)
                      in zip(range(0, 6),
                             (
                              "middle",   # MIDDLE_MASK = 1 << 1
                              "genus",    # GENUS_MASK = 2 << 1
                              "token",    # TOKEN_MASK
                              "year",     # YEAR_MASK
                              "vicinity", # VICINITY_MASK
                              "epithet",  # EPITHET_MASK
                             ))
                      if (things & (1<<i)) != 0))
  word = "comparison"
  if False:
    classified = classify_comparison(comparison)
    if classified == HOMOTYPIC:
      word = "homotypic"
    elif classified == HETEROTYPIC:
      word = "heterotypic"
    else:
      word = "review"

  (misses, hits) = comparison
  return ("%s. hits(%s) misses(%s)" %
          (word, explode(hits), explode(misses)))
    
# -----------------------------------------------------------------------------

# distance is thresholded, so it only really matters whether it's small
# or large

# For mammals, tip to root is expected to be about 13... 
# For 2M species, tip to root is expected to be about 20... 

def compute_distance(u, v, get_estimate):
  assert separated(u, v)
  half1 = compute_half_distance(u, v, get_estimate)
  half2 = compute_half_distance(v, u, get_estimate)
  return int((half1 + half2)//2)

def compute_half_distance(u, v, get_estimate):
  assert separated(u, v)
  u_est = get_estimate(u, None)
  # u < u1 <= (v1 < m > v)
  dist = simple.distance_in_checklist(get_outject(v), get_outject(u_est.record))
  return dist

# Convenience.  Phase this out?  Or rename it?

def get_link(u, default=-19):
  uf = really_get_typification_uf(u, None)
  if uf:
    (_, u2, v) = uf.payload()
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
  parts = get_parts(get_outject(u))
  return parts.epithet or parts.genus 

def blorb(u):
  if is_accepted(u):
    return get_scientific(u)
  else:
    return get_scientific(u) + "*"
