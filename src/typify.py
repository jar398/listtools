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
from workspace import isinA, isinB, local_accepted
from simple import simple_le, distance_in_checklist

def endo_typifications(AB, subproblems):
  find_extra_homotypics(AB)
  log("* Matching within the A checklist:")
  for (key, (us, vs)) in subproblems.items():
    find_subproblem_endohomotypics(us)
  log("* Matching within the B checklist:")
  for (key, (us, vs)) in subproblems.items():
    find_subproblem_endohomotypics(vs)

def find_extra_homotypics(AB):
  log("* Scanning for declared homotypic synonyms:")
  def scan(AB):
    for x in all_records(AB.A):
      if 'homotypic' in get_nomenclatural_status(x, ''):
        # GBIF has 400,000 of these
        p = AB.in_left(get_accepted(x))
        u = AB.in_left(x)
        if get_parts(p).epithet != get_parts(q).epithet:
          log("# Homotypic synonym %s ~ %s" % (blurb(u), blurb(p)))
        equate_typifications(u, p)
  scan(AB)
  scan(AB.swap())

def find_subproblem_endohomotypics(us):
  for i in range(0, len(us)):
    u1 = us[i]
    others = []
    for j in range(i+1, len(us)):
      u2 = us[j]
      result = endohomotypic(u1, u2)     # caches result
      if False and monitor(u1):
        log("# Endohomotypic: %s %s => %s" % (blorb(u1), blorb(u2), result))

# Given two records in a workspace coming from the same checklist,
# are their concepts homotypic?

def endohomotypic(u, v):
  AB = get_workspace(u)
  a = local_accepted(AB, u)
  b = local_accepted(AB, v)

  def hom(u, v, thresh):
    if same_typification(u, v): return True      # peephole optimization
    comparison = compare_records(u, v)
    if classify_comparison(comparison) < thresh:
      return False
    x = get_outject(u); y = get_outject(v)
    if simple_le(x, y) or simple_le(y, x):
      equate_typifications(u, v)
      return True
    return False

  result = (hom(u, a, MOTION) and
            hom(a, b, HOMOTYPIC) and
            hom(b, v, MOTION))
  return result

# This can be configured to run once or run twice.  'last' means we're
# on the last pass.

def find_typifications(AB, subproblems, get_estimate, last):
  # This sets the 'typification_uf' property of ... some ... records.

  n = 1
  for (key, (us, vs)) in subproblems.items():
    if n % 1000 == 0 or PROBE in key:
      log("# Subproblem %s %s %s %s %s" % (n, len(us), len(vs), blurb(us[0]), blurb(vs[0])))
    n += 1
    all_movers = []
    if (any(monitor(u) for u in us) or
        any(monitor(v) for v in vs)):
      log("# Subproblem: '%s'" % key)
    for i in range(0, len(us)):
      u = us[i]
      if monitor(u):
        log("# Subproblem row: '%s' '%s'" % (key, blorb(u)))
      # If it's a synonym, see if it matches any accepted in same checklist?
      winner = None
      targets = []
      for j in range(0, len(vs)):
        v = vs[j]
        if get_estimate:
          dist = compute_distance(u, v, get_estimate)
        else:
          dist = None
        comparison = compare_records(u, v, dist) # shows!
        classified = classify_comparison(comparison)
        if classified == HOMOTYPIC:
          if True:  #same_vicinity(u, v):
            if winner:
              if not same_typification(winner, v):
                log("# Induced homotypism: %s\n#   %s ~ %s" % (blorb(u), blorb(v), blorb(winner)))
            equate_typifications(u, v)
            winner = v
            # break
        elif classified == HETEROTYPIC:
          if False and (monitor(u) or monitor(v)):
            log("# Heterotypic: '%s' '%s'" % (blorb(u), blorb(v)))
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
          if monitor(u) or monitor(v):
            log("* Review: '%s' '%s'" % (blorb(u), blorb(v)))
      # end j loop

      if winner:
        pass
      elif len(targets) > 0 and get_rank(u) == 'species':
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
            u_by_vicinity[v.id] = True
          u_by_vicinity[v.id] = u
          if u.id in v_by_vicinity:
            v_by_vicinity[u.id] = True
          v_by_vicinity[u.id] = v
      for (u_id, v) in v_by_vicinity.items():
        if v == True:
          log("# Ambiguous motion %s -> ..." %
              (blorb(u),))
        else:
          u = us_by_id[u_id]
          u_rev = u_by_vicinity[v.id]
          if u_rev == True:
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
    log("# No family: %s" % blorb(u0))
    return False
  g = get_family(v0)
  if not g:
    log("# No family: %s" % blorb(v0))
    return False
  return f == g

# Not a match:
#  Maxomys pagensis (Miller, 1903) -> Macaca pagensis (G. S. Miller, 1903)
#  (families Muridae, Cercopithecidae)

def get_family(u):              # returns canonical name
  while get_rank(u) != 'family':
    sup = local_sup(get_source(u), u)
    if sup == None:
      log("# No sup: %s %s" % (get_primary_key(u), blorb(u)))
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

# u and v are in workspace but may or may not be from same checklist

def equate_typifications(u, v):     # opposite checklists. u might be species
  if u is not v:
    equate_typification_ufs(get_typification_uf(u), get_typification_uf(v))
    if monitor(u) or monitor(v):
      log("# Unifying exemplar(%s) = exemplar(%s)" % (blorb(u), blorb(v)))
  return u

def equate_typification_ufs(uf, vf):
  if same_typification_ufs(uf, vf): return

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

def same_typification_ufs(uf1, uf2):
  return uf1.find() is uf2.find()

def same_typification(u, v):
  return same_typification_ufs(get_typification_uf(u),
                               get_typification_uf(v))

# Only workspace nodes have uf records

def get_typification_uf(u):
  #log("# Thinking about typification for %s" % blorb(u))
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
  
def get_exemplar_record(z):
  uf = really_get_typification_uf(z, None)
  if uf:
    r = uf.payload()
    (_, u, v) = r
    if u and v:
      return r
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
    #log("# Exemplar %s: (%s) <-> (%s)" % (xid, blorb(u), blorb(v)))
  return xid

# This chooses the 'best' record to represent a given exemplar

def pick_better_record(v1, v2):
  # See whether v2 is an improvement over v1 (lower unimportance value)
  if v1 == None: return v2
  if v2 == None: return v1
  g1 = unimportance(v1)
  g2 = unimportance(v2)
  if g1 < g2:                   # 1 is less important than 2
    return v2
  elif g1 > g2:
    return v1
  else:
    assert False
    
# More important -> lower number, earlier in sequence

def unimportance(u):
  parts = get_parts(u)
  if parts.epithet == MISSING: imp = 4      # Foo
  elif parts.middle == parts.epithet: imp = 1     # Foo bar bar
  elif parts.middle == None or parts.middle == '':  imp = 2     # Foo bar
  else: imp = 3                         # Foo bar baz
  x = get_outject(u)
  return (1 if is_accepted(x) else 2,
          imp,
          get_scientific(x, None),
          get_primary_key(x))


HOMOTYPIC   = 10
MOTION      = 6
REVIEW      = 5
HETEROTYPIC = 0

# Compare potentially homotypic taxa.
# distance is distance (perhaps estimated) between u and v in combined model.
# If distance (in tree) is close, give a pass on genus mismatch.

def compare_records(u, v, distance=None):
  return compare_parts(get_parts(get_outject(u)),
                       get_parts(get_outject(v)),
                       distance)

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
  dist = distance_in_checklist(get_outject(v), get_outject(u_est.record))
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
