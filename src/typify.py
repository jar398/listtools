from parse import PROBE, duplicate_parts

import util
import property as prop
import simple

from util import log, MISSING, UnionFindable
from rcc5 import EQ, NEQ, OVERLAP, DISJOINT
from checklist import get_parts, monitor, get_superior, get_children, \
  is_accepted, blurb, get_scientific, get_primary_key, get_rank, \
  get_canonical, get_source_tag, get_nomenclatural_status, \
  get_children, get_synonyms
from checklist import get_duplicate_from, set_duplicate_from
from workspace import separated, get_outject, get_workspace, local_sup, get_source
from workspace import isinA, isinB, local_accepted, all_records
from simple import simple_le, distance_in_checklist

ENDOHOMOTYPIC = 10
HOMOTYPIC     = 9
MOTION        = 6
REVIEW        = 5
HETEROTYPIC   = 0

def endo_typifications(AB, subprobs):
  log("* Matching within the A checklist:")
  find_endohomotypics(AB, subprobs, lambda xy: xy[0])
  log("* Matching within the B checklist:")
  find_endohomotypics(AB, subprobs, lambda xy: xy[1])
  find_extra_homotypics(AB)

def find_extra_homotypics(AB):
  log("* Scanning for declared homotypic synonyms:")
  def scan(AB):
    for x in all_records(AB.A):
      if 'homotypic' in get_nomenclatural_status(x, ''):
        # GBIF has 400,000 of these
        p = AB.in_left(get_accepted(x))
        u = AB.in_left(x)
        if get_parts(p).epithet != get_parts(u).epithet:
          log("# Homotypic synonym %s ~ %s" % (blurb(u), blurb(p)))
        equate_typifications(u, p)
  scan(AB)
  scan(AB.swap())

def find_endohomotypics(AB, subprobs, getit):
  dups = []

  for both in subprobs.values():
    # both = (us, vs)
    find_subproblem_endohomotypics(AB, getit(both), dups)

def find_subproblem_endohomotypics(AB, us, dups):
  dups = []
  for i in range(0, len(us)):
    u1 = us[i]
    x1 = get_outject(u1)
    for j in range(i+1, len(us)):
      u2 = us[j]
      x2 = get_outject(u2)

      classified = endohomotypic(u1, u2)     # caches classified
      if monitor(u1) or monitor(u2):
        log("# %s: '%s' ~ '%s'" % (explain_classified(classified), blorb(u1), blorb(u2)))
      if classified >= ENDOHOMOTYPIC:
        assert same_typification(u1, u2)
      if duplicates(u1, u2):    # Is u1 a duplicate taxon of u2?
        if is_accepted(x1) and not is_accepted(x2):
          log("# Keeping accepted %s, flushing synonym" % blurb(x1))
          pass
        elif not get_duplicate_from(x1, None):
          set_duplicate_from(x1, x1)
        set_duplicate_from(x2, x1)
        dups.append((u1, u2))

  dups.sort(key=lambda u1u2: (get_primary_key(get_outject(u1u2[0])),
                              get_primary_key(get_outject(u1u2[1]))))
  for (u1, u2) in dups:
    x1 = get_outject(u1)
    x2 = get_outject(u2)
    n = inferior_count(x2)
    note = " (%s children)" % n if n > 0 else ""
    log("** %s: Duplicate of %s: %s%s" %
        (get_primary_key(x2), get_primary_key(x1), blorb(u2), note))

def duplicates(u, v):
  return (duplicate_parts(get_parts(u), get_parts(v)) and
          (not get_rank(u) or
           not get_rank(v) or
           get_rank(u) == get_rank(v)))

# Given two records in a workspace coming from the same checklist,
# are their concepts homotypic?

def endohomotypic(u1, u2):
  AB = get_workspace(u1)
  a = local_accepted(AB, u1)
  b = local_accepted(AB, u2)

  def hom(u1, u2, thresh):
    if same_typification(u1, u2): return ENDOHOMOTYPIC      # peephole optimization
    comparison = compare_records(u1, u2)
    classified = classify_comparison(comparison)
    if classified >= thresh:
      x = get_outject(u1); y = get_outject(u2)
      if simple_le(x, y) or simple_le(y, x):
        # Might be duplicate records e.g. MDD1.10 Acomys johannis
        equate_typifications(u1, u2)  # Cache it
        classified = ENDOHOMOTYPIC
    return classified

  answer = min(hom(u1, a, MOTION),
               hom(a, b, HOMOTYPIC), 
               hom(b, u2, MOTION))
  return answer

def inferior_count(u):
  return len(get_children(u, ())) + len(get_synonyms(u, ()))

# This can be configured to run once or run twice.  'last' means we're
# on the last pass.

def find_typifications(AB, subprobs, get_estimate, last):
  # This sets the 'typification_uf' property of ... some ... records.

  n = 1
  for (key, (us, vs)) in subprobs.items():
    if n % 1000 == 0 or PROBE in key:
      log("# Subproblem %s %s %s %s %s" % (n, len(us), len(vs), blurb(us[0]), blurb(vs[0])))
    n += 1
    all_movers = []
    if (any(monitor(u) for u in us) or
        any(monitor(v) for v in vs)):
      log("# Subproblem: '%s'" % key)
    for i in range(0, len(us)):
      u = us[i]
      x = get_outject(u)
      if monitor(u):
        log("# Subproblem row: '%s' '%s'" % (key, blorb(u)))
      # If it's a synonym, see if it matches any accepted in same checklist?
      winner = None
      targets = []
      for j in range(0, len(vs)):
        v = vs[j]
        y = get_outject(v)
        if get_estimate:
          dist = compute_distance(u, v, get_estimate)
        else:
          dist = None
        if same_typification(u, v): break
        comparison = compare_records(u, v, dist) # shows!
        classified = classify_comparison(comparison)
        if classified >= HOMOTYPIC:
          if not winner:
            equate_typifications(u, v)
            winner = v
          elif same_typification(winner, v):
            pass
          else:
            x3 = simple.mrca(get_outject(winner), get_outject(v))
            y = get_outject(v)
            if (local_accepted(AB, v) is local_accepted(AB, winner)
                and get_children(y, ()) == ()):
              # Siblings without children can safely be equated??
              # log("# Homotypism induced by %s, allowed, mrca = %s:\n#   %s ~ %s" % (blorb(u), blorb(x3), blorb(v), blorb(winner)))
              equate_typifications(winner, v)
            else:
              pass
              #log("# Ambiguous matches for %s, mrca = %s:\n#   %s ~ %s" % (blorb(u), blorb(x3), blorb(v), blorb(winner)))
        elif classified == MOTION:
          if last:
            if get_rank(u) == 'species' and get_rank(v) == 'species':
              if same_vicinity(u, v):
                targets.append(v)
              else:
                if monitor(u) or monitor(v):
                  log("# Too far apart: %s -> %s" %
                      (blorb(u), blorb(v)))
        elif classified == HETEROTYPIC:
          pass
        else:                   # REVIEW
          if monitor(u) or monitor(v):
            log("* Review: '%s' '%s'" % (blorb(u), blorb(v)))
      # end j loop

      if winner:
        pass
      elif len(targets) > 0:
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
          # Prefer to match the duplicate that has children
          # (or more children)
          -len(get_children(x, ())),
          -len(get_synonyms(x, ())),
          get_scientific(x, None),
          get_primary_key(x))


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

  (misses, hits) = comparison
  return ("hits(%s) misses(%s)" %
          (explode(hits), explode(misses)))
    
def explain_classified(classified):
  if classified == ENDOHOMOTYPIC:
    word = "endohomotypic"
  elif classified == HOMOTYPIC:
    word = "homotypic"
  elif classified == MOTION:
    word = "motion"           # kinetypic?
  elif classified == REVIEW:
    word = "review"
  elif classified == HETEROTYPIC:
    word = "heterotypic"
  else:
    word = str(classified)
  return word

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
  if u == None:
    return "None"
  x = get_outject(u)
  if get_workspace(u):
    prefix = get_source_tag(x) + "."
  else:
    prefix = get_source_tag(x) + "."
  if is_accepted(x):
    name = get_scientific(x)
  else:
    name = get_scientific(x) + "*"
  return prefix + name
