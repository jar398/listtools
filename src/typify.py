from parse import PROBE, duplicate_parts

import util
import property as prop
import simple

from util import log, MISSING, UnionFindable
from rcc5 import EQ, NEQ, OVERLAP, DISJOINT
from checklist import get_parts, monitor, \
  get_superior, get_children, get_synonyms, \
  is_accepted, blurb, blorb, get_scientific, get_primary_key, get_rank, \
  get_canonical, get_source_tag, get_nomenclatural_status, \
  get_children, get_synonyms, \
  get_duplicate_from, set_duplicate_from, all_records
from workspace import separated, get_outject, get_workspace, local_sup, get_source
from workspace import isinA, isinB, local_accepted, all_records
from simple import simple_le, distance_in_checklist
from specimen import sid_to_epithet, same_specimens, \
  equate_typifications, equate_specimens, \
  get_typification, get_specimen_id, \
  sid_to_specimen, get_typifies

HOMOTYPIC     = 10
MOTION        = 6
REVIEW        = 5
HETEROTYPIC   = 0

def find_endohomotypics(AB):
  def process_checklist(AB):
    def process_inf(x, p):          # p = parent of x
      process(x)
      comparison = compare_records(AB.in_left(x), AB.in_left(p))
      return classify_comparison(comparison) >= MOTION
    def process(x):
      u = AB.in_left(x)
      candidates = []
      for c in get_children(x, ()):
        if process_inf(c, x):
          # c nominates itself as the type of x
          candidates.append(c)
      n = len(candidates)
      if n == 1:
        # maybe return a score, and sort ???
        equate_typifications(AB.in_left(candidates[0]), u)
      elif n > 1:
        # no candidate, or ambiguity
        log("# Putatively homotypic children of %s: %s, %s" %
            (blurb(x), blurb(candidates[0]), blurb(candidates[1])))
      for c in get_synonyms(x, ()):
        # N.b. synonyms do not have descendants
        if process_inf(c, x):
          # c is a homotypic synonym of x
          equate_typifications(AB.in_left(c), u)
        elif 'homotypic' in get_nomenclatural_status(c, ''):
          equate_typifications(AB.in_left(c), u)
    process(AB.A.top)
  process_checklist(AB)

# TBD
def duplicates(u, v):
  x = get_outject(u)
  y = get_outject(v)
  return (duplicate_parts(get_parts(x), get_parts(y)) and
          (not get_rank(x) or
           not get_rank(y) or
           get_rank(x) == get_rank(y)))

# This can be configured to run once or run twice.  'last' means we're
# on the last pass so if there was a first pass, distances are available.

def find_typifications(AB, subprobs, get_estimate, last):
  # This sets the 'typification_uf' property of ... some ... records.

  n = 1
  for (key, (us, vs)) in subprobs.items():  # For each subproblem
    if n % 1000 == 0 or PROBE in key:
      log("# Subproblem %s %s %s %s %s" %
          (n, len(us), len(vs), blurb(us[0]), blurb(vs[0])))
    n += 1
    u_matches = {}    # u_sid -> (class, [v_spec, ...])
    v_matches = {}    # v_sid -> (class, [u_spec, ...])
    for i in range(0, len(us)):
      u = us[i]
      if monitor(u): log("# Subproblem row: '%s' '%s'" % (key, blorb(u)))
      u_spec = get_typification(u)
      for j in range(0, len(vs)):
        v = vs[j]
        v_spec = get_typification(v)
        comparison = compare_records(u, v) # shows!
        classified = classify_comparison(comparison)
        if classified == MOTION and not near_enough(u, v):
          if monitor(u) or monitor(v):
            log("# Too far apart: %s -> %s" %
                (blorb(u), blorb(v)))
          classified = REVIEW
        if classified >= REVIEW:
          observe_match(u_spec, v_spec, u_matches, classified)
          observe_match(v_spec, u_spec, v_matches, classified)
      # end j loop
    # end i loop
    # u_matches : u_sid -> (v_clas, v_specs)
    for (u_sid, (v_clas, v_specs)) in u_matches.items():
      u_spec = sid_to_specimen(AB, u_sid)
      if v_clas < REVIEW:       # HETEROTYPIC
        pass
      elif v_clas < MOTION:
        log("# Review 1 %s" %         # or, make a note of it for review
            (blorb(get_typifies(u_spec)),))
      elif len(v_specs) > 1:
        log("# B ambiguity %s %s -> %s, %s" %
            (explain_classified(v_clas),
             blorb(get_typifies(u_spec)),
             blorb(get_typifies(v_specs[0])),
             blorb(get_typifies(v_specs[1]))))
      else:
        v_spec = v_specs[0]
        v_sid = get_specimen_id(v_spec)
        # v_matches : v_sid -> (u_clas, u_specs)
        results = v_matches.get(v_sid)
        if results:
          (u_clas, u_specs) = results
          if u_clas < REVIEW:
            pass
          elif u_clas < MOTION:    # REVIEW
            log("# Review 2" %         # or, make a note of it for review
                (blorb(get_typifies(v_spec)),))
          elif len(u_specs) > 1:
            log("# A ambiguity %s %s -> %s, %s" %
                (explain_classified(u_clas),
                 blorb(get_typifies(v_spec)),
                 blorb(get_typifies(u_specs[0])),
                 blorb(get_typifies(u_specs[1]))))
          elif not u_specs[0] is u_spec:
            # How can this happen ???
            log("# Vee %s ->\n  %s -> %s" %
                (blorb(get_typifies(u_spec)),
                 blorb(get_typifies(v_spec)),
                 blorb(get_typifies(u_specs[0])))),
          else:
            if False and (u_clas == MOTION or v_clas == MOTION):
              log("# Genus change %s %s" % (explain_classified(u_clas),
                                            explain_classified(v_clas)))
            equate_specimens(u_spec, v_spec)
        else:
          log("# No return match: %s -> %s" %   # No return match for v - shouldn't happen
              (blorb(get_typifies(u_spec)),
               blorb(get_typifies(v_spec))))
      # End u_sid loop.

      # Nope.  This is not a match:
      #   Tylonycteris malayana Chasen, 1940 -> Euroscaptor malayanus (Chasen, 1940)
      # end subproblem loop

  equate_typifications(AB.in_left(AB.A.top),
                       AB.in_right(AB.B.top))

# u_matches : u_sid -> (v_clas, v_specs)

def observe_match(u_spec, v_spec, u_matches, classified):
  u_sid = get_specimen_id(u_spec)

  have = u_matches.get(u_sid)   # (v_clas, [v_spec, ...])
  if have:
    (v_clas, v_specs) = have
    if classified < v_clas:
      # Not so good as what we have already
      pass
    elif classified > v_clas:
      # Replace specs with new and shiny
      u_matches[u_sid] = (classified, [v_spec])
    else:
      v_spec = v_spec.find()
      if v_spec in v_specs:
        pass
      else:
        # Tie - add to list
        v_specs.append(v_spec)
  else:
    u_matches[u_sid] = (classified, [v_spec])

# Are u0 and v0 near enough that v0 might be a genus change of u0 or
# vice versa?

def near_enough(u, v):
  f = get_family(get_outject(u))
  if not f:
    log("# No family: %s" % blorb(u))
    return False
  g = get_family(get_outject(v))
  if not g:
    log("# No family: %s" % blorb(v))
    return False
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
  return get_canonical(x)

# More important -> lower number, earlier in sequence

def unimportance(u):
  x = get_outject(u)
  parts = get_parts(x)
  if parts.epithet == MISSING: imp = 4      # Foo
  elif parts.middle == parts.epithet: imp = 1     # Foo bar bar
  elif parts.middle == None or parts.middle == '':  imp = 2     # Foo bar
  else: imp = 3                         # Foo bar baz
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
    pep = p.genus.lower() if p.epithet == '' and p.genus else p.epithet
    qep = q.genus.lower() if q.epithet == '' and q.genus else q.epithet
    if pep == qep: hits |= EPITHET_MASK
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
  if classified == HOMOTYPIC:
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
  uf = maybe_get_typification(u, None)
  if uf:
    (_, u2, v) = uf.payload()
    return v if (v and separated(u, v)) else u2
  return None
