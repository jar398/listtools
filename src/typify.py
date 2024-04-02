from parse import PROBE, duplicate_parts

import util
import property as prop
import simple

from util import log, MISSING, UnionFindable
from rcc5 import EQ, NEQ, OVERLAP, DISJOINT
from checklist import get_parts, monitor, \
  get_superior, get_children, get_synonyms, \
  is_accepted, blurb, blorb, get_scientific, get_primary_key, get_rank, \
  get_canonical, get_tag, get_nomenclatural_status, \
  get_children, get_synonyms, get_accepted, \
  get_duplicate_from, set_duplicate_from, all_records
from workspace import separated, get_outject, get_workspace, local_sup, get_source
from workspace import isinA, isinB, local_accepted, all_records
from simple import simple_le, distance_in_checklist
from specimen import sid_to_epithet, same_specimens, \
  equate_typifications, equate_specimens, \
  get_typification, get_specimen_id, \
  sid_to_specimen, get_typifies

HOMOTYPIC     = 10
CORRECTION    = 8
MOVED        = 6
REVIEW        = 5               # Not used any more?
DISTANT       = 1
HETEROTYPIC   = 0

def explain_classified(typy):
  if typy == HOMOTYPIC:
    word = "homotypic"
  elif typy == CORRECTION:
    word = "correction"
  elif typy == MOVED:
    word = "moved"           # kinetypic?
  elif typy == REVIEW:
    word = "review"
  elif typy == DISTANT:
    word = "distant"
  elif typy == HETEROTYPIC:
    word = "heterotypic"
  else:
    word = str(typy)
  return word

# Assess degree of homotypy between u and v.

def ws_homotypy(u, v):
  return homotypy(get_outject(u), get_outject(v))

def homotypy(x, y):
  p = get_parts(x)
  q = get_parts(y)
  comparison = compare_parts(p, q)
  clas = classify_comparison(comparison)
  if clas == MOVED and not near_enough(x, y):
    clas = DISTANT
  return clas

# Match quality {mismatch, correction, match}
# Mobility {mismatch, moved, match}

def classify_comparison(comp):
  (misses, hits) = comp

  # No epithet or genus match -> mismatch
  if (hits & (EPITHET_MASK | GENUS_MASK)) == 0:
    return HETEROTYPIC          # ??
  elif misses == 0:
    return HOMOTYPIC
  elif misses == YEAR_MASK:
    return CORRECTION          # ??
  elif misses == TOKEN_MASK:
    return CORRECTION          # ??
  elif misses == GENUS_MASK:
    return MOVED
  else:
    return HETEROTYPIC          # ??

# Duplicates to skip are indicated by get_duplicate_from

def find_endohomotypics(AB):
  check_for_duplicates(AB)      # Doesn't really belong here
  log("# Designating type specimens for %s checklist" % get_tag(AB.A))

  def process_checklist(AB):
    def process(x):
      u = AB.in_left(x)
      candidates = []
      for c in get_children(x, ()):
        process(c)
        k = AB.in_left(c)
        if homotypy(c, x) >= HOMOTYPIC:
          # c nominates itself as the type of x
          # Don't let any dups engage in homotypy
          if not get_duplicate_from(c, None):
            candidates.append(c)
      n = len(candidates)
      if n == 1:
        k = AB.in_left(candidates[0])
        equate_typifications(k, u)
      elif n > 1:
        # maybe compute a score, and sort ???
        log("** Multiple children of %s are homotypic to it: %s" %
            (blurb(x), "; ".join(map(blurb, candidates))))
      for c in get_synonyms(x, ()):
        # N.b. synonyms do not have descendants
        process(c)              # Unnecessary because no inferiors
        k = AB.in_left(c)
        if (homotypy(c, k) >= MOVED or
            'homotypic' in get_nomenclatural_status(c, '')):
          # c is a homotypic synonym of x
          if not get_duplicate_from(k, None):
            equate_typifications(k, u)
    process(AB.A.top)
  process_checklist(AB)

# Doesn't really belong in this file.  Where does it go?
# Filter by near_enough (family)?

def check_for_duplicates(AB):
  log("* Checking for duplicates in %s checklist" % get_tag(AB.A))
  index_by_name = {}    # (epithet, token, year)
  for x in all_records(AB.A):
    p = get_parts(x)
    key = (p.genus, p.middle, p.epithet, p.token, p.year)
    if key in index_by_name:
      index_by_name[key].append(x)
    else:
      index_by_name[key] = [x]
  for (key, xs) in index_by_name.items():
    accept = []
    synos = []
    xs.sort(key=lambda x: (not is_accepted(x),    # accepted < synonym
                           -len(get_children(x, ())), # more < fewer
                           get_primary_key(x)))       # lower < higher
    accept_c = []   # Accepted, have children: keep or merge
    accept_nc = []  # Accepted, have no chidren: pick one or discard all
    synos = []      # Synonyms: pick one or discard all
    leader = x[0]
    # I wish this could handle corrections.  The following isn't right:
    # xs = (leader,) + tuple(x for x in xs if homotypy(x, leader) >= CORRECTION)
    for x in xs:
      if is_accepted(x):
        if len(get_children(x, ())) > 0:
          accept_c.append(x)
        else:
          accept_nc.append(x)
      else:
        synos.append(x)
    if (len(synos) == 0 and len(accept_nc) == 0 and
        all((simple_le(x, leader) or simple_le(leader, x)
             for x in accept_c))):
      # subgenus, ugh.
      continue

    blurbs = ["** Ambiguous (%s):" % blurb(leader)]
    if len(accept_c) > 0:
      if len(accept_c) > 1:
        action = "keep or merge"
        blurbs.append("%s %s" %
                      (action,
                       "; ".join(map(lambda x:get_primary_key(x), accept_c))))
    if len(accept_nc) > 0:
      if len(accept_c) > 0:
        # flush all the non-child-bearing accepteds
        action = "discard accepteds"   # list the taxonIDs
      elif len(accept_nc) == 1:
        # flush all the synonyms... this will be done below
        action = None
      else:
        action = "choose one accepted dup"
      if action:
        blurbs.append("%s %s" %
                      (action,
                       "; ".join(map(lambda x:get_primary_key(x), accept_nc))))
    if len(synos) > 0:
      if len(accept_c) + len(accept_nc) > 0:
        if len(synos) == 1:
          action = "discard dup with accepted"
        else:
          action = "discard dups with accepteds"
      else:
        action = "choose one accepted among"
      for x in synos: set_duplicate_from(x, leader)
      blurbs.append("%s %s" %
                    (action, "; ".join(map(lambda x:blurb(get_accepted(x)),
                                           synos))))
    log(" ".join(blurbs))

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
      assert not get_duplicate_from(u, None)
      if monitor(u): log("# Subproblem row: '%s' '%s'" % (key, blorb(u)))
      u_spec = get_typification(u)
      for j in range(0, len(vs)):
        v = vs[j]
        assert not get_duplicate_from(v, None)
        v_spec = get_typification(v)
        classified = ws_homotypy(u, v)
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
      elif len(v_specs) > 1:
        log("# B ambiguity (%s): %s -> %s, %s" %
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
          elif len(u_specs) > 1:
            log("# A ambiguity (%s): %s -> %s, %s" %
                (explain_classified(u_clas),
                 blorb(get_typifies(v_spec)),
                 blorb(get_typifies(u_specs[0])),
                 blorb(get_typifies(u_specs[1]))))
          elif not u_specs[0] is u_spec:
            # How can this happen ???
            log("# Wayward (%s, %s): %s ->\n  %s -> %s" %
                (explain_classified(u_clas),
                 explain_classified(v_clas),
                 blorb(get_typifies(u_spec)),
                 blorb(get_typifies(v_spec)),
                 blorb(get_typifies(u_specs[0])))),
          elif u_clas < MOVED:    # REVIEW
            log("# Review: %s -> %s" %         # or, make a note of it for review
                (blorb(get_typifies(u_spec)),
                 blorb(get_typifies(v_spec)),))
          elif u_clas == CORRECTION:
            if False:    # there are way too many of these
             log("# Correction: %s -> %s" %         # or, make a note of it for review?
                 (blorb(get_typifies(u_spec)),
                  blorb(get_typifies(v_spec)),))
            equate_specimens(u_spec, v_spec)
          else:
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

def ws_near_enough(u, v):
  return near_enough(get_outject(u), get_outject(v))

def near_enough(x, y):
  f = get_family(x)
  if not f:
    log("# No family: %s" % blorb(x))
    return False
  g = get_family(y)
  if not g:
    log("# No family: %s" % blorb(y))
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
  return get_canonical(x, None)

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
    if p.genus == q.genus: hits |= GENUS_MASK
    else: misses |= GENUS_MASK

    # elif p.protonymp == True and q.protonymp == True: ...
    # If both are protonyms, the genera have to match

  if p.middle != None and q.middle != None:
    pmid = p.epithet if p.middle == '' else p.middle
    qmid = q.epithet if q.middle == '' else q.middle
    if pmid == qmid: hits |= MIDDLE_MASK
    else: misses |= MIDDLE_MASK

  return (misses, hits)

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
  v1 = get_estimate(u)
  # If m is MRCA of u and v, then u <= v1 <= m >= v
  dist = distance_in_checklist(get_outject(v), get_outject(v1.record))
  return dist

# Convenience.  Phase this out?  Or rename it?

def get_link(u, default=-19):
  uf = maybe_get_typification(u, None)
  if uf:
    (_, u2, v) = uf.payload()
    return v if (v and separated(u, v)) else u2
  return None
