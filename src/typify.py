from parse import PROBE, duplicate_parts

import util
import property as prop
import simple

from util import log, MISSING
from rcc5 import EQ, NEQ, OVERLAP, DISJOINT
from checklist import get_parts, monitor, \
  get_superior, get_children, get_synonyms, get_inferiors, \
  is_accepted, blurb, blorb, get_scientific, get_primary_key, get_rank, \
  get_canonical, get_source_tag, get_nomenclatural_status
from workspace import separated, get_outject, get_workspace, local_sup, get_source
from workspace import isinA, isinB, local_accepted, \
  swap
from simple import simple_le, compare_per_checklist
from specimen import sid_to_epithet, same_specimens, \
  equate_type_ufs, equate_specimens, \
  get_type_uf, get_specimen_id, \
  sid_to_specimen, get_typifies, same_specimens
from proximity import near_enough

HOMOTYPIC     = 10
MOTION        = 8
CORRECTION    = 7

AUTHORLESS    = 5   # no token or year
REVIEW        = 2   # includes CORRECTION + MOTION
HETEROTYPIC   = 0

# Problem perhaps: This creates specimen objects even when they aren't matched.

def find_endohomotypics(AB):
  find_homotypics_in_checklist(AB)
  find_homotypics_in_checklist(swap(AB))

def find_homotypics_in_checklist(AB):
  def process(x):

    epithets = {}
    def consider(c):
      c_ep = get_parts(c).epithet
      if c_ep:
        if c_ep in epithets:
          equate_type_ufs(AB.in_left(c), epithets[c_ep])
        else:
          epithets[c_ep] = AB.in_left(c)

    # Match siblings
    for c in get_inferiors(x):
      process(c)
      classified = relate_records(AB.in_left(c), AB.in_left(x))
      if classified >= MOTION:
        consider(c)

    # Match parent/child
    consider(x)

  process(AB.A.top)

# This finds homotypic siblings.  See exemplar.py for use.
# Compare redundant.py.
# Equations have to happen in a workspace, even when they don't
# span checklists.
# !!! TO DO:  Find both kinds of equations (by epithet):
#   (a) inferior/superior
#   (b) siblings

def find_endohomotypics_original(AB):
  def process_checklist(AB):
    def process(x):
      u = AB.in_left(x)
      # Find parent/child homotypies
      # TO DO: find sibling homotypies
      candidates = []   # those matching x ...
      for c in get_children(x, ()):
        if process_inferior(c, x):
          # c nominates itself as the type of x
          candidates.append(c)
      n = len(candidates)
      if n == 1:
        # maybe return a score, and sort ???
        equate_type_ufs(AB.in_left(candidates[0]), u)
      elif n > 1:
        # no candidate, or ambiguity
        log("# Putatively homotypic children of %s: %s, %s" %
            (blurb(x), blurb(candidates[0]), blurb(candidates[1])))
      for c in get_synonyms(x, ()):
        # N.b. synonyms do not have descendants
        if process_inferior(c, x):
          # c is a homotypic synonym of x
          equate_type_ufs(AB.in_left(c), u)
        elif 'homotypic' in get_nomenclatural_status(c, ''):
          # synonym inferior of accepted
          equate_type_ufs(AB.in_left(c), u)
    def process_inferior(x, p):          # p = parent of x
      process(x)
      # uhh I don't think this is right
      classified = relate_records(AB.in_left(x), AB.in_left(p)) # Type subspecies?
      return classified >= MOTION
    process(AB.A.top)
  process_checklist(AB)
  process_checklist(swap(AB))

# This can be configured to run once or run twice.  'last' means we're
# on the last pass so if there was a first pass, proximity is available.

# Unify the types of A with the types of B.

def find_type_ufs(AB, subprobs, get_estimate, last):
  # This sets the 'type_uf' property of ... some ... records.

  n = 1
  for (key, (us, vs)) in subprobs.items():  # For each subproblem
    if n % 1000 == 0 or PROBE in key:
      log("# Subproblem %s %s %s %s %s" %
          (n, len(us), len(vs), blurb(us[0]), blurb(vs[0])))
    # us and vs are sorted by 'unimportance' (i.e. most important first)
    n += 1
    if any(map(monitor, us)):
      log("# Doing subproblem %s" % key)
      log("#  us = %s" % list(map(blurb, us),))
      log("#  vs = %s" % list(map(blurb, vs),))

    u_matches = {}    # u_sid -> (class, [v_spec, ...])
    v_matches = {}    # v_sid -> (class, [u_spec, ...])
                      # spec = (sid, u, v)
                      # Didn't count on multiple specs per sid!
    for i in range(0, len(us)):
      u = us[i]
      if monitor(u): log("# Subproblem row: '%s' '%s'" % (key, blorb(u)))
      for j in range(0, len(vs)):
        v = vs[j]
        classified = relate_records(u, v) # shows!
        # relate_records checks proximity
        observe_match(u, v, u_matches, classified)
        observe_match(v, u, v_matches, classified)
        if monitor(u):
          log("# observe %s => %s = %s" %
              (blurb(u), blurb(v), explain_classified(classified)))
      # end j loop
    # end i loop
    # u_specs and v_specs are sorted by 'unimportance'
    # u_matches : u_sid -> (v_clas, v_specs)
    # clas means a classification
    for (u_sid, (v_clas, v_specs)) in u_matches.items():
      u_spec = sid_to_specimen(AB, u_sid)
      # if monitor(get_typifies(u_spec)): ...
      if v_clas == MOTION:
        # Maybe log AUTHORLESS and/or REVIEW
        # Put this in the report somehow ??
        u0 = get_typifies(u_spec)
        v0 = get_typifies(v_specs[0])
        log("# %s: %s -> %s" %  # or, make a note of it for review
            (explain_classified(v_clas), blorb(u0), blorb(v0)))
      if v_clas < MOTION:
        pass
      elif len(v_specs) > 1:
        # Problem here
        v0 = get_typifies(v_specs[0]) # record that has given specimen
        v1 = get_typifies(v_specs[1])
        if is_accepted(v0) and is_accepted(v1):
          if compare_per_checklist(get_outject(v0), get_outject(v1)) \
             .relationship & DISJOINT != 0:
            # Another kind of REVIEW
            log("# B ambiguity: %s %s -> %s %s, %s %s" %
                (explain_classified(v_clas),
                 blorb(get_typifies(u_spec)),
                 get_primary_key(v0),
                 blorb(v0),
                 get_primary_key(v1),
                 blorb(v1)))
          else:
            equate_specimens(u_spec, v_specs[0])
            equate_specimens(u_spec, v_specs[1])
      else:
        v_spec = v_specs[0]
        v_sid = get_specimen_id(v_spec)
        # v_matches : v_sid -> (u_clas, u_specs)
        results = v_matches.get(v_sid)
        if results:
          (u_clas, u_specs) = results
          if u_clas == MOTION:
            # also maybe AUTHORLESS, REVIEW ?
            log("# Review 2 %s -> %s" %         # or, make a note of it for review
                (blorb(get_typifies(u_spec)),
                 blorb(get_typifies(v_spec)),))
          if u_clas < MOTION:    # e.g. REVIEW
            pass
          elif len(u_specs) > 1:
            # Select type A b b from u_specs if possible...
            u0 = get_typifies(u_specs[0])
            u1 = get_typifies(u_specs[1])
            if compare_per_checklist(get_outject(u0), get_outject(u1)) \
               .relationship & DISJOINT != 0:
              log("# A ambiguity %s: %s, %s -> %s" %
                  (explain_classified(u_clas),
                   blorb(u0),
                   blorb(u1),
                   blorb(get_typifies(v_spec))))
            else:
              # Can be simultaneously compatible with both.
              # (TBD: deal with 3-way or more ambiguity.)
              equate_specimens(u_specs[0], v_spec)
              equate_specimens(u_specs[1], v_spec)
          elif not same_specimens(u_specs[0], u_spec):
            # How can this happen ???
            log("# Vee %s ->\n  %s -> %s" %
                (blorb(get_typifies(u_spec)),
                 blorb(get_typifies(v_spec)),
                 blorb(get_typifies(u_specs[0])))),
          else:
            # WRONG? - only do this for 'tipward' nodes?
            equate_specimens(u_spec, v_spec)
        else:
          log("# No return match: %s -> %s" %   # No return match for v - shouldn't happen
              (blorb(get_typifies(u_spec)),
               blorb(get_typifies(v_spec))))
      # End u_sid loop.
      # end subproblem loop

  equate_type_ufs(AB.in_left(AB.A.top),
                       AB.in_right(AB.B.top))

# u_matches : u_sid -> (v_clas, v_specs)
#   'spec' is short for 'specimen'

# Ideally, only one match per specimen id.
# 'Classified' might be a failure e.g. HETEROTYPIC

def observe_match(u, v, u_matches, classified):
  u_spec = get_type_uf(u)
  v_spec = get_type_uf(v)

  u_sid = get_specimen_id(u_spec)
  v_sid = get_specimen_id(v_spec)

  have = u_matches.get(u_sid)   # (v_clas, [v_speq, ...])
  if have:
    (v_clas, v_speqs) = have
    if classified > v_clas:
      # Discard previously seen lower value matches.
      if monitor(u):
        log("* Discarding previously seen lower value matches: %s" %
            blurb(u))
      u_matches[u_sid] = (classified, [v_spec])
    elif classified != v_clas:
      if monitor(u):
        log("* Discarding lower value match: %s" %
            blurb(u))
      pass
    elif any(map(lambda v_speq: same_specimens(v_spec, v_speq),
                 v_speqs)):
      # Don't add same specimen twice
      if monitor(u):
        log("* Discarding redundant match: %s" %
            blurb(u))
      pass
    else:
      if monitor(u):
        log("* Adding match: %s" % blurb(u))
      v_speqs.append(v_spec)    # ????
  else:
    u_matches[u_sid] = (classified, [v_spec])

# Compare potentially homotypic taxa in same workspace.

def relate_records(u, v):
  classified = compare_parts(u, v)

  if classified == MOTION:
    if not near_enough(u, v):
    # Prouimity within the hierarchy is essential if no genus match
      return REVIEW

  if monitor(u) or monitor(v):
    log("# Compare %s, %s = %s" %
        (blurb(u), blurb(v), explain_classified(classified)))
  return classified

EPITHET_MASK = 32
VICINITY_MASK = 16
YEAR_MASK = 8
TOKEN_MASK = 4
GENUS_MASK = 2
MIDDLE_MASK = 1

NEAR_THRESHOLD = 30
FAR_THRESHOLD = 60

# Compare potentially homotypic names by examining their syntactic parts.
# Compare Macropus robustus, Amblysomus robustus = heterotypic NOT
# Could be in same or different checklists?  Fix this.

def compare_parts(u, v):
  (misses, hits) = parts_comparison_detail(get_parts(u), get_parts(v))
  return classify_comparison_details(misses, hits)

# Given a comparison of parts, classify it as appropriate

def classify_comparison_details(misses, hits):
  # Epithets must match (perhaps '' ?)
  if ((hits & EPITHET_MASK) == 0):
    return HETEROTYPIC          # Not a hit

  mask = YEAR_MASK | TOKEN_MASK
  if misses & mask == mask:
    return HETEROTYPIC

  elif hits & mask == mask:
    if hits & GENUS_MASK != 0:
      return HOMOTYPIC
    else:
      return MOTION

  elif hits & mask != 0:
    # Single-field error?
    if hits & GENUS_MASK != 0:
      return CORRECTION
    else:
      return REVIEW  # !!????

  elif misses & mask != 0:
    return AUTHORLESS

  return REVIEW

# Compare potentially homotypic names.  Returns (m1, m2) where m1 and
# m2 are integer masks, m1 for differences and m2 for similarities.

def parts_comparison_detail(p, q):
  hits = misses = 0

  if p.epithet != None and q.epithet != None:
    # Treat 'Foo' like 'Foo foo'
    pep = p.genus.lower() if p.epithet == '' and p.genus else p.epithet
    qep = q.genus.lower() if q.epithet == '' and q.genus else q.epithet
    if pep == qep: hits |= EPITHET_MASK
    else: misses |= EPITHET_MASK

  if p.token != None and q.token != None:
    if p.token == q.token: hits |= TOKEN_MASK
    else: misses |= TOKEN_MASK

  if p.year != None and q.year != None:
    if p.year == q.year: hits |= YEAR_MASK
    else: misses |= YEAR_MASK

  if p.genus != None and q.genus != None:
    if p.genus == q.genus: hits |= GENUS_MASK
    else: misses |= GENUS_MASK

  # Something about if both are protonyms then miss ?

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
    
def explain_classified(classified):
  if classified == HOMOTYPIC:
    word = "homotypic"
  elif classified == ERROR:
    word = "error"
  elif classified == MOTION:
    word = "motion"           # kinetypic?
  elif classified == AUTHORLESS:
    word = "authorless"
  elif classified == REVIEW:
    word = "review"
  elif classified == HETEROTYPIC:
    word = "heterotypic"
  else:
    word = str(classified)
  return word

# Convenience.  Phase this out?  Or rename it?

def get_link(u, default=-19):
  uf = maybe_get_type_uf(u, None)
  if uf:
    (_, u2, v) = uf.payload()
    return v if (v and separated(u, v)) else u2
  return None
