from parse import PROBE, duplicate_parts

import util
import property as prop
import simple

from util import log, MISSING
from rcc5 import EQ, NEQ, OVERLAP, DISJOINT
from checklist import get_parts, monitor, \
  get_superior, get_children, get_synonyms, get_inferiors, \
  is_accepted, blurb, blorb, get_scientific, get_primary_key, get_rank, \
  get_canonical, get_source_tag, get_nomenclatural_status, \
  all_records
from workspace import separated, get_outject, get_workspace, local_sup, get_source
from workspace import isinA, isinB, local_accepted, all_records, \
  swap
from simple import simple_le
from specimen import sid_to_epithet, same_specimens, \
  equate_typifications, equate_specimens, \
  get_typification, get_specimen_id, \
  sid_to_specimen, get_typifies
from proximity import near_enough

HOMOTYPIC     = 10
MOTION        = 6
REVIEW        = 5
HETEROTYPIC   = 1
HETEROTYPIC2  = 0

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
          equate_typifications(AB.in_left(c), epithets[c_ep])
        else:
          epithets[c_ep] = AB.in_left(c)

    # Match siblings
    for c in get_inferiors(x):
      process(c)
      classified = relate_records(c, x)
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
        equate_typifications(AB.in_left(candidates[0]), u)
      elif n > 1:
        # no candidate, or ambiguity
        log("# Putatively homotypic children of %s: %s, %s" %
            (blurb(x), blurb(candidates[0]), blurb(candidates[1])))
      for c in get_synonyms(x, ()):
        # N.b. synonyms do not have descendants
        if process_inferior(c, x):
          # c is a homotypic synonym of x
          equate_typifications(AB.in_left(c), u)
        elif 'homotypic' in get_nomenclatural_status(c, ''):
          # synonym inferior of accepted
          equate_typifications(AB.in_left(c), u)
    def process_inferior(x, p):          # p = parent of x
      process(x)
      # uhh I don't think this is right
      classified = relate_records(x, p) # Type subspecies?
      return classified >= MOTION
      # was comparison = ws_compare_records(AB.in_left(x), AB.in_left(p))
      # return classify_comparison(comparison) >= MOTION
    process(AB.A.top)
  process_checklist(AB)
  process_checklist(swap(AB))

# This can be configured to run once or run twice.  'last' means we're
# on the last pass so if there was a first pass, proximity is available.

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
                      # spec = (sid, u, v)
                      # Didn't count on multiple specs per sid!
    for i in range(0, len(us)):
      u = us[i]
      if monitor(u): log("# Subproblem row: '%s' '%s'" % (key, blorb(u)))
      u_spec = get_typification(u)
      for j in range(0, len(vs)):
        v = vs[j]
        classified = ws_relate_records(u, v) # shows!
        v_spec = get_typification(v)
        observe_match(u_spec, v_spec, u_matches, classified)
        observe_match(v_spec, u_spec, v_matches, classified)
        if monitor(u):
          log("# observe %s => %s = %s" %
              (blurb(u), blurb(v), explain_classified(classified)))
      # end j loop
    # end i loop
    # u_matches : u_sid -> (v_clas, v_specs)
    # clas means a classification
    for (u_sid, (v_clas, v_specs)) in u_matches.items():
      u_spec = sid_to_specimen(AB, u_sid)
      # if monitor(get_typifies(u_spec)): ...
      if v_clas < REVIEW:       # HETEROTYPIC
        pass
      elif v_clas < MOTION:
        log("# Review 1 %s" %         # or, make a note of it for review
            (blorb(get_typifies(u_spec)),))
      elif len(v_specs) > 1:
        # Problem here
        v0 = get_typifies(v_specs[0])
        v1 = get_typifies(v_specs[1])
        log("# B ambiguity %s %s -> %s %s, %s %s" %
            (explain_classified(v_clas),
             blorb(get_typifies(u_spec)),
             get_primary_key(v0),
             blorb(v0),
             get_primary_key(v1),
             blorb(v1)))
        # Two records with the same type specimen ?
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
            log("# Review 2 %s" %         # or, make a note of it for review
                (blorb(get_typifies(v_spec)),))
          elif len(u_specs) > 1:
            log("# A ambiguity %s %s -> %s, %s" %
                (explain_classified(u_clas),
                 blorb(get_typifies(v_spec)),
                 blorb(get_typifies(u_specs[0])),
                 blorb(get_typifies(u_specs[1]))))
          elif get_specimen_id(u_specs[0]) != u_sid:
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
#   'spec' is short for 'specimen'

# Ideally, only one match per specimen id.
# 'Classified' might be a failure e.g. HETEROTYPIC

def observe_match(u_spec, v_spec, u_matches, classified):
  u_sid = get_specimen_id(u_spec)
  v_sid = get_specimen_id(v_spec)

  have = u_matches.get(u_sid)   # (v_clas, [v_speq, ...])
  if have:
    (v_clas, v_speqs) = have
    if classified > v_clas:
      # Discard previously seen lower value matches.
      u_matches[u_sid] = (classified, [v_spec])
    elif classified == v_clas:
      if False and v_sid in map(get_specimen_id, v_speqs):
        # Two records with the same type specimen
        log("# Adding v_spec twice %s %s" %
            (explain_classified(v_clas),
             blorb(get_typifies(v_spec))))
      v_speqs.append(v_spec)    # ????
  else:
    u_matches[u_sid] = (classified, [v_spec])

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


# Compare potentially homotypic taxa in same workspace.

def ws_relate_records(u, v):
  return relate_records(get_outject(u), get_outject(v))

# In same checklist or in different checklists ... ?

def relate_records(x, y):
  comparison = compare_parts(get_parts(x),
                             get_parts(y))
  classified = classify_comparison(comparison)
  if monitor(x) or monitor(y):
    log("# Compare %s, %s = %s" %
        (blurb(x), blurb(y), explain_classified(classified)))
    (misses, hits) = comparison
    log("#  hits %s misses %s\n#  x %s\n#  y %s" %
        ((hits & TOKEN_MASK, hits & YEAR_MASK,), 
         (misses & TOKEN_MASK, misses & YEAR_MASK,),
         get_parts(x),
         get_parts(y),))

  if classified == MOTION:
    if not near_enough(x, y):
    # Proximity within the hierarchy is essential if no genus match
      return REVIEW
  return classified

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

def compare_parts(p, q):
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

# Compare Macropus robustus, Amblysomus robustus = heterotypic NOT

def classify_comparison(comp):
  (misses, hits) = comp

  # Epithets must match (perhaps '')
  if ((hits & EPITHET_MASK) == 0):
    return HETEROTYPIC          # Not a hit

  # If year and token BOTH differ then no match, yes?
  m = misses & (YEAR_MASK | TOKEN_MASK)
  if m == (YEAR_MASK | TOKEN_MASK):
    return HETEROTYPIC2          # Both misses

  # If EITHER year or token differs then requires review
  if m != 0:
    return REVIEW

  # Genus is only difference -> motion
  if (misses & GENUS_MASK) != 0:
    return MOTION

  return HOMOTYPIC

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
    word = "heterotypic (epithet)"
  elif classified == HETEROTYPIC2:
    word = "heterotypic (authority)"
  else:
    word = str(classified)
  return word

# Convenience.  Phase this out?  Or rename it?

def get_link(u, default=-19):
  uf = maybe_get_typification(u, None)
  if uf:
    (_, u2, v) = uf.payload()
    return v if (v and separated(u, v)) else u2
  return None
