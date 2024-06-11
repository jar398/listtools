from parse import PROBE, duplicate_parts

import util
import property as prop
import simple

from util import log, MISSING, UnionFindable
from rcc5 import EQ, NEQ, OVERLAP, DISJOINT
from checklist import get_parts, monitor, \
  get_superior, get_children, get_synonyms, get_inferiors, \
  has_inferiors, \
  is_accepted, blurb, blorb, get_scientific, get_primary_key, \
  get_tag, get_source_tag, get_nomenclatural_status, \
  get_children, get_synonyms, get_accepted, \
  get_suppressed, set_suppressed
from workspace import separated, get_outject, get_source
from specimen import sid_to_epithet, same_specimens, \
  equate_typifications, equate_specimens, same_typifications, \
  get_typification, get_specimen_id, \
  sid_to_specimen, get_typifies
from proximity import near_enough

def match_typifications(AB):
  A_index = prepare(AB.A, lambda x: AB.in_left(x))
  B_index = prepare(AB.B, lambda y: AB.in_right(y))
  subproblems = collate_subproblems(A_index, B_index)
  log("* Matching types between checklists:")
  match_in_subproblems(subproblems)
  if False:
    find_estimates(AB)            # for distance calculations
    log("* Refining type matches using distance calculations:")
    match_in_subproblems(subproblems)

# -----------------------------------------------------------------------------
# STEPS 1 and 2.

# Prepare a single checklist for matching.
# Doesn't feel like prepare really belongs in this file.  Where does it go?

def prepare(checklist, inject):
  log("* Preparing checklist %s" % get_tag(checklist))
  basic_typifications(checklist, inject)
  def get_key(x):
    return get_subproblem_key(x)
  index = index_subtree(checklist,
                        get_key,
                        inject)
  suppress_extras(index)
  log("")
  return index

# Each subproblem covers a single epithet (or name, if higher taxon)
# z is in AB

def get_subproblem_key(x):
  parts = get_parts(x)
  ep = parts.epithet            # stemmed
  key = ep if ep else parts.genus
  if key:
    pass
  else:
    log("** %s: Name missing or ill-formed: %s" %
        (get_primary_key(x), parts,))
    key = '?' + get_primary_key(x)
  return key

# Returns a dict key -> list of record

def index_subtree(checklist, get_key, inject):
  index = {}
  def traverse(x):
    key = get_key(x)
    u = inject(x)
    have = index.get(key, None) # part
    if have:
      have.append(u)
    else:
      index[key] = [u]
    for i in get_inferiors(x): traverse(i)
  traverse(checklist.top)
  for (key, us) in index.items():
    us.sort(key=unimportance)
  return index

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

# -----------------------------------------------------------------------------
# STEP 1.
# Find typification identities based on topology,
# starting from a given root node

def basic_typifications(checklist, inject):
  log("# Designating type specimens for %s checklist" %
      get_source_tag(checklist.top))
  equate=equate_typifications

  def process(x):
    winner = x
    candidates = []           # homotypic children
    for c in get_children(x, ()):
      process(c)
      if homotypy(c, x) >= HOMOTYPIC:
        # c nominates itself as the type of x
        candidates.append(c)
    n = len(candidates)
    if n == 1:
      winner = candidates[0]
      # Prefer child with same epithet
      equate(inject(winner), inject(x))
    for s in get_synonyms(x, ()):
      # N.b. synonyms do not have descendants
      process(s)              # Unnecessary because no inferiors
      if (homotypy(s, x) >= MOVED or
          'homotypic' in get_nomenclatural_status(s, '')):
        # s is a homotypic synonym of x (because of name)
        equate(inject(winner), inject(s))
  process(checklist.top)

# -----------------------------------------------------------------------------
# STEP 2.
# Find homotypic records and apparent errors within a single checklist.
# Duplicates (to skip) are indicated by get_suppressed.
# This probably needs to be rewritten.

def suppress_extras(index):
  n = 1
  for (key, us) in index.items():
    if n % 1000 == 0 or PROBE in key:
      log("# Suppressing extras %s %s %s" %
          (n, len(us), blurb(us[0])))
    n += 1

    # For each high-importance record in the given subproblem, find
    # other records that are putatively homotypic to it, and suppress
    # them.

    # *** Go by specimens or by records? ***
    # Well, since we end up suppressing records, not specimens, we should
    # go by records.

    seen = set()               # set of record primary keys ...
    u_extras = {}              # u_pk -> (u, class, [extra, ...])
    for i in range(0, len(us)):
      u = us[i]
      u_pk = get_primary_key(u)
      if u_pk in seen: continue
      seen.add(u_pk)
      for j in range(i+1, len(us)):
        u2 = us[j]
        if same_typifications(u, u2): continue
        classified = homotypy(u, u2)
        if classified >= REVIEW:
          # Report on this homotypy assessment
          if classified > MOVED and not has_inferiors(u2):
            u2_pk = get_primary_key(u2)
            if u2_pk in seen: continue
            seen.add(u2_pk)
            if get_superior(get_outject(u, u)).record is get_superior(get_outject(u2, u2)).record:
              if False:
                log("# Homotypic sibs: %s: %s ~ %s" %
                    (explain_classified(classified), designate(u), designate(u2)))
              equate_typifications(u, u2)
            else:
              log("# Suppressing: %s: %s ~ %s" %
                  (explain_classified(classified), designate(u), designate(u2)))
            set_suppressed(u2, u)
          else:
            log("# Review: %s: %s ~ %s" %
                (explain_classified(classified), designate(u), designate(u2)))
      # end j loop
    # end i loop

# -----------------------------------------------------------------------------

# STEP 3. Find specimen matches between checklists

def collate_subproblems(A_index, B_index):
  subprobs = {}
  for (key, us) in A_index.items():
    assert key != MISSING, blurb(us[0])
    vs = B_index.get(key, None)
    if vs != None:
      subprobs[key] = (us, vs)
  log("* There are %s subproblems." % len(subprobs))
  return subprobs

# Records are matched between checklists by updating their 'typification_uf' 
# settings.

def match_in_subproblems(subprobs):
  n = 1
  for (key, (us, vs)) in subprobs.items():  # For each subproblem
    if n % 1000 == 0 or PROBE == key:
      log("# Subproblem %s %s %s %s %s" %
          (n, len(us), len(vs), blurb(us[0]), blurb(vs[0])))
    n += 1
    match_in_subproblem(us, vs)

# Specimen/specimen matches are induced from record/record matches.
# We need to filter out ambiguities...
# But... it is hoped most or all ambiguities are filtered out at the
# 'suppression' stage.  So ideally this does nothing.
# candidates_for_u : u_sid -> (u_spec, v_clas, v_specs)

def match_in_subproblem(us, vs):
  (candidates_for_u, candidates_for_v) = \
    candidates_in_subproblem(us, vs, get_typification)
  half_match_in_subproblem(candidates_for_u, candidates_for_v)
  half_match_in_subproblem(candidates_for_v, candidates_for_u)

def half_match_in_subproblem(candidates_for_u, candidates_for_v):
  # candidates_for_u : u_sid -> (u_spec, class, [v_spec, ...])
  # candidates_for_v : v_sid -> (v_spec, class, [u_spec, ...])
  # Test every candidate to see if it is a unique match.
  for (u_spec, v_clas, v_specs) in candidates_for_u.values():
    if len(v_specs) > 1:
      # Get rid of synonyms
      valids = []
      for v_spec in v_specs:
        v = get_typifies(v_spec)
        if is_accepted(get_outject(v)):
          valids.append(v_spec)
        else:
          log("# Ambiguity via synonym (%s): %s -> %s" %
              (explain_classified(v_clas),
               designate(get_typifies(u_spec)),
               designate(get_typifies(v_spec))))
      v_specs = valids          # Flush synonyms
    # Cf. THIS_STUFF
    if len(v_specs) > 1:
      for v_spec in v_specs:
        v = get_typifies(v_spec)
        #if not has_inferiors(v):
        log("# Ambiguity via accepted (%s): %s -> %s" %
              (explain_classified(v_clas),
               designate(get_typifies(u_spec)),
               designate(get_typifies(v_spec))))
    if len(v_specs) == 1:
      v_spec = v_specs[0]
      u_specs = candidates_for_v.get(get_specimen_id(v_spec))
      assert u_specs
      if len(u_specs) == 1:
        assert u_specs[0] is u_spec
        equate_specimens(u_spec, v_spec)


# We are matching specimens to specimens; any such match induces some
# number of record/record matches, via get_specimen and get_typifies.

# A candidate match (v_spec) for a record (u) is a homotypy match to
# any of that match's records.

def candidates_in_subproblem(us, vs, get_specimen):
  candidates_for_u = {}    # u_sid -> (u_spec, class, [v_spec, ...])
  candidates_for_v = {}    # v_sid -> (v_spec, class, [u_spec, ...])
  for i in range(0, len(us)):
    u = us[i]
    if get_suppressed(u, None): continue
    if monitor(u): log("# Subproblem row: '%s'" % (blorb(u),))
    for j in range(0, len(vs)):
      v = vs[j]
      if get_suppressed(v, None): continue
      classified = homotypy(u, v)
      if classified > REVIEW:
        u_spec = get_specimen(u)
        v_spec = get_specimen(v)
        observe_homotypy(u_spec, v_spec, candidates_for_u, classified)
        observe_homotypy(v_spec, u_spec, candidates_for_v, classified)
      elif classified == REVIEW:
        # TBD: Make a note of it!
        pass
    # end j loop
  # end i loop
  return (candidates_for_u, candidates_for_v)

# candidates_for_u : u_sid -> (u_spec, v_clas, v_specs)

def observe_homotypy(u_spec, v_spec, candidates_for_u, classified):
  u_sid = get_specimen_id(u_spec)

  have = candidates_for_u.get(u_sid)   # (v_clas, [v_spec, ...])
  if have:
    (_, v_clas, v_specs) = have
    if classified < v_clas:
      # Not so good as what we have already
      pass
    elif classified > v_clas:
      # Replace specs with new and shiny
      candidates_for_u[u_sid] = (u_spec, classified, [v_spec])
    else:
      v_spec = v_spec.find()    # uf
      if get_specimen_id(v_spec) in map(get_specimen_id, v_specs):
        pass
      else:
        # Tie - add to list
        v_specs.append(v_spec)
  else:
    candidates_for_u[u_sid] = (u_spec, classified, [v_spec])

# _____________________________________________________________________________

# This is not a match:
#   Tylonycteris malayana Chasen, 1940 -> Euroscaptor malayanus (Chasen, 1940)

HOMOTYPIC     = 10
CORRECTION    = 8
MOVED         = 6
REVIEW        = 5               # Not used any more?
DISTANT       = 1               # too far away to be MOVED
HETEROTYPIC   = 0

def explain_classified(typy):
  if typy == HOMOTYPIC:    word = "homotypic"
  elif typy == CORRECTION: word = "correction"
  elif typy == MOVED:      word = "mobile"           # kinetypic?
  elif typy == REVIEW:     word = "review"
  elif typy == DISTANT:    word = "distant"
  elif typy == HETEROTYPIC: word = "heterotypic"
  else: word = str(typy)
  return word

# x and y are in same checklist or same workspace

def homotypy(x, y):
  # TBD: what if x is y?  if (x is y) return IDENTICAL
  # TBD: consider matched existing typification if any?
  px = get_parts(x)
  py = get_parts(y)
  (misses, hits) = compare_parts(px, py)

  # No epithet or genus match -> mismatch
  if (hits & (EPITHET_MASK | GENUS_MASK)) == 0:
    clas = HETEROTYPIC          # ??
  elif misses == 0:
    clas = HOMOTYPIC
  elif misses == YEAR_MASK:
    clas = CORRECTION          # ??
  elif misses == TOKEN_MASK:
    clas = CORRECTION          # ??
  elif misses == GENUS_MASK:
    if not near_enough(get_outject(x, x),
                       get_outject(y, y)):
      clas = DISTANT
    elif px.protonymp and py.protonymp:
      # Both are protonyms - no go, but maybe we should review
      clas = REVIEW
    else:
      clas = MOVED
  else:
    clas = HETEROTYPIC          # ??
  return clas

# -----------------------------------------------------------------------------

EPITHET_MASK = 32
YEAR_MASK = 16
TOKEN_MASK = 8
GENUS_MASK = 4
PROTONYMP_MASK = 2
MIDDLE_MASK = 1

# Compare potentially homotypic names.  Returns (m1, m2) where m1 and
# m2 are integer masks, m1 for differences and m2 for similarities.

def compare_parts(p, q):
  hits = misses = 0

  if p.epithet != None and q.epithet != None:
    # If not epithet, use genus
    pep = p.genus.lower() if p.epithet == '' and p.genus else p.epithet
    qep = q.genus.lower() if q.epithet == '' and q.genus else q.epithet
    if pep == qep: hits |= EPITHET_MASK
    else: misses |= EPITHET_MASK

  if p.year != None and q.year != None:
    if p.year == q.year: hits |= YEAR_MASK
    else: misses |= YEAR_MASK

  if p.token != None and q.token != None:
    if p.token == q.token: hits |= TOKEN_MASK
    else: misses |= TOKEN_MASK

  if p.genus != None and q.genus != None:
    if p.genus == q.genus: hits |= GENUS_MASK
    else: misses |= GENUS_MASK

  if p.protonymp != None and q.protonymp != None:
    if p.protonymp == q.protonymp: hits |= PROTONYMP_MASK
    else: misses |= PROTONYMP_MASK

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
                      in zip(range(0, 5),
                             (
                              "middle",   # MIDDLE_MASK = 1 << 1
                              "genus",    # GENUS_MASK = 2 << 1
                              "token",    # TOKEN_MASK
                              "year",     # YEAR_MASK
                              "epithet",  # EPITHET_MASK
                             ))
                      if (things & (1<<i)) != 0))

  (misses, hits) = comparison
  return ("hits(%s) misses(%s)" %
          (explode(hits), explode(misses)))
    
# -----------------------------------------------------------------------------

# Convenience.  Phase this out?  Or rename it?
# See link.py, risk.py

def get_link(u, default=-19):
  uf = maybe_get_typification(u, None)
  if uf:
    (_, u2, v) = uf.payload()
    return v if (v and separated(u, v)) else u2
  return None

def designate(x):
  return "%s %s" % (blorb(x), get_primary_key(get_outject(x, x)))

