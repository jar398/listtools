# Name-based homotypy estimation

#  (A homotypy is an implication of a shared protonym, and protonyms
#  are 1-1 with type specimens)

import util
import simple

from util import log, MISSING, VERSION
from checklist import get_parts, get_rank, get_inferiors, \
  get_canonical, blurb, blorb, monitor, get_outject, \
  get_taxonomic_status
from rcc5 import DISJOINT
from specimen import equate_typicals, equate_specimens
from proximity import near_enough
from workspace import separated


# Name comparison outcomes.  Only the ordering matters.

# All in
DECLARED      = 11
HOMOTYPIC     = 10
CORRECTION    = 9   # token or year change, otherwise OK
MIDDLE        = 8.5
MOTION        = 8   # genus change: depends on topology etc.
CLASH         = 7   # genus change not annotated as such with '(' authority ')'
ALL_IN        = 6   # ---- all in above

# Equivocal
BINOMIAL      = 5   # no token and no year (binomial) (was AUTHORLESS)

# No way
NO_WAY        = 4   # ---- no way below
REVIEW        = 2   # token+year agree but genera don't
DISAGREE      = 1   # token and year both disagree
HETEROTYPIC   = 0
TOPO          = -1  # cannot match due to topological constraints

# Identify specimens in A checklist (swap to get B checklist).
# Several taxa (synonyms, or one a descendant of the other) might share
# a specimen.
# Called from find_endohomotypics in typify.py.

# It might be possible to make these matches by indexing by epithet.
# But this seems easier.

def find_homotypics_in_checklist(AB):
  if True:
    find_homotypics_in_checklist_v1(AB)
  else:
    find_homotypics_in_checklist_v2(AB.A)
    find_homotypics_in_checklist_v2(AB.B)

def find_homotypics_in_checklist_v1(AB):
  def process(x, epithets):     # x in A
    if get_rank(x, None) == "genus":
      epithets = {}             # epithet -> taxon in checklist
    if epithets != None:
      for c in get_inferiors(x):
        # assume children / accepted are found first
        z1 = AB.in_left(c)
        ep = get_parts(c).epithet or get_canonical(c)  # ???
        if ep in epithets:    # Seen before?
          # Previous time epithet has been encountered - same?
          c2 = epithets[ep]
          z2 = AB.in_left(c2)
          rel = simple.compare_per_checklist(c, c2)
          if rel.relation == DISJOINT:
            #log("# Keeping %s apart from %s; mrca = %s" %
            #    (blurb(z1), blurb(z2), blurb(simple.mrca(c, c2))))
            pass
            # Later-encountered does not get an exemplar; type specimens
            # are not unified.
          elif compare_record_protonyms(z1, z2) >= MOTION:
            # Found via tree traversal, so motion is OK ... ?
            equate_typicals(z1, z2)
        else:
          epithets[ep] = c
    for c in get_inferiors(x):
      process(c, epithets)
  process(AB.A.top, None)

# Not in workspace

def find_homotypics_in_checklist_v2(C):
  def process(x, prev):            # b is x's type taxon, so far
    b = prev or x
    for c in get_inferiors(x):  # sorted ?
      a = process(c, b)
      # Is c a "better" type for x than b?
      if compare_record_protonyms(c, x) >= ALL_IN:   # excludes disjoint
        if b is x:
          b = c                 # b is better, deeper
        else:
          # b is no better and maybe worse, stick with c.
          # Different taxon, should have same protonym as c
          assert compare_record_protonyms(c, b) >= ALL_IN   # not assured
        equate_types(c, b)   # b forwards to c
    return b
  process(C.top, None)

# ---------- compare_record_protonyms
# 
# Make an effort to determine how the two records' protonyms would
# compare, if we only knew them.

# Return value >= ALL_IN corresponds to INTERSECT.
# Return value <= NO_WAY corresponds to DISJOINT.

def compare_record_protonyms(u, v):
  # TBD:  If in same checklist, use simple.compare to filter.
  #  (that's probably useless)
  # TBD: use synonym/accepted status for some decisions.
  # TBD: 'homotypic synonym'
  x = get_outject(u)
  y = get_outject(v)
  c1 = compare_names_1(get_parts(x), get_parts(y))
  if VERSION <= 1:
    return c1
  else:
    c2 = None
    if not separated(u, v):     # In same checklist?
      # Check against Darwin Core semantics
      if not simple.simple_intersect(x, y):
        c2 = TOPO
      elif objective_synonym_of(x, y) or objective_synonym_of(y, x):
        c2 = DECLARED

    if c2 == None:
      c2 = compare_names_2(get_parts(x), get_parts(y))
    if c2 != c1:
      if c1 >= ALL_IN and c2 >= ALL_IN:
        pass # concur
      elif c1 <= NO_WAY and c2 <= NO_WAY:
        pass
      else:
        log("# Change from v1 to v2, %s -> %s\n  %s\n  %s" %
              (c1, c2, blorb(u), blorb(v)))
    return c2

# Wait, if u and v are in different checklists, this works only when y is
# x's superior's equivalent.  And we have no way to do that.  Foo.

def objective_synonym_of(x, y):
  return ("objective" in get_taxonomic_status(x, MISSING) and
          get_superior(x).record is y)

WILD = None # cf. parse.py

# Compare potentially homotypic names by examining their syntactic parts.
# Compare Macropus robustus, Amblysomus robustus = heterotypic NOT
# u and v are in the workspace, from different checklists.

# VERSION > 1

def compare_names_2(u_parts, v_parts):
  mismatches = 0     # Number of reasons to think they're het
  wilds = 0          # How musch wildcard reliance in authority?  0, 1, 2

  # EPITHET: mismatch implies taxon mismatch, deal breaker.
  u_ep = u_parts.epithet
  v_ep = v_parts.epithet
  if u_ep != WILD and v_ep != WILD:
    if u_ep != v_ep:
      # mismatches += 1
      return HETEROTYPIC
    else:
      wilds += 1

  # -- Canonicalize epithet and middle
  u_ge = u_parts.genus
  v_ge = v_parts.genus
  # Hack to force treating 'G' like 'G g', so that 'G' doesn't 
  # match both 'G x' and 'G y'
  if u_ep == MISSING and u_ge != WILD: u_ep = u_ge.lower()
  if v_ep == MISSING and v_ge != WILD: v_ep = v_ge.lower()

  # Make 'Foo baz' match 'Foo baz baz' but not 'Foo bar baz' ...right?
  # Or both?  Complicated.
  u_mid = u_parts.middle
  v_mid = v_parts.middle
  # Hack to force treating 'G e' like 'G e e' - is this right??
  # Depends on whether x and y are synonyms, yes?
  # Tricky logic here, be careful.
  if u_mid == MISSING and u_ep != WILD: u_mid = u_ep  # was WILD
  if v_mid == MISSING and v_ep != WILD: v_mid = v_ep  # was WILD


  # -- First peel off the no-ways --

  # AUTHOR
  u_tok = u_parts.token
  v_tok = v_parts.token
  if u_tok != WILD and v_tok != WILD:
    if u_tok != v_tok:
      mismatches += 1
  else:
    wilds += 1

  # YEAR
  u_yr = u_parts.year
  v_yr = v_parts.year
  if u_yr != WILD and v_yr != WILD:
    if u_yr != v_yr:
      mismatches += 1
  else:
    wilds += 1

  binomial = (wilds >= 2)       # No authority field

  # Mild to good agreement on the authority, look at the polynomial

  mid_change = False
  pre_epithet_change = False

  # MIDDLE - complicated, play around with this.
  if u_mid != WILD and v_mid != WILD:
    if u_mid != v_mid:
      mid_change = True
      pre_epithet_change = True
  else:
    pass

  # GENUS: Protonym mismatch implies mismatch

  # Three genus comparison cases: clash < motion < same
  clash = motion = False
  if u_ge != WILD and v_ge != WILD:
    if u_ge != v_ge:
      if u_parts.protonymp and v_parts.protonymp:   # No parens
        clash = True                                # even worse
      else:
        motion = True           # Foo x Smith -> Bar x (Smith) etc.
      pre_epithet_change = True
  else:
    pass

  if pre_epithet_change:
    mismatches += 1

  # Genus+epithet match, auth+year no mismatch...
  if mismatches > 1:
    return DISAGREE

  # --- End no-ways, now distinguish special situations from all-ins
  # Often we have combinations of these conditions

  # Risky
  if binomial: return BINOMIAL  # Foo x could be Foo x Smith or Foo x Jones
  if clash: return CLASH        # motion from Foo Smith to Bar (Smith) etc.

  # Not so bad
  if motion: return MOTION      # motion from Foo Smith to Bar Smith etc.

  # Downright benign
  if mid_change: return MIDDLE

  # Not bad at all
  if mismatches > 0: return CORRECTION

  return HOMOTYPIC     # All parts match


# VERSION == 1

def compare_names_1(name1, name2):
  (misses, hits) = parts_comparison_detail(name1, name2)
  return classify_comparison_details(misses, hits)

# Given a comparison of parts, classify it as appropriate

# Parts masks

EPITHET_MASK = 32
VICINITY_MASK = 16
YEAR_MASK = 8
TOKEN_MASK = 4
GENUS_MASK = 2
MIDDLE_MASK = 1

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
      return REVIEW  # v1 !!????

  elif misses & mask == 0:
    return BINOMIAL

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

# VERSION > 1



# not used?

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
  elif classified == CORRECTION:
    word = "correction"
  elif classified == MOTION:
    word = "motion"           # kinetypic?
  elif classified == BINOMIAL:
    word = "authorless"
  elif classified == REVIEW:
    word = "review"
  elif classified == HETEROTYPIC:
    word = "heterotypic"
  else:
    word = str(classified)
  return word

# --------------------

# Compare potentially homotypic taxa in same workspace.

def relate_records(u, v):
  classified = compare_record_protonyms(u, v)

  if classified == MOTION:
    if not near_enough(u, v):
    # Proximity within the hierarchy is essential if no genus match
      return REVIEW

  if monitor(u) or monitor(v):
    log("# Compare %s, %s = %s" %
        (blurb(u), blurb(v), explain_classified(classified)))
  return classified

