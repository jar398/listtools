# Name-based homotype estimation

import util
import simple

from util import log
from checklist import get_parts, get_rank, get_inferiors, \
  get_canonical
from rcc5 import DISJOINT
from specimen import equate_type_ufs


# Name comparison outcomes

HOMOTYPIC     = 10
CORRECTION    = 9
MOTION        = 8

AUTHORLESS    = 5   # no token or year
REVIEW        = 2   # includes CORRECTION + MOTION
HETEROTYPIC   = 0

# Parts masks

EPITHET_MASK = 32
VICINITY_MASK = 16
YEAR_MASK = 8
TOKEN_MASK = 4
GENUS_MASK = 2
MIDDLE_MASK = 1

# Identify specimens in A checklist (swap to get B checklist).
# Several taxa (synonyms, or one a descendant of the other) might share a specimen.

def find_homotypics_in_checklist(AB):
  def process(x, epithets):     # x in A
    if get_rank(x, None) == "genus":
      epithets = {}
    if epithets != None:
      for c in get_inferiors(x):
        ep = get_parts(c).epithet or get_canonical(c)  # ???
        if ep in epithets:
          # Previous time epithet has been encountered - same?
          c2 = epithets[ep]
          z = AB.in_left(c)
          z2 = AB.in_left(c2)
          if simple.compare_per_checklist(c, c2) == DISJOINT:
            log("# Keeping %s apart from %s" %
                (blurb(z), blurb(z2)))
            pass
          else:
            # if obviously distinct then ... ?  e.g. year+token
            # compare_parts(z, z2) >= MOTION ...
            equate_type_ufs(z, z2)
        else:
          epithets[ep] = c
    for c in get_inferiors(x):
      process(c, epithets)
  process(AB.A.top, None)


# ---------- compare_parts

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

  elif misses & mask == 0:
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
