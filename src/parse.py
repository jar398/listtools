#!/usr/bin/env python3

# Combine existing info from start.py output with parsed info from gnparse
# to obtain meaningful parts.
# Intended mainly for protonym sameness checks, which are for exemplar selection.

# Synthesize a string ('unparse') from parts in inverse operation.

import util, regex
from typing import NamedTuple, Any
from util import MISSING, log

class Parts(NamedTuple):
  scientific : Any              # what got parsed!
  canonical  : Any              # genus + middle + unstemmed
  genus      : Any
  epithet    : Any              # stem
  authorship : Any              # includes initials, all authors, etc
  token      : Any              # just Lastname of first author
  year       : Any
  protonymp  : Any              # False means moved, None means don't know.

# -----------------------------------------------------------------------------
# Parsing: string -> Parts

# Another way to split would be using the canonicalName and authorship
# columns of the original DwC - if they're present.

def parse_name(verbatim,
               gn_full=MISSING, gn_stem=MISSING, gn_auth=MISSING, \
               canonical=MISSING, authorship=MISSING):
  verbatim = verbatim.replace('  ', ' ')     # two -> one
  canonical = canonical.strip()
  authorship = authorship.strip()
  # If authorship is ill-formed, discard it e.g. Depéret, 1895 (Gervais, 1847)
  if (authorship.endswith(')') and
      not authorship.startswith('(')):
    authorship = MISSING
  if canonical != MISSING and authorship != MISSING:
    (canonical0, auth0) = (canonical, authorship)
  else:
    (canonical0, auth0) = split_name(verbatim)
  if auth0 != None:
    auth0 = auth0.strip(', ')
  assert auth0 == None or not auth0.startswith(',')
  auth_n = gn_auth.count(' ') + 1 if gn_auth != MISSING else 0
  if gn_full != MISSING:
    # Interesting to look at.
    #  e.g. 'Sigmodon dalquesti Stangl, Jr., 1992' != 'Sigmodon dalquesti' + 'Stangl & Jr. 1992'

    if False and gn_full.count(' ') + auth_n != verbatim.count(' '):
      log("** chunk count mismatch: %s + %s != %s" %
          (gn_full, gn_auth, verbatim))
    if gn_stem == MISSING: gn_stem = gn_full # Space saving kludge
    # Epithet will be stemmed
    if is_polynomial(gn_full):
      (canonical, g, e) = recover_canonical(gn_full, gn_stem, gn_auth, canonical0)
    else:
      (canonical, g, e) = (gn_full, gn_full, MISSING)
    auth = gn_auth
  else:
    canonical = canonical0
    auth = auth0
    # Epithet will not be stemmed
    if is_polynomial(canonical) and ' ' in canonical:
      chunks = canonical.split(' ')
      g = chunks[0]
      e = chunks[-1]              # do our own stemming !? no
    else:
      g = canonical
      e = MISSING

  if e == '?': e = None
  if g == '?': g = None         # cf. use_gnparse.py, but careful MSW3
  (t, y, protonymp) = analyze_authorship(auth)
  (t2, y2, protonymp2) = analyze_authorship(auth0)
  if not t: t = t2              # ???
  if t and t2 and t != t2:
    pass
    # log("# parse: tokens differ: gn %s / ad hoc %s" % (t, t2))
    # parse: tokens differ: De DeKay
    # parse: tokens (by different methods) differ: La LaVal
    # parse: tokens (by different methods) differ: Fitz FitzGibbon
  if not y: y = y2              # ???
  if y and y2 and y != y2:
    if False:
      log("# parse: gnparse year %s [%s] differs from ad hoc year %s [%s]" %
          (y, auth, y2, auth0))
    y = min(y, y2)
  # ugh, tangled code
  if protonymp == None: protonymp = protonymp2
  if protonymp != None and protonymp2 != None and \
     protonymp != protonymp2:
    log("# parse: protonymps differ: gn %s / ad hoc %s" %
        (auth, auth0))
  return Parts(verbatim, canonical, g, e, auth, t, y, protonymp)

# This ought to be trivial but the problem is that gnparser drops
# tokens that occur after the last alphabetic epithet.  So we have to
# recover them from the original non-GN scientific name ('verbatim').
# Returns (canonical, genus, epithet)

def recover_canonical(gn_full, gn_stem, gn_auth, hack_canonical):
  # Recover parts that gnparse stripped off, from verbatim
  hack_canonical = hack_canonical.strip()
  hack_canonical_chunks = hack_canonical.split(' ')
  n_hack_canonical_chunks = len(hack_canonical_chunks)

  n_full_chunks = gn_full.count(' ')
  # gnparser "Cryptotis avia  G. M. Allen, 1923 as C. thomasi & C. avia."
  if n_full_chunks > n_hack_canonical_chunks:
    # ** Ill formed canonical name: Homo sapiens × Rattus norvegicus fusion cell line
    log("** Ill formed canonical name: full %s\n / hack %s" % (gn_full, hack_canonical_chunks))
    extra = []
  else:
    # TBD: Should interpolate chunks that are in full but missing
    # from stemmed into stemmed  (wait, why?  have worked around this)
    extra = hack_canonical_chunks[n_full_chunks:]
  stem_chunks = gn_stem.split(' ')
  if len(extra) > 0:
    e = extra[-1]
    c = hack_canonical
  elif ' ' in gn_stem:
    e = stem_chunks[-1]
    c = gn_full
  else:
    e = ''
    c = gn_full
  g = stem_chunks[0]
  return (c, g, e)

# Split complete into genus, middle, epithet

def analyze_authorship(auth):
  if not auth: return (None, None, None) # moved? don't know
  if auth.endswith(').'): auth = auth[0:-1] # for MDD

  # GBIF sometimes has (yyyy) where it should have yyyy
  m = parenyear_re.search(auth)
  if m:
    auth = auth.replace(m[0], m[0][1:-1])

  # Find (...) in auth, assume it means this isn't a protonym, and flush parents
  if '(' in auth and ')' in auth:
    auth = auth.replace('(','')
    auth = auth.replace(')','')
    protonymp = False
  else:
    protonymp = True

  t_match = token_re.search(auth)
  tok = t_match[0] if t_match else None
  if tok == 'Someone': tok = None         # Wild card
  y_match = year_re.search(auth)
  year = y_match[1] if y_match else None
  if year == '1111': year = None       # Wild card
  return (tok, year, protonymp)

# This may not agree with gnparse exactly...
# [1] is complete canonical; [2] is paren or empty; [3] is authorship part (name(s) + maybe year)
LP = "\\("
RP = "\\)"
split_re = \
  regex.compile(u"(.+?) (((%s)?)\p{Uppercase_Letter}[\p{Letter}.'-]+.*)$" % LP)
token_re = regex.compile(u"\p{Uppercase_Letter}[\p{Letter}{3}-]+")
# Loses with Van Beneden & Gervais, 1868-79
year_re_string = '\\b([12][0-9]{3})\\b'
year_re = regex.compile(year_re_string)
#year_re = regex.compile(' ([12][0-9]{3})\)?$')
starts_auth_re = \
  regex.compile(u"((%s)?)\p{Uppercase_Letter}[\p{Letter}.'-]" % LP)
parenyear_re = \
  regex.compile("%s%s%s" % (LP, year_re_string, RP))

# Split verbatim using regex.  This is a parsing strategy independent
# of what gnparse does.
# Assume that there is no authority for hybrids!

def split_name(verbatim):
  if verbatim.endswith(')') and '(' in verbatim:
    i = verbatim.index('(')
    hack_canonical = verbatim[0:i].strip()
    hack_auth = verbatim[i:]
  else:
    m = split_re.search(verbatim)
    if m:
      hack_canonical = m[1].strip()
      hack_auth = m[2].strip()
    else:
      hack_canonical = verbatim
      hack_auth = None      # Not present
  return (hack_canonical, hack_auth)

# parsable into [genus] [middle] [epithet]

def is_polynomial(name):
  return not (' × ' in name or ' x ' in name)

# -----------------------------------------------------------------------------
# Unparsing = inverse of parsing = string from parts

# '' and None are both falsish, so be a bit careful

# Umm... really ought to just combine canonical and authorship

def unparse_parts(parts):
  (v, c, g, e, a, tok, y, protonymp) = parts
  assert g

  # Species name
  if ep == None:
    # Unknown genus
    canon = "%s ?" % g
  elif ep == '':
    # Higher taxon, genus or above
    canon = g  # Genus epithet
  else:
    canon = "%s %s" % (g, ep)  # Genus epithet

  # Authorship
  # y and tok are never '' but can be None
  if y == None and tok == None:
    auth = None
  else:
    if tok == None:
      tok = "Someone"
    if y == None:
      auth = tok                  # Jones
    else:
      auth = "%s, %s" % (tok, y) # Jones, 2088
    if protonymp == False: auth = "(%s)" % auth
    assert protonymp != None    # ???

  all = "%s %s" % (canon, auth) if auth != None else canon
  return all

# -----------------------------------------------------------------------------
# Parts - description of names or whatever with wildcard and so on -

def unify_parts(p1, p2):
  tup = tuple((unify_values(x, y) for (x, y) in zip(p1, p2)))
  return None if False in tup else tup

def unify_values(x, y):
  if x and y:
    return x if x == y else False
  else:
    return x or y               # None if neither

if __name__ == '__main__':
  import sys
  def qu(x): return "?" if x == None else "'%s'" % x

  parts = parse_name(sys.argv[1])
  (v, c, g, ep, a, tok, y, protonymp) = parts
  print("canonical %s" % qu(c))
  print(" genus %s" % qu(g))
  print(" epithet %s" % qu(ep))
  print("authorship %s" % qu(a))
  print(" token %s" % qu(tok))
  print(" year %s" % qu(y))
  print(" protonymp %s" % ('?' if protonymp == None else protonymp))
  print("recovered namestring '%s'" % unparse_parts(parts))
