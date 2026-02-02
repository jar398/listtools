#!/usr/bin/env python3

PROBE='zzzz'

# Combine existing info from start.py output with parsed info from gnparse
# to obtain meaningful parts.
# Intended mainly for protonym sameness checks, which are for exemplar selection.

# Synthesize a string ('unparse') from parts in inverse operation.

import util, regex
from util import MISSING, log, VERSION
from typing import NamedTuple, Any

class Parts(NamedTuple):
  scientific : Any              # what got parsed!
  canonical  : Any              # genus + middle + unstemmed
  genus      : Any
  middle     : Any              # stuff between genus and epithet
  epithet    : Any              # stem
  authorship : Any              # includes initials, all authors, etc
  token      : Any              # just Lastname of first author
  year       : Any
  protonymp  : Any              # False = moved (...), True = apparent protonym

# A missing (unknown) part is represented as None
WILD = None

# -----------------------------------------------------------------------------
# Parsing: string -> Parts

# Another way to split would be using the canonicalName and authorship
# columns of the original DwC - if they're present.

def parse_name(verbatim,
               gn_full=None, gn_stem=None, gn_auth=None, 
               canonical=None, authorship=None):
  verbatim = verbatim.replace('  ', ' ')     # two -> one
  if canonical == MISSING:
    canonical = None
  elif canonical:
    canonical = canonical.strip()
  if authorship == MISSING:
    authorship = gn_auth       # ?
  elif authorship:
    authorship = authorship.strip()
    if authorship.endswith('.'):
      authorship = authorship[0:-1] # for MDD

  # If authorship is ill-formed, discard it e.g. Depéret, 1895 (Gervais, 1847)
  if (authorship != None and
      authorship.endswith(')') and
      not authorship.startswith('(')):
    authorship = None
  if canonical != None and authorship != None:
    (canonical0, auth0) = (canonical, authorship)
  else:
    (canonical0, auth0) = split_name(verbatim)
  if auth0 != None:
    auth0 = auth0.strip(', ')

  # canonical0 and/or auth0 could be None
  assert auth0 == None or not auth0.startswith(',')
  if gn_auth == MISSING or gn_auth == None:
    gn_auth = None
    auth_n = 0
  else:
    auth_n = gn_auth.count(' ') + 1
  if gn_stem == MISSING:
    gn_stem = None
  if gn_full == MISSING:
    gn_full = None              # ???
  if gn_full == None:
    canonical = canonical0
    auth = auth0
    # Epithet will not be stemmed
    if is_polynomial(canonical) and ' ' in canonical:
      chunks = canonical.split(' ')
      g = chunks[0]
      mid = ' '.join(chunks[1:-1]) # maybe ''
      e = chunks[-1]              # do our own stemming !? no
    else:
      g = canonical
      # mid = '' if VERSION <= 1 else mid = WILD   ???
      mid = ''                  # no middle stuff
      e = ''                    # no epithet
    if PROBE in verbatim:
      log("# Parsed canonical '%s' as '%s' '%s'" % (canonical, g, e))
  else:
    # Epithet will be stemmed
    if is_polynomial(gn_full):
      (canonical, g, mid, e) = recover_canonical(gn_full, gn_stem, canonical0)
    else:
      (canonical, g, mid, e) = (gn_full, gn_full, '', '')
    auth = gn_auth
    if PROBE in verbatim:
      log("# Parsed stem '%s' as '%s' '%s'" % (gn_stem, g, e))


  if g == '?': g = WILD         # cf. use_gnparse.py, but careful MSW3
  if mid == '?': mid = WILD
  if e == '?': e = WILD
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
  # We never know something is a protonym, do we?
  # If authorship is missing altogether, it should be WILD.
  if protonymp == WILD: protonymp = protonymp2
  if protonymp != WILD and protonymp2 != WILD and \
     protonymp != protonymp2:
    # parse: protonymps differ: gn Hinton 1919 / ad hoc (Fitzinger, 1867) wroughtoni Hinton, 1919
    if False:   #too much of this in msw3
     log("# parse: protonymps differ: gn %s / ad hoc %s" %
         (auth, auth0))
    protonymp = protonymp2
  parts = Parts(verbatim, canonical, g, mid, e, auth, t, y, protonymp)
  if PROBE in verbatim: log("# Parts: %s" % (parts,))
  return parts

def duplicate_parts(p1, p2):
  return (p1.genus == p2.genus and
          p1.middle == p2.middle and
          p1.epithet == p2.epithet and
          p1.token == p2.token and
          p1.year == p2.year)

# Recover the canonical name and its part from the gnparser result.
# N.b. the returned epithet will be stemmed.

# This ought to be trivial but the problem is that gnparser drops
# tokens that occur after the last alphabetic epithet.  So we have to
# recover them from the original non-GN scientific name ('verbatim').
# Returns (canonical, genus, epithet).  Missing epithet is ''.

def recover_canonical(gn_full, gn_stem, hack_canonical):
  # Recover parts that gnparse stripped off, from verbatim
  if hack_canonical == None:
    hack_canonical_chunks = []
  else:
    hack_canonical = hack_canonical.strip()
    hack_canonical_chunks = hack_canonical.split(' ')
    n_hack_canonical_chunks = len(hack_canonical_chunks)

  n_full_chunks = gn_full.count(' ') + 1
  # gnparser "Cryptotis avia  G. M. Allen, 1923 as C. thomasi & C. avia."
  if n_full_chunks >= n_hack_canonical_chunks and gn_full and gn_stem:
    if n_full_chunks > n_hack_canonical_chunks:
      # ** Ill formed canonical name: Homo sapiens × Rattus norvegicus fusion cell line
      # This could be the symptom of an ill-formed authority
      # ** Using gn_full 'Sorex ? megalotis' not ad hoc 'Sorex ?'
      # ** Using gn_full 'Taphozous swirae sudani' not ad hoc 'Taphozous swirae'
      log("** canonicalName '%s' is truncated; gn_full is '%s'" %
          (hack_canonical, gn_full))
    c = gn_full
    # Assume gn_stem is non-None
    stem_chunks = gn_stem.split(' ')
    g = stem_chunks[0]
    if len(stem_chunks) > 1:
      mid = ' '.join(stem_chunks[1:-1]) # ?
      e = stem_chunks[-1]
    else:
      # Uninomial (e.g. genus, family)
      mid = ''
      e = ''
    if PROBE in e:
      log("# Recovered GN '%s' '%s' '%s'" % (g, mid, e))
    assert e or g, ("good", stem_chunks)
  else:
    # gnparser has dropped stuff off the end or out of middle.  Do not use.
    # This is not a good parse if name is not neolatinate.
    c = hack_canonical          # could be None
    g = hack_canonical_chunks[0]
    mid = ' '.join(hack_canonical_chunks[1:-1])
    e = hack_canonical_chunks[-1]
    # TBD: make use of gnparse's stemming
    if PROBE in hack_canonical:
      log("# Recovered ad hoc '%s' '%s' '%s'" % (g, mid, e))
      log("#  because '%s' '%s'" %
          (gn_full, hack_canonical))
    assert e or g, "bad"
  return (c, g, mid, e)

# Split authorship into token, year 
#  and note whether protonymp (not parenthesized)

def analyze_authorship(auth):
  if not auth: return (WILD, WILD, WILD) # moved? don't know

  # GBIF sometimes has (yyyy) where yyyy would be more helpful
  m = parenyear_re.search(auth)
  if m:
    auth = auth.replace(m[0], m[0][1:-1])

  # Find (...) in auth, assume it means this isn't a protonym, and flush parents
  if '(' in auth and ')' in auth:
    # Moved, not a protonym
    auth = auth.replace('(','')
    auth = auth.replace(')','')
    protonymp = False
  else:
    # Not moved, probably a protonym
    protonymp = True

  # Change W. E. Smith to just Smith.
  # Should also take care of von, etc.
  initial_match = initial_re.search(auth)
  if initial_match:
    auth2 = auth[len(initial_match[0]):]
    if False and auth2 != auth:
      log("# Trim %s -> %s" % (auth, auth2))
    auth = auth2

  t_match = token_re.search(auth)
  tok = t_match[0] if t_match else WILD
  if tok == 'Someone': tok = WILD         # Wild card, *, ?, Someone
  y_match = year_re.search(auth)
  year = y_match[1] if y_match else WILD
  if year == '1111': year = WILD       # Wild card

  if protonymp == WILD and tok != WILD and year != WILD:
    protonymp = True

  return (tok, year, protonymp)

# This may not agree with gnparse exactly...
# [1] is complete canonical; [2] is paren or empty; [3] is authorship part (name(s) + maybe year)
LP = "\\("
RP = "\\)"
split_re = \
  regex.compile(u"(.+?) (((%s)?)\\p{Uppercase_Letter}[\\p{Letter}.'-]+.*)$" % LP)

# Get token?  Consider 'Say in James, 1823' ~ 'Say, 1823'
initial_re = \
  regex.compile(u"(\\p{Uppercase_Letter}[.] )*")

# Pull out author's first capitalized (last) name part.
# Should actually be the last capitalized in first sequence.  or something.
# except if it's a Chinese or Vietnamese etc. name.
# 'Muntiacus jiangkouensis Gu Yonghe & Xu Longhui, 1998'
# '(Ho Hsi J., 1935)'
if VERSION <= 1:                # Four letter minimum
  token_re = regex.compile(u"\\p{Uppercase_Letter}[\\p{Letter}-]{3,}")
else:                           # Two letter minimum
  token_re = regex.compile(u"\\p{Uppercase_Letter}[\\p{Letter}-]{1,}")

# Loses with Van Beneden & Gervais, 1868-79
year_re_string = '\\b([12][0-9]{3})\\b'
year_re = regex.compile(year_re_string)
#year_re = regex.compile(' ([12][0-9]{3})\\)?$')
starts_auth_re = \
  regex.compile(u"((%s)?)\\p{Uppercase_Letter}[\\p{Letter}.'-]" % LP)
parenyear_re = \
  regex.compile("%s%s%s" % (LP, year_re_string, RP))

# Split verbatim using regex.  This is a parsing strategy independent
# of what gnparse does.
# Assume that there is no authority for hybrids!
# Missing parts are returned as WILD.

def split_name(verbatim):
  if verbatim.endswith(')') and '(' in verbatim:
    i = verbatim.index('(')
    hack_canonical = verbatim[0:i]
    hack_auth = verbatim[i:]
  else:
    m = split_re.search(verbatim)
    if m:
      hack_canonical = m[1]
      hack_auth = m[2]
    else:
      hack_canonical = verbatim
      hack_auth = None      # Not present
  if hack_canonical != None:
    hack_canonical = hack_canonical.strip()
    if len(hack_canonical) == 0: hack_canonical = None
  if hack_auth != None:
    hack_auth = hack_auth.strip()
    if len(hack_auth) == 0: hack_auth = None
  return (hack_canonical, hack_auth)

# parsable into [genus] [middle] [epithet]

def is_polynomial(name):
  return not (' × ' in name or ' x ' in name)

# -----------------------------------------------------------------------------
# Unparsing = inverse of parsing = string from parts

# '' (MISSING) and None (WILD) are both falsish, so be a bit careful

# Umm... really ought to just combine canonical and authorship

def unparse_parts(parts):
  (v, c, g, mid, ep, a, tok, y, protonymp) = parts
  assert g

  # Species name
  if ep == WILD:
    # Unknown genus
    canon = "%s ?" % g
  elif ep == '':
    # Higher taxon, genus or above
    canon = g  # Genus epithet
  elif mid == '':
    canon = "%s %s" % (g, ep)  # Genus epithet
  else:
    canon = "%s %s %s" % (g, mid, ep)  # Genus middle epithet

  # Authorship
  # y and tok are never '' but can be WILD
  if y == WILD and tok == WILD:
    auth = WILD
  else:
    if tok == WILD:
      tok = "?"
    if y == WILD:
      auth = tok                  # Jones
    else:
      auth = "%s, %s" % (tok, y) # Jones, 2088
    if protonymp == False: auth = "(%s)" % auth

  all = "%s %s" % (canon, auth) if auth != None else canon
  return all

if __name__ == '__main__':
  import sys
  def qu(x): return "?" if x == WILD else "'%s'" % x

  parts = parse_name(sys.argv[1])
  (v, c, g, mid, ep, a, tok, y, protonymp) = parts
  print("canonical %s" % qu(c))
  print(" genus %s" % qu(g))
  print(" middle %s" % qu(mid))
  print(" epithet %s" % qu(ep))
  print("authorship %s" % qu(a))
  print(" token %s" % qu(tok))
  print(" year %s" % qu(y))
  print(" protonymp %s" % ('?' if protonymp == None else protonymp))
  print("recovered namestring '%s'" % unparse_parts(parts))
