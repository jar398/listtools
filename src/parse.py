#!/usr/bin/env python3

# Combine existing info from start.py output with parsed info from gnparse
# to obtain protonym parts.
# Intended mainly for protonym sameness checks, which are for exemplar selection.

# Synthesize a string ('unparse') from protonym parts in inverse operation.

import util, regex
from util import MISSING, log

class Protonymic(NamedTuple):
  epithet : Any                 # None means taxon is a genus
  genus   : Any                 # None means moved
  token   : Any
  year    : Any
  currentGenus : Any            # might be same as .genus
  # maybe others


# -----------------------------------------------------------------------------
# Parsing: string -> Protonymic

def parse_protonymic(verbatim, **more):     # returns a 'proto'
  (c, a) = split_protonymic(verbatim, **more)
  (g, _, e) = analyze_complete(c)
  (t, y, moved) = analyze_authorship(a)
  cg = g
  if moved: g = None                # moved
  return Protonymic(g, e, t, y, cg)

# This may not agree with gnparse exactly...
# [1] is complete canonical; [2] is paren or empty; [3] is authorship part (name(s) + maybe year)
LP = "\\("
split_re = \
  regex.compile(u"(.+?) (((%s)?)\p{Uppercase_Letter}[\p{Letter}-]+.*)$" % LP)
token_re = regex.compile(u"(\p{Uppercase_Letter}[\p{Letter}-]+)")
year_re = regex.compile(' ([12][0-9]{3})\)?$')
starts_auth_re = \
  regex.compile(u"((%s)?)\p{Uppercase_Letter}\p{Letter}" % LP)

# Returns (complete, authorship)

def split_protonymic(verbatim, gn_full=None, gn_authorship=None, gn_stemmed=None):
  assert verbatim

  # There are two ways to split the verbatim name into complete + authorship.
  # One is by using the gnparse output, filling in missing in-between
  # segment from verbatim.  The other is by using a regular expression
  # to parse verbatim directly.

  # 1. Provisional verbatim split using regex
  m = split_re.search(verbatim)
  if m:
    hack_complete = m[1].strip()
    hack_authorship = m[2].strip()
  else:
    hack_complete = verbatim
    hack_authorship = None      # Not present

  # For low-quality names, gnparser may throw away authorship entirely.

  # 2. Provisional verbatim split using gnparse (preferred)
  if gn_full != None and gn_authorship != None:

    # Recover parts that gnparse stripped off, from verbatim
    n_full_parts = gn_full.count(' ')
    verbatim = verbatim.strip()
    verb_parts = verbatim.split(' ')
    n_verb_parts = len(verb_parts)
    n_auth_parts = gn_auth.count(' ')
    assert n_full_parts + n_auth_parts <= n_verb_parts
    middle = verb_parts[n_full_parts: auth_pos]

    # TBD: Should map from full to stemmed
    if gn_stemmed:
      n_stemmed_parts = gn_stemmed.count(' ')
      if n_full_parts == n_stemmed_parts:
        stemmed = gn_full
      else:
        assert "elaborating canonicalStem from canonicalFull NYI", gn_full

    gn_complete = ' '.join([stemmed] + middle)
    assert auth_pos == n_verb_parts - n_auth_parts
    assert starts_auth_re.match(parts[auth_pos])

  assert (gn_authorship != None) == (hack_authorship != None)

  if gn_full != None:
    return (gn_complete, gn_authorship)
  else:
    return (hack_complete, hack_authorship)

# Split complete into genus, middle, epithet

def analyze_complete(complete, stemmed=None):
  parts = complete.split(' ')
  if len(parts) > 1:
    g = parts[0]
    if g == '?' or g == 'Nil': g = None
    return (g, ' '.join(parts[1:-1]), parts[-1])
  else:
    return (g, None, parts[0])       # could be genus, family etc.

def analyze_authorship(authorship):
  if not authorship: return (None, None, False)
  moved = (authorship[0] == '(' and authorship[-1] == ')')
  if moved:
    authorship = authorship[1:-1]
  t = token_re.search(authorship)
  token = t[1] if t else None
  if token == "Someone": token = None
  y = year_re.search(authorship)
  year = y[1] if y else None
  return (token, year, moved)

# -----------------------------------------------------------------------------
# Unparsing = inverse of parsing = string from protonymic

# '' and None are both falsish, so be a bit careful

def unparse_proto(proto):
  (g, ep, tok, y) = proto
  assert ep

  # Species name
  assert g != None
  if ep == None:
    canon = g
  else:
    assert ep
    canon = "%s %s" % (g, ep)  # Genus epithet

  # Authorship
  # y and tok are never '' but can be None
  if y == None and tok == None:
    auth = None
  else:
    if tok == None:
      tok = "Someone"
    if y == None:
      auth = t                  # Jones
    else:
      auth = "%s, %s" % (tok, y) # Jones, 2088

  all = "%s %s" % (canon, auth) if auth != None else canon
  return "(%s)" if g == None else all

# -----------------------------------------------------------------------------
# Protonymic - description of protonyms with wildcard and so on -

def unify_protos(p1, p2):
  tup = tuple((unify_values(x, y) for (x, y) in zip(p1, p2)))
  return None if False in tup else tup

def unify_values(x, y):
  if x and y:
    return x if x == y else False
  else:
    return x or y               # None if neither

if __name__ == '__main__':
  import sys
  proto = parse_protonymic(sys.argv[1])
  (g, ep, tok, y) = proto
  print(" genus '%s'" % (g or '?'))
  print(" epithet '%s'" % ep)
  print(" token '%s'" % (tok or 'Someone'))  
  print(" year '%s'" % (y or "no year"))
  print("protonymic '%s'" % unparse_proto(proto))
