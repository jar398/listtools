#!/usr/bin/env python3

# Combine existing info from start.py output with parsed info from gnparse
# to obtain protonym parts.
# Intended mainly for protonym sameness checks, which are for exemplar selection.

# Synthesize a string ('unparse') from protonym parts in inverse operation.

import util, regex
from util import MISSING, log

class Parts(NamedTuple):
  complete  : Any               # for comparisons  ... ?
  genus     : Any
  epithet   : Any
  authority : Any               # for comparisons ... ?
  token     : Any
  year      : Any
  moved     : Any

class Protonymic(NamedTuple):
  epithet : Any               # None means taxon is a genus
  seenInGenus : Any           # might be same as .genus
  token   : Any
  year    : Any
  # maybe others
  protoGenus   : Any          # Genus of original coinage

# -----------------------------------------------------------------------------
# Parsing: string -> Protonymic

def parse_name(verbatim, **more):
  (c, a) = split_name(verbatim, **more)
  (g, _, e) = analyze_canonical(c)
  (t, y, moved) = analyze_authorship(a)
  return Parts(c, g, e, a, t, y, moved)

def parse_protonymic(verbatim, **more):     # returns a 'proto'
  (_, g, e, _, t, y, moved) = parse_name(verbatim, **more)
  cg = g
  if moved: g = None                # moved
  return Protonymic(g, e, t, y, cg)

# Returns (complete, authorship)

def split_name(verbatim, gn_full=None, gn_authorship=None, gn_stemmed=None):
  assert verbatim

  # There are two ways to split the verbatim name into complete + authorship.
  # One is by using the gnparse output, filling in missing in-between
  # segment from verbatim.  The other is by using a regular expression
  # to parse verbatim directly.

  (hack_canonical, hack_auth) = hack_split(verbatim)

  # For low-quality names, gnparser may throw away authorship entirely.

  # 2. Provisional verbatim split using gnparse (preferred)
  if gn_full != None and gn_authorship != None:

  (gn_canonical, gn_auth) = gn_split(verbatim, gn_full, gn_stemmed, gn_auth)

  # May need to revise this.
  assert (gn_auth != None) == (hack_auth != None)

  if gn_full != None:
    return (gn_canonical, gn_auth)
  else:
    return (hack_canonical, hack_auth)

def hack_split(verbatim):
  # 1. Provisional verbatim split using regex
  m = split_re.search(verbatim)
  if m:
    hack_canonical = m[1].strip()
    hack_auth = m[2].strip()
  else:
    hack_canonical = verbatim
    hack_auth = None      # Not present
  return (hack_canonical, hack_auth)

def gn_split(verbatim, gn_full, gn_stemmed, gn_auth):

  # Recover parts that gnparse stripped off, from verbatim
  n_full_parts = gn_full.count(' ')
  verbatim = verbatim.strip()
  verbatim_parts = verbatim.split(' ')
  n_verbatim_parts = len(verbatim_parts)
  n_auth_parts = gn_auth.count(' ')
  assert n_full_parts + n_auth_parts <= n_verbatim_parts
  middle = verbatim_parts[n_full_parts: auth_pos]
  # TBD: Should map from full to stemmed
  if gn_stemmed:
    n_stemmed_parts = gn_stemmed.count(' ')
    if n_full_parts == n_stemmed_parts:
      stemmed = gn_full
    else:
      assert "elaborating canonicalStem from canonicalFull NYI", gn_full

  gn_canonical = ' '.join([stemmed] + middle)
  assert auth_pos == n_verbatim_parts - n_auth_parts
  assert starts_auth_re.match(parts[auth_pos])
  return (gn_canonical, gn_auth)

# Split complete into genus, middle, epithet

def analyze_canonical(complete, stemmed=None):
  parts = complete.split(' ')
  g = parts[0]
  if len(parts) > 1:
    if g == '?' or g == 'Nil': g = None
    return (g, ' '.join(parts[1:-1]), parts[-1])
  else:
    return (g, None, parts[0])       # could be genus, family etc.

def analyze_authorship(auth):
  if not auth: return (None, None, True)
  moved = (auth[0] == '(' and auth[-1] == ')')
  if moved:
    auth = auth[1:-1]
  t = token_re.search(auth)
  token = t[1] if t else None
  if token == "Someone": token = None
  y = year_re.search(auth)
  year = y[1] if y else None
  return (token, year, moved)

# This may not agree with gnparse exactly...
# [1] is complete canonical; [2] is paren or empty; [3] is authorship part (name(s) + maybe year)
LP = "\\("
split_re = \
  regex.compile(u"(.+?) (((%s)?)\p{Uppercase_Letter}[\p{Letter}-]+.*)$" % LP)
token_re = regex.compile(u"(\p{Uppercase_Letter}[\p{Letter}-]+)")
year_re = regex.compile(' ([12][0-9]{3})\)?$')
starts_auth_re = \
  regex.compile(u"((%s)?)\p{Uppercase_Letter}\p{Letter}" % LP)

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
  print(" seenInGenus '%s'" % (y or "no year"))
  print("protonymic '%s'" % unparse_proto(proto))
