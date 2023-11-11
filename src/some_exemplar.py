#!/usr/bin/env python3

import math
import property as prop
import checklist, workspace, simple
import linkage
#import exemplar

from util import log, UnionFindable
from checklist import *
from workspace import *

# Exemplars (representing individual specimens or occurrences with
# known classification in BOTH checklists) are chosen heuristically
# via name matches.  Not every name match qualifies; only those that
# are 'tipward' do.  'Tipward' means that at least one of the two taxa
# has no tipward accepted descendant taxon.

# See AB.exemplar_ufs for map {xid: uf, ...}
# Each exemplar has a representative in A and one in B.

# find_exemplars is defined in exemplar.py.  Does 2 passes.
# find_some_exemplars is just one pass; invoked twice.

def find_some_exemplars(AB, subproblems, get_pre_estimate, allow_nontipward):
  log("# Finding some exemplars")
  linkage.really_find_links(AB, subproblems, get_pre_estimate)

  def do_exemplars(CD):
    def traverse(x):      # species is an ancestor of x, in A
      tipward = True
      for c in get_inferiors(x):
        if traverse(c): tipward = False
      if tipward:
        u = CD.in_left(x)
        v = get_link(u, None)
        if v:
          u2 = get_link(v, None)
          if u2 is u:
            # But only if tipward? or overlap??
            equate_exemplars(u, v)

    traverse(CD.A.top)

  do_exemplars(AB)
  do_exemplars(swap(AB))
  equate_exemplars(AB.in_left(AB.A.top),
                   AB.in_right(AB.B.top))

# Issues: 
#   1. choice of taxa (both sides)
#   2. choice of name

def equate_exemplars(u, v):     # opposite checklists. u might be species
  if u != v:
    # assert u descends from v if in same checklist?
    equate_exemplar_ufs(get_exemplar_uf(u), get_exemplar_uf(v))
  return u

def equate_exemplar_ufs(uf, vf):
  uf = uf.find()
  vf = vf.find()
  (i1, u1, v1) = uf.payload()
  (i2, u2, v2) = vf.payload()
  assert i1 == None or i2 == None or i1 == i2
  assert u1 or v1
  assert u2 or v2
  ef = uf.absorb(vf)
  assert ef is uf
  r = ef.payload()
  r[1] = pick_better_record(u1, u2)
  r[2] = pick_better_record(v1, v2)
  assert r[1] or r[2]
  return ef

# u1 and u2 are in workspace, in same checklist

def pick_better_record(u1, u2):
  if not u2: return u1
  if not u1: return u2
  x1 = get_outject(u1)
  x2 = get_outject(u2)
  rel = simple.compare_per_checklist(x1, x2)
  if rel.relationship == LT:
    return u1                   # Prefer more tipward
  if rel.relationship == GT:
    return u2                   # Prefer more tipward
  # The following is probably not necessary, but it's tidy?
  # if necessary, should work harder to canonicalize, yes?
  parts1 = get_parts(x1)
  parts2 = get_parts(x2)
  if is_accepted(x2):
    if False and (not parts1.protonymp and parts2.protonymp and parts2.genus != None):
      log("# preferring protonym %s to %s" % (blurb(u2), blurb(u1)))
      return u2                   # Prefer protonym - NOT
    elif not is_accepted(x1):
      return u2
    elif simple.descends_from(x2, x1):
      return u2
  return u1

# Only workspace nodes have uf records

def get_exemplar_uf(u):
  #log("# Thinking about exemplar for %s" % blurb(u))
  probe = really_get_exemplar_uf(u, None)
  if probe: return probe
  AB = get_workspace(u)
  exem = [None, u, None] if isinA(AB, u) else [None, None, u]
  uf = UnionFindable(exem)
  assert exem[1] or exem[2]
  set_exemplar_uf(u, uf)
  return uf

# Union-find nodes start out with no xid, then get xid when exemplars are identified

(really_get_exemplar_uf, set_exemplar_uf) = \
  prop.get_set(prop.declare_property("exemplar_uf"))


# Returns exemplar record (xid, u, v) or None.

def get_exemplar(z):
  uf = really_get_exemplar_uf(z, None)
  if uf:
    uf = uf.find()
    r = uf.payload()
    (xid, u, v) = r
    if u and v:
      if xid == None:
        # Create exemplar id (for set operations) on demand
        ws = get_workspace(u)
        xid = fresh_exemplar_id(ws)
        r[0] = xid
        ws.exemplar_ufs[xid] = uf
        #log("# Exemplar %s: (%s) <-> (%s)" % (xid, blurb(u), blurb(v)))
      return r
  return None

def get_bare_exemplar(z):
  uf = really_get_exemplar_uf(z, None)
  if uf:
    r = uf.payload()
    (xid, u, v) = r
    if u and v:
      return r
  return None
  

# -----------------------------------------------------------------------------
# Find tipward record matches (TRMs)
# Not used at present

(get_tipward, set_tipward) = prop.get_set(prop.declare_property("tipward"))

# Postorder iteration over one of the two summands.
# Sets the 'tipward' property meaning we want to use the node to represent its
# exemplar.

def analyze_tipwards(AB):
  def traverse(AB, x):
    seen = seen2 = 0
    for c in get_children(x, ()):
      seen += traverse(AB, c)
    for c in get_synonyms(x, ()):
      seen2 += traverse(AB, c)
    if seen == 0:
      # not seen means that this node could be tipward
      u = AB.in_left(x)            # u is tipward...
      v = get_link(u, None)        # same name -> same exemplar(s?)
      if monitor(v):
        log("# exemplar: Tipward link %s -> %s" % (blurb(u), blurb(v)))
      if v != None:                     # not None, not False = ambiguous
        if v != False:
          assert separated(u, v)
          seen = 1
        # u is tipward and has link, but doesn't necessarily link tipward ...
        set_tipward(u, v)
        if monitor(u):
          log("# exemplar: Set tipward %s -> %s" % (blurb(u), blurb(v)))
    return seen + seen2
  traverse(AB, AB.A.top)
  traverse(swap(AB), AB.B.top)
