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

# See AB.exemplars for map {id: [id, v, w], ...}
# Each exemplar has a representative in A and one in B.

# find_exemplars is defined in exemplar.py.  Does 2 passes.
# find_some_exemplars is just one pass; invoked twice.

CARE_TIPWARD = False

def find_some_exemplars(AB, get_pre_estimate):
  linkage.really_find_links(AB, get_pre_estimate)
  if CARE_TIPWARD:
    analyze_tipwards(AB)

  def do_exemplars(CD):
    def traverse(x, species):      # species is an ancestor of x, in A
      u = CD.in_left(x)
      
      if not CARE_TIPWARD:
        # All reciprocal links become exemplars
        v = get_link(u, None)
        if v:
          u2 = get_link(v, None)
          if u2 == u:
            equate_exemplars(u, v)
        for c in get_inferiors(x):
          traverse(c, species)

      else:
        # Only reciprocal tipward links, and species, become exemplars

        v = get_tipward(u, None)
        # u has priority over v
        if v:                     # not None or False
          # 1. Mutual tipward matches mean we have an exemplar
          u_back = get_tipward(v, None)
          if u_back is u:         # Mutual?
            if monitor(u) or monitor(v):
              log("# exemplar: Mutual: %s -> %s" % (blurb(u), blurb(v)))
            equate_exemplars(u, v) # Normal case.
          elif u_back == None:
            if monitor(u) or monitor(v):
              log("# exemplar: Nonmutual: %s -> %s -> (none)" % 
                  (blurb(u), blurb(v)))
            equate_exemplars(u, v)
          elif u_back == False:
            if monitor(u) or monitor(v):
              log("# exemplar: Ambiguous: %s -> %s -> (2 or more)" % 
                  (blurb(u), blurb(v)))
            equate_exemplars(u, v)
          else: # u_back != None:     doesn't seem to happen?
            log("# exemplar: Returns to different tipward %s -> %s -> %s" % 
                (blurb(u), blurb(v), blurb(u_back)))
        else:
          if not species and get_rank(u, None) == 'species':
            species = u

          found_species = False
          for c in get_inferiors(x):
            found_species |= traverse(c, species)

          # 2. A matched species with descendents, having an epithet different
          # from any descendent, is also an exemplar. (but matched to what??)


          v = get_link(u, None)   # None, False, node

          # WORK IN PROGRESSS

          species = u
          if not found_species and species and linkage.homotypic(u, species): # same epithet stem
            equate_exemplars(u, species)    # ? what is this for ?
            found_species = True

          # TBD: Handle multiple subspecifics with same epithet...

          return found_species

    traverse(CD.A.top, None)

  do_exemplars(AB)
  do_exemplars(swap(AB))
  equate_exemplars(AB.in_left(AB.A.top),
                   AB.in_right(AB.B.top))
  report_on_exemplars(AB)

def report_on_exemplars(AB):
  exemplars = set()             # for counting distinct exemplars
  for x in preorder_records(AB.A):
    u = AB.in_left(x)
    ex = get_exemplar(u)
    if ex:
      exemplars.add(ex[0])      # id only
  log("-- Exemplars: %s " % len(exemplars))

# Issues: 
#   1. choice of taxa (both sides)
#   2. choice of name

def equate_exemplars(u, v):     # opposite checklists. u might be species
  if u != v:
    # assert u descends from v if in same checklist?
    uf = get_exemplar_uf(u)
    vf = get_exemplar_uf(v)
    merge_representatives(uf.find(), vf.find())
    uf.absorb(vf)

def merge_representatives(uf, vf): # absorb into uf
  r = uf.payload()
  (i1, u1, v1) = r
  (i2, u2, v2) = vf.payload()
  r[1] = unify_records(u1, u2)
  r[2] = unify_records(v1, v2)
  if r[1] and r[2]:
    #log("# Same exemplar: e(%s) = e(%s)" % (blurb(r[1]), blurb(r[2])))
    pass
  if i1 and i2:
    r[0] = min(i1, i2)          # wait, what is this for?
  else:
    r[0] = i1 or i2

# x1 and u2 are in same checklist (not in workspace)

def unify_records(x1, x2):
  if not x2: return x1
  if not x1: return x2
  rel = simple.compare_per_checklist(x1, x2)
  if rel.relationship == LT:
    return x1                   # Prefer more tipward
  if rel.relationship == GT:
    return x2                   # Prefer more tipward
  # The following is probably not necessary, but is tidy?
  # if necessary, should work harder to canonicalize, yes?
  parts1 = get_parts(x1)
  parts2 = get_parts(x2)
  if (not parts1.protonymp and parts2.protonymp and parts2.genus != None):
    log("# preferring protonym %s to %s" % (blurb(x2), blurb(x1)))
    return x2                   # Prefer protonym
  else:
    return x1

# Only workspace nodes have uf records

def get_exemplar_uf(u):
  #log("# Thinking about exemplar for %s" % blurb(u))
  probe = really_get_exemplar_uf(u, None)
  if probe: return probe
  AB = get_workspace(u)
  uf = UnionFindable([None, get_outject(u), None] if isinA(AB, u) else [None, None, get_outject(u)])
  set_exemplar_uf(u, uf)
  return uf

# Union-find nodes start out with no id, then get id when exemplars are identified

def init_exemplar(AB, exemplars, id, x, y):
  assert get_source(x) is AB.A
  assert get_source(y) is AB.B
  ex = (id, x, y)
  uf = UnionFindable(ex)
  set_exemplar_uf(u, uf)
  set_exemplar_uf(v, uf)
  assert not id in exemplars, id
  exemplars[id] = ex
  return ex

(really_get_exemplar_uf, set_exemplar_uf) = \
  prop.get_set(prop.declare_property("exemplar_uf"))


# Returns exemplar record or None.  TO BE USED ONLY AFTER
# analyze_exemplars HAS FINISHED ITS WORK.

def get_exemplar(u):
  ws = get_workspace(u)
  uf = really_get_exemplar_uf(u, None)
  if uf:
    r = uf.payload()
    (id, x, y) = r
    if x and y:
      if not id:
        # Create id (for set operations) on demand
        id = fresh_exemplar_id(ws)
        r[0] = id
        ws.exemplars[id] = uf
        #log("# Exemplar %s: e(%s) = (%s)" % (id, blurb(x), blurb(y)))
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
