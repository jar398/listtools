#!/usr/bin/env python3

import math
import property as prop
import checklist, workspace, simple
#import exemplar

from util import log, UnionFindable
from checklist import *
from workspace import *
from linkage import get_link, really_find_links
from linkage import pick_better_record
from typify import equate_typifications

# Exemplars (representing individual specimens or occurrences with
# known classification in BOTH checklists) are chosen heuristically
# via name matches.  Not every name match qualifies; only those that
# are 'tipward' do.  'Tipward' means that at least one of the two taxa
# has no tipward accepted descendant taxon.

# See AB.exemplar_ufs for map {xid: uf, ...}
# Each exemplar has a representative in A and one in B.

# find_exemplars is defined in exemplar.py.  Does 2 passes.
# find_some_exemplars is just one pass; invoked twice.

def find_some_exemplars(AB, subproblems, get_pre_estimate):
  log("# Finding some exemplars")
  really_find_links(AB, subproblems, get_pre_estimate)

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
            equate_typifications(u, v)

    traverse(CD.A.top)

  do_exemplars(AB)
  do_exemplars(swap(AB))
  equate_typifications(AB.in_left(AB.A.top),
                       AB.in_right(AB.B.top))

# Issues: 
#   1. choice of taxa (both sides)
#   2. choice of name
  

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
