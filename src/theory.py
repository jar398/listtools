#!/usr/bin/env python3

# Prepare for doing within-tree MRCA operations
#  cache_levels(C, roots)

import sys, csv, argparse
import functools
import util
import property as prop
import checklist

from rcc5 import *



# Reading and writing articulation sets, and doing inference on them

# Proposing alignments that create spanning trees (maximally
# consistent, ><-free theories)

# Read a theory from 4 files (A, B, seed, record matches)

def load_theory(a_iterator, b_iterator,
                seed_iterator=None, rm_sum_iterator=None):
  if rm_sum_iterator == None:
    # Need to copy the iterators!
    (a_iterator, a_iterator_copy) = dup_iterator(a_iterator)
    (b_iterator, b_iterator_copy) = dup_iterator(b_iterator)
    rm_sum_iterator = match_records.match_records(a_iterator_copy, b_iterator_copy)
  A = checklist.load_source(a_iterator, "A")
  B = checklist.load_source(b_iterator, "B")
  # Load record matches
  rm_sum = load_sum(rm_sum_iterator, A, B)
  # Create 'theory object'
  AB = checklist.make_sum(A, B, "", rm_sum=rm_sum)
  init_source(A)
  init_source(B)


def init_source(A):
  (A.get_level, A.set_level) = \
    prop.get_set(prop.get_property("level"), context=A)


def dup_iterator(iter):
  clunk = list(iter)
  return ((x for x in clunk), (x for x in clunk))


# Read/write an articulation set

# Finally need to sort this out...

"""
a way to identify a record in either source;
hmmmm...

need a nice way to get unique names for records.
this will be hard, but can use .... ?

canonical / scientific / dbid ...

hmm.  maybe do it at the taxon id level?  need tooling to create these
things, but much more reliable and concise that way.

so, taxonID_A and taxonID_B columns as before.

Oh, the merged taxonomy is a very different artifact!!!

So, at least 2 different outputs.

SUGGESTED:
  links inherited from sources  ?  goes without saying
  every = articulation  ('peers')
  consequences of =.  equations between children
  child/parent links in monotypic ladders, for completeness
    ('altParent', 'altAccepted')
  dropped children, to avoid conflict and change >< to > or < or =
    'splice out'  or 'synonymize'

INFERRED:
  'estimates'  (individual or as [a, z] estimates)
  child/parent and synonym/accepted links, many of them

So... I think this is not too bad.  An alignment has
articulations
  - at most one = per node (I think)
  - at most parent or accepted link per node
  - either can override ... how??
  - also, deleted or deprecated conflicting nodes
    (conversion to synonym)
  - mutex combination

OK a list with pairs of TNUs and a single relation seems about right.
But the relationship will always be = or < (or <= for synonyms).
a = b
a < c
x <= c
with reasons:
  tipward match on canonicalName,
  similar extension and match on canonicalName,
  coextensional,
  conflict resolution,

"""

# -----------------------------------------------------------------------------
# RCC5 decision procedure
#   1. for nodes in same tree
#   2. for nodes in different trees

# This seems wrong somehow
# Things would be simpler if only source nodes had levels

def in_same_tree(C, x, y):
  return source_checklist(x) == source_checklist(y)

def le(C, x, y):
  return tree_le(C, x, y) if in_same_tree(C, x, y) else cross_le(C, x, y)

def lt(C, x, y):
  return lt(C, x, y) and not eq(C, x, y)

def eq(C, x, y):
  return tree_eq(C, x, y) if in_same_tree(C, x, y) else cross_le(C, x, y)

# -----------------------------------------------------------------------------
# Decide about a single source tree

# E.g. level(x) >= level(y)  should imply  x <= y  if not disjoint

def tree_eq(C, x, y):
  return x == y

def tree_le(C, x, y):
  y1 = x     # proceeds down to y
  # if {level(x) >= level(y1) >= level(y)}, then x <= y1 <= y or disjoint
  stop = get_level(y)
  while get_level(y1) > stop:
    y1 = get_superior(y1)    # y1 > previously, level(y1) < previously
  # level(y) = level(y1) <= level(x)
  return y1 == y    # Not > and not disjoint

def tree_lt(x, y): return tree_le(x, y) and x != y

def find_peers(C, x, y):
  while C.get_level(x) < C.get_level(y):
    y = get_superior(y)
  while C.get_level(x) > C.get_level(y):
    x = get_superior(x)
  return (x, y)

# MRCA within the same tree

def mrca(x, y):
  if x == BOTTOM: return y
  if y == BOTTOM: return x
  (x, y) = find_peers(x, y)
  while not (x is y):
    x = get_superior(x)
    y = get_superior(y)
  return x

def relation(x, y):             # Within a single tree
  (x1, y1) = find_peers(x, y)    # Decrease levels as needed
  assert not x1.checklist is y1.checklist
  if x1 == y1:
    if x == y:
      return EQ
    elif x1 == x:
      return LT     # y > y1 = x
    elif y1 == y:
      return GT     # x > x1 = y
    else:
      assert False
  else:
    x2 = get_accepted(x1, x1)
    y2 = get_accepted(y1, x1)
    if x2 == y2:                                  # Same accepted
      if x1 != x2 and y1 != y2: return NOINFO    # synonym ? synonym
      if x1 != x2: return LE                      # synonym <= accepted
      else:        return GE                      # accepted >= synonym
    else:
      return DISJOINT

# -----------------------------------------------------------------------------
# Decide about relations between two trees

# Returns a pair (overlap, outlier)
#   overlap is in y and x {< = > ><} y (but not !)
#   outlier is in y and x {< >< !} y (but not {> =})
#   either can be None

def analyze(AB, x, y):
  (over, out) = (None, None)
  for d in AB.get_children(y):
    (over2, out2) = analyze(x, d)
    over = over or over2; out = out or out2
    if over and out: return (over, out)
  for d in AB.get_synonyms(y):
    (over2, out2) = analyze(x, d)
    over = over or over2; out = out or out2
    if over and out: return (over, out)
  if over or out:      # Exclude peripherals
    return (over, out)
  m = AB.record_match(y)
  j = AB.join(m, y)
  return (j, None) if m else (None, j)

def outlier(AB, x, y):
  (over, _) = analyze(x, y)
  return out

# Cache this.

def get_shadow(AB, x, other):
  m = functools.reduce(mrca,
                       (get_shadow(c) for c in get_inferiors(x)),
                       None)
  if m == None:
    return record_match(x, other) or None
  return m

# RCC5 decision procedure, where x and y are in different sources.

def cross_ge(AB, x, y):
  return not outlier(AB, x, y)

def cross_eq(AB, x, y):
  # see if they're peers ??
  return cross_ge(x, y) and cross_ge(y, x)

def gt(AB, x, y):
  return cross_ge(x, y) and not cross_ge(y, x)

# -----------------------------------------------------------------------------
# Cache every node's level (distance to root)
#   simple recursive descent from roots

# Level is contravariant with RCC5: x < y implies level(x) > level(y)

def cache_levels(S, roots):
  def cache(x, n):
    S.set_level(x, n)
    for c in S.get_inferiors(x):
      cache(c, n+1)
  for root in roots:
    cache(root, 1)



if __name__ == '__main__':
  src = load_source(csv.reader(sys.stdin), "A")

  props = (prop.get_property(label)
           for label in ("taxonID", "canonicalName", "parentNameUsageID"))
  writer = csv.writer(sys.stdout)
  for row in emit_rows(src, props): writer.writerow(row)

