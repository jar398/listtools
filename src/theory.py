#!/usr/bin/env python3

import types, functools
import property as prop, checklist, workspace

from util import log
from checklist import *
from workspace import *
from match_records import match_records
from rcc5 import rcc5_relationship, reverse_relationship

def isinA(AB, z):
  return AB.case(z, lambda x: True, lambda x: False)

def isinB(AB, z):
  return AB.case(z, lambda x: False, lambda x: True)

#-----------------------------------------------------------------------------

# Store equivalents on nodes of AB

def find_equivalents(AB):
  v = AB.in_left(AB.A.top)
  w = AB.in_right(AB.B.top)
  set_equivalent(v, relation(EQ, w, "top"))
  set_equivalent(w, relation(EQ, v, "top"))
  count = 1
  for y in checklist.postorder_records(AB.B):
    w = AB.in_right(y)
    rel2 = get_matched(w)
    if rel2:
      v = rel2.record
      b1 = get_block(v, BOTTOM_BLOCK)
      b2 = get_block(w, BOTTOM_BLOCK)
      if same_block(b1, b2):    # Perhaps empty
        rel1 = get_matched(v)
        if rel1:
          assert rel1.record == w
          set_equivalent(v, relation(EQ, w, rel1.status, rel1.note))
          set_equivalent(w, relation(EQ, v, rel2.status, rel2.note))
          count += 1
    else:
      v = get_cross_mrca(w, BOTTOM)
      if v != BOTTOM:
        u = get_cross_mrca(v, BOTTOM)
        if u == w:
          log("by subtension %s = %s" % (blurb(v), blurb(w)))
          set_equivalent(v, relation(EQ, w, None, "by subtension"))
          set_equivalent(w, relation(EQ, v, None, "by subtension"))
          count += 1
  log("# found %s B nodes with A equivalents" % count)

(get_equivalent, set_equivalent) = prop.get_set(prop.declare_property("equivalent"))

#-----------------------------------------------------------------------------

# This is the implementation of the RCC-5 theory of AB (A+B).

# v and w are Records in ... AB, not in A or B
# Returns a Relative to w

# TBD: set status to "accepted" or "synonym" as appropriate

def cross_relationship(AB, v, w):
  assert isinA(AB, v) != isinA(AB, w)

  e = get_equivalent(v, None)
  if e and e.record == w: return e

  b1 = get_block(v, BOTTOM_BLOCK)
  b2 = get_block(w, BOTTOM_BLOCK)
  ship = block_relationship(b1, b2)
  if ship != EQ:
    return relation(ship, w, None, "MTRM set comparison")

  if True:

    # Look for a record match for v or w, and punt to simple case
    rel1 = rel2 = None
    n = get_matched(v)              # v, n, n.record in AB (from B)
    if n:
      # NO - the relative.record needs to be w
      rel1 = simple_relationship(get_outject(n.record), get_outject(w))

    m = get_matched(w)      # w, m, m.record in AB
    if m:
      rel2 = simple_relationship(get_outject(v), get_outject(m.record))

    # There's valuable information in the note... keep it?

    rel = rel1 or rel2
    if rel:
      if (not rel1 or not rel2 or
          rel1.relationship == rel2.relationship):
        return relation(rel.relationship, w, rel.status, rel.note)
      return relation(NOINFO, w, None, "lineage out of order??")

  if is_empty_block(b1):
    return relation(DISJOINT, w, None, "no overlap")

  # No equivalent either way, but same block, so OVERLAP.
  # TBD: If no child or parent of v or w is in same block, then EQ.
  return relation(OVERLAP, w, None, "same block but no record match")

# RCC-5 relationship across the two checklists
# x and y are in AB

def cross_lt(AB, v, w):
  return cross_relationship(AB, v, w).relationship == LT

# -----------------------------------------------------------------------------

# Returns a Relative y to x such that x = y.  x and y are in AB, not A or B.
# TBD: should cache this.

def find_equivalent(AB, v):
  y = None
  rel1 = get_matched(v)            # Record match
  if rel1:
    y = rel1.record
  else:
    w = AB.get_cross_mrca(v)
  rel2 = cross_relationship(AB, v, w)
  if rel2.relationship == EQ:
    return rel2                 # hmm. combine rel1 and rel2 ?
  return None

# -----------------------------------------------------------------------------
# Least node in opposite (B) checklist that's definitely greater than w (in A).

def cross_superior(AB, w):
  if isinB(AB, w):
    AB = AB.swap()
  if monitor(w): log("cross_superior %s" % blurb(w))
  w0 = w
  if is_empty_block(get_block(w, BOTTOM_BLOCK)):
    # increase w until it overlaps with B
    while True:
      x = get_outject(w)
      if is_top(x): return None
      rel = get_superior(x, None)
      assert rel
      ra = rel.record
      w = AB.in_left(ra)        # in A
      if is_empty_block(get_block(w, BOTTOM_BLOCK)):
        # w0 < w < ... overlaps with B ...
        if monitor(w0): log("cross_superior None empty %s" % blurb(w))
      else:
        # w0 < w = overlaps with B 
        break
    q = get_cross_mrca(w, BOTTOM)
    if monitor(w0): log("cross_superior q %s" % blurb(q))
    return q # ???
  else:
    # increase q until x is less than it
    q = get_cross_mrca(w, BOTTOM)     # candidate in AB
    assert isinB(AB, q)
    if monitor(w0): log("xsup 2 %s %s" % (blurb(x), blurb(q),))
    while True:
      ship = cross_relationship(AB, w, q).relationship
      if ship == LT:
        # Satisfactory parent candidate.
        return q
      qb = get_outject(q)
      if monitor(w0): log("xsup 3 %s %s" % (blurb(q), blurb(qb)))
      rel = get_superior(qb, None)
      if rel:
        sb = rel.record
        assert sb, blurb(sb)
        q = AB.in_right(sb)
        if ship == EQ:
          return q
      else:
        if monitor(w): log("cross_superior None 2 %s" % blurb(x))
        return None             # Off top

# -----------------------------------------------------------------------------
# Find mutual tipward record matches (MTRMs)

(get_tipward, set_tipward) = prop.get_set(prop.declare_property("tipward"))

def analyze_tipwards(AB):
  find_tipwards(AB.A, AB.in_left)
  find_tipwards(AB.B, AB.in_right)

# Postorder iteration over one of the two summands.

def find_tipwards(A, in_left):
  def traverse(x):
    seen = False
    for c in get_inferiors(x):
      seen = traverse(c) or seen
    if not seen and not is_top(x):
      z = in_left(x)            # z is tipward...
      rel = get_matched(z)  # rel is a Relative...
      if monitor(x): log("tipwards %s %s" % (blurb(x), blurb(rel)))
      if rel:
        # z is tipward and has a match, not necessarily tipward ...
        set_tipward(z, rel)
        seen = x
    return seen
  traverse(A.top)

# argument is in AB, result is Relative in AB (from opposite summands)
# returns Relative - a record match to the argument

def get_mtrm(z):
  rel = get_tipward(z, None)
  #if monitor(z): log("get_mtrm 1 %s %s" % (blurb(z), blurb(rel)))
  if rel:
    q = rel.record              # in AB
    assert q
    rel2 = get_tipward(q, None)
    #if monitor(z): log("get_mtrm 2 %s %s" % (blurb(q), blurb(rel2)))
    if rel2 and rel2.record == z:
      return rel
    else:
      if monitor(z):
        log("not mutual trm %s %s" % (blurb(q), blurb(rel2)))
  return None

# -----------------------------------------------------------------------------
# cross_mrca (in AB) satisfies: 
#   1. if y = cross_mrca(x) then x ~<= y
#   2. if furthermore x' = cross_mrca(y), then
#      (a) if x' <= x then x ≅ y, 
#      (b) if x' > x then x > y.  [why?]
#
# i.e.   x ≲ y ≲ x' ≤ x   →   x ≅ y   (same 'block')
# i.e.   x ≲ y ≲ x' > x   →   x < y
#
# Cached in AB nodes under the 'cross_mrca' property.
# Needed for equivalent and cosuperior calculations

def compute_cross_mrcas(AB):
  def crosses(AB):
    def traverse(x):            # arg in A, result in B
      m = BOTTOM                # identity for mrca
      z = AB.in_left(x)         # in AB
      probe = get_mtrm(z)       # in AB
      if probe:
        m = AB.case(probe.record,   # in B
                    lambda x:2/(3-3),
                    lambda y:y)
      else:
        for c in get_inferiors(x):
          q = traverse(c)       # in B
          if q:
            m = mrca(q, m)      # in B
      if m != BOTTOM:
        set_cross_mrca(AB.in_left(x), AB.in_right(m))
        return m
      return None
    traverse(AB.A.top)
  crosses(AB)
  crosses(swap(AB))

def swap(AB):
  BA = AB.swap()
  BA.A = AB.B
  BA.B = AB.A
  return BA

(get_cross_mrca, set_cross_mrca) = \
  prop.get_set(prop.declare_property("cross_mrca"))

# -----------------------------------------------------------------------------
# Precompute 'blocks' (tipe sets, implemented in one of various ways).
# A block is represented as either a set or TOP_BLOCK.
# Blocks are stored on nodes in AB.

def compute_blocks(AB):
  def traverse(x, in_left):
    # x is in A or B
    if monitor(x): log("computing block for %s" % (blurb(x),))
    v = in_left(x)
    e = BOTTOM_BLOCK
    for c in get_inferiors(x):  # inferiors in A/B
      e = combine_blocks(e, traverse(c, in_left))
      if monitor(c): log("got subblock %s -> %s" % (blurb(c), len(e),))
    e2 = possible_mtrm_block(AB, v)
    if (not is_empty_block(e) and not is_empty_block(e2)
        and monitor(x)):
      log("# non-mutual TRM: %s" % blurb(x))
    e = combine_blocks(e2, e)
    if monitor(x): log("got mtrm %s -> %s" % (len(e2), len(e),))
    if is_top(x):
      e = TOP_BLOCK
    if e != BOTTOM_BLOCK:
      set_block(v, e)
    #log("# block(%s) = %s" % (blurb(x), e))
    return e
  traverse(AB.A.top, AB.in_left)
  traverse(AB.B.top, AB.in_right)

# z is in AB

def possible_mtrm_block(AB, z):
  rel = get_mtrm(z)             # a Relative
  return mtrm_block(get_outject(z), get_outject(rel.record)) if rel else BOTTOM_BLOCK

# -----------------------------------------------------------------------------
# Implementation of blocks as Python sets of 'tipes.'

# RCC-5 relationship between two blocks

def block_relationship(e1, e2):   # can assume overlap
  if e1 == e2: return EQ          # same block
  elif e1 == TOP_BLOCK: return GT
  elif e2 == TOP_BLOCK: return LT
  elif e1.issubset(e2): return LT
  elif e2.issubset(e1): return GT
  elif e1.isdisjoint(e2): return DISJOINT
  else: return CONFLICT

def same_block(e1, e2):
  return e1 == e2

# Lattice join (union) of two blocks

BOTTOM_BLOCK = set()
TOP_BLOCK = True
def combine_blocks(e1, e2):
  if e1 == BOTTOM_BLOCK: return e2
  if e2 == BOTTOM_BLOCK: return e1
  if e1 == TOP_BLOCK: return e1
  if e2 == TOP_BLOCK: return e2
  return e1 | e2
def same_block(e1, e2): return e1 == e2
def is_empty_block(e): return e == BOTTOM_BLOCK

def mtrm_block(x, y):           # x = y
  return {min(x.id, y.id)}

(get_block, set_block) = prop.get_set(prop.declare_property("block"))

# -----------------------------------------------------------------------------
# Same-tree relations

#    x1 ? y1     'larger'
#   /       \
#  x         y   'smaller'

# Returns a Relative explaining justification

def simple_relationship(x, y):             # Within a single tree
  (x, synx) = nip_synonym(x)
  (y, syny) = nip_synonym(y)
  if x == y:
    if synx or syny:
      if synx and syny:
        return relation(NOINFO, y, None, "sibling synonyms") # blah
      elif synx:
        return relation(LE, y, "synonym", "synonym <= accepted")
      else:
        return relation(GE, y, None, "accepted >= synonym")
    else:
      return relation(EQ, y, None, "accepted = accepted")
  (x1, y1) = find_peers(x, y)    # Decrease levels as needed
  if x1 == y1:     # x <= x1 = y1 >= y
    if x1 == x:     # x = x1 = y1 > y, x > y
      return relation(GT, y, None, "in-checklist")
    elif y1 == y:
      return relation(LT, y, None, "in-checklist")
    else:
      assert False  # can't happen
  else:
    return relation(DISJOINT, y, "in-checklist")

# Assumes already 'nipped'

def find_peers(x, y):
  if get_level(x, None) == None or get_level(x, None) == None:
    clog("# No level for one of these:", x, y)
    return NOINFO
  i = get_level(x)
  while i < get_level(y):
    y = get_parent(y)
  j = get_level(y)
  while get_level(x) > j:
    x = get_parent(x)
  return (x, y)

# MRCA within the same tree

def mrca(x, y):
  if x == BOTTOM: return y
  if y == BOTTOM: return x
  if x == y: return x

  (x, _) = nip_synonym(x)
  (y, _) = nip_synonym(y)
  (x, y) = find_peers(x, y)
  while not (x is y):
    x = get_parent(x)
    y = get_parent(y)
  return x

BOTTOM = None

def simple_le(x, y):
  # Is x <= y?  Scan from x upwards, see if we find y
  if x == y: return True
  (y, _) = nip_synonym(y)
  x1 = x
  stop = get_level(y)

  while get_level(x1) > stop:
    x1 = get_parent(x1)    # y1 > previously, level(y1) < previously

  return x1 == y

def simple_lt(x, y):
  return simple_le(x, y) and x != y

def in_same_tree(AB, x, y):
  return (AB.case(x, lambda x:1, lambda x:2) ==
          AB.case(y, lambda x:1, lambda x:2))

(really_get_level, set_level) = prop.get_set(prop.declare_property("level", inherit=False))

# ----- Levels maintenance

# 2nd result is True if x is a synonym, False otherwise
#
def nip_synonym(x):
  synx = False
  while True:
    sup = get_superior(x, None)
    if not sup or sup.relationship == ACCEPTED:
      break
    if sup.relationship != EQ:
      synx = True
    x = sup.record
  assert get_level(x, None) != None, blurb(x)
  return (x, synx)

def get_level(x, default=None):
  i = really_get_level(x, None)
  if i == None:
    set_level(x, "cycle")
    p = get_parent(x)
    if not p:
      i = 1                 # shouldn't happen.  top is level 1
    else:
      assert p != x
      lev = get_level(p, None)
      assert lev != "cycle", \
        ("**** cycle detected while traversing ancestor chain: %s->%s" %
            (blurb(x), blurb(p)))
      i = lev + 1
    set_level(x, i)
    #if monitor(x):
    #  log("# level(%s) = %s" % (blurb(x), i))
  return i

def get_parent(x):
  rp = get_superior(x, None)
  if rp:
    if False and rp.relationship != ACCEPTED:    ##????
      return get_parent(rp.record)
    else:
      return rp.record
  else: return None

def ensure_levels(S):
  def cache(x, n):
    set_level(x, n)
    for c in get_inferiors(x):
      # This isn't right -- get_level works better
      cache(c, n+1)
  cache(S.top, 1)

# -----------------------------------------------------------------------------
# load_matches is in file checklist.py

# Load/dump a set of provisional matches (could be either record match
# or taxonomic matches... but basically, record matches).  The matches are stored 
# as Relations under the 'match' property of nodes in AB.

# x and y are in AB

def load_matches_redundant_definition(row_iterator, AB):
  log("# Loading matches")

  header = next(row_iterator)
  plan = prop.make_plan_from_header(header)
  match_count = 0
  miss_count = 0
  for row in row_iterator:
    # row = [matchID (x id), rel, taxonID (y id), remark]
    match = prop.construct(plan, row)
    x = y = None
    xkey = get_match_key(match, None)
    if xkey:
      x_in_A = look_up_record(AB.A, xkey)
      if x_in_A:
        x = AB.in_left(x_in_A)
    ykey = get_primary_key(match, None)
    if ykey:
      y_in_B = look_up_record(AB.B, ykey)
      if y_in_B:
        y = AB.in_right(y_in_B) 

    ship = rcc5_relationship(get_match_relationship(match)) # EQ, NOINFO
    note = get_basis_of_match(match, MISSING)
    # x or y might be None with ship=NOINFO ... hope this is OK
    if y:
      set_match(y, relation(reverse_relationship(ship), x, None,
                            reverse_note(note)))
      match_count += 1
    if x:
      set_match(x, relation(ship, y, None, note))
      match_count += 1
    if x and y: match_count += 1
    else: miss_count += 1

  log("# loaded %s matches" % (match_count))

def reverse_note(note):
  if ' ' in note:
    return "↔ " + note            # tbd: deal with 'coambiguous'
  else:
    return note

# x should be in A or B ... ?

def is_accepted(x):
  sup = get_superior(x, None)
  return (not sup) or sup.relationship != SYNONYM
