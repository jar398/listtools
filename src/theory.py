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

def find_equivalents(AB):
  for y in checklist.postorder_records(AB.B):
    rel = get_match(y, None)
    if rel and rec.record:
      x = rel.record
      b1 = get_block(x, BOTTOM_BLOCK)
      b2 = get_block(y, BOTTOM_BLOCK)
      if same_block(b1, b2):    # Perhaps empty
        rel2 = get_match(x, None)
        if rel2 and rel2.record:
          set_equivalent(y, rel)
          set_equivalent(x, rel2)

(get_equivalent, set_equivalent) = prop.get_set(prop.declare_property("equivalent"))

#-----------------------------------------------------------------------------

# This is the implementation of the RCC-5 theory of AB (A+B).

# x and y are Records in ... AB, not in A or B
# Returns a Relative to y

# TBD: set status to "accepted" or "synonym" as appropriate

def cross_relationship(x, y):
  e = get_equivalent(x, None)
  if e and e.record == y: return e

  b1 = get_block(x, BOTTOM_BLOCK)
  b2 = get_block(y, BOTTOM_BLOCK)
  ship = block_relationship(b1, b2)
  if ship != EQ:
    return relation(ship, y, None, "MTRM set comparison")

  if False:

    # Look for a record match for x or y, and punt to simple case
    rel1 = rel2 = None
    n = get_match(x, None)              # n.record is in AB (from B)
    if n and n.record and n.relationship == EQ:
      # NO - the relative.record needs to be y
      rel1 = simple_relationship(outject(n.record), outject(y))

    m = get_match(y, None)
    if m and m.record and m.relationship == EQ:
      rel2 = simple_relationship(outject(x), outject(m.record))

    rel = rel1 or rel2
    if rel:
      if (not rel1 or not rel2 or
          rel1.relationship == rel2.relationship):
        return relation(rel.relationship, y, rel.status, rel.note)
      return relation(NOINFO, y, None, "lineage out of order??")

  if is_empty_block(b1):
    return relation(DISJOINT, y, None, "no overlap")

  # No equivalent either way, but same block, so OVERLAP.
  # TBD: If no child or parent of x or y is in same block, then EQ.
  return relation(OVERLAP, y, None, "same block but no record match")

# RCC-5 relationship across the two checklists
# x is in A and y is in B, or vice versa

def cross_lt(x, y):
  return cross_relationship(x, y).relationship == LT

# -----------------------------------------------------------------------------

# Returns a Relative y to x such that x = y.  x and y are in AB, not A or B.
# TBD: should cache this.

def find_equivalent(AB, x):
  y = None
  rel1 = get_match(x)            # Record match
  if rel1 and rel1.record and rel1.relationship == EQ:
    y = rel1.record
  else:
    y = AB.get_cross_mrca(x)
  rel2 = cross_relationship(x, y)
  if rel2.relationship == EQ:
    return rel2                 # hmm. combine rel1 and rel2 ?
  return None

def get_injected_superior(AB, a):
  return AB.case(z,
                 lambda x: AB.in_left(get_superior(x)),
                 lambda y: AB.in_right(get_superior(y)))

# -----------------------------------------------------------------------------
# Least node in opposite checklist that's definitely greater than x.
# x is in A (or B), not AB, and result is in B (or A), not in AB.
# Returns a Relative to y in B, or None

def cross_superior(AB, x):
  q = AB.get_cross_mrca(x, None) # in B
  # increase q until x is less than it
  i = 0
  while q and not cross_lt(x, q):
    rel = get_superior(q, None)
    if rel:
      q = rel.record
      assert len(checklist.get_source_name(q)) == 1
    else: break
  return q

# -----------------------------------------------------------------------------
# Find mutual tipward record matches (MTRMs)

(get_tipward, set_tipward) = prop.get_set(prop.declare_property("tipward"))

def analyze_tipwards(AB):
  find_tipwards(AB.A, AB.in_left)
  find_tipwards(AB.B, AB.in_right)

# Postorder iteration over one of the two summands.

def find_tipwards(A, in_left):
  ensure_inferiors_indexed(A)  # children and synonyms properties
  def traverse(x):
    seen = False
    for c in get_inferiors(x):
      seen = traverse(c) or seen
    if not seen and not is_top(x):
      z = in_left(x)            # z is tipward...
      rel = get_match(z, None)  # rel is a Relative...
      if rel:
        # z is tipward and has a match, not necessarily tipward ...
        set_tipward(z, rel)
        seen = x
    return seen
  traverse(A.top)

# argument is in AB, result is in AB (from opposite summands)
# returns Relative - a record match to the argument

def get_mtrm(z):
  rel = get_tipward(z, None)
  if rel:
    q = rel.record              # in AB
    rel2 = get_tipward(q, None)
    if rel2 and rel2.record == z:
      return rel
  return None

# -----------------------------------------------------------------------------
# cross_mrca satisfies: 
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
  (get_cross_mrca, set_cross_mrca) = \
      prop.get_set(prop.declare_property("cross_mrca"))
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
        set_cross_mrca(x, m)
        return m
      return None
    traverse(AB.A.top)
  AB.get_cross_mrca = get_cross_mrca
  crosses(AB)
  crosses(swap(AB))

def swap(AB):
  BA = AB.swap()
  BA.get_cross_mrca = AB.get_cross_mrca
  BA.A = AB.B
  BA.B = AB.A
  return BA

# -----------------------------------------------------------------------------
# Precompute 'blocks' (tipe sets, implemented in one of various ways).
# A block is represented as either a set or TOP_BLOCK

def compute_blocks(AB):
  def traverse(x, in_left):
    # x is in A or B
    e = BOTTOM_BLOCK
    for c in get_inferiors(x):  # inferiors in A/B
      e = combine_blocks(e, traverse(c, in_left))
    e2 = possible_mtrm_block(AB, in_left(x))
    if (not is_empty_block(e) and not is_empty_block(e2)
        and monitor(x)):
      log("# non-mutual TRM as tipe: %s" % blurb(x))
    e = combine_blocks(e2, e)
    if is_top(x):
      e = TOP_BLOCK
    if e != BOTTOM_BLOCK:
      set_block(x, e)
    #log("# block(%s) = %s" % (blurb(x), e))
    return e
  traverse(AB.A.top, AB.in_left)
  traverse(AB.B.top, AB.in_right)

def possible_mtrm_block(AB, z):
  rel = get_mtrm(z)             # a Relative
  return mtrm_block(z, rel.record) if rel else BOTTOM_BLOCK

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
      return (LT, y, None, "in-checklist")
    else:
      assert False  # can't happen
  else:
    return (DISJOINT, y, "in-checklist")

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
# Load/dump a set of provisional matches (could be either record match
# or taxonomic matches... but basically, record matches).  The matches are stored 
# as Relations under the 'match' property of nodes in AB.

def load_matches(row_iterator, AB):
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
    if x:
      set_match(x, relation(ship, y, None, note))
    if x and y: match_count += 1
    else: miss_count += 1

  #log("# %s matches, %s misses" % (match_count, miss_count))

def reverse_note(note):
  if ' ' in note:
    return "↔ " + note            # tbd: deal with 'coambiguous'
  else:
    return note

def record_match(x):
  return get_matched(x)

