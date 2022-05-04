#!/usr/bin/env python3

import types, functools
import property as prop, checklist, workspace

from util import log
from checklist import *
from workspace import *
from match_records import match_records
from rcc5 import rcc5_relationship, reverse_relationship

# Satisfies: 
#   1. if y = cross_mrca(x) then x ~<= y
#   2. if furthermore x' = cross_mrca(y), then
#      (a) if x' <= x then x ≅ y, 
#      (b) if x' > x then x > y.  [why?]
#
# i.e.   x ≲ y ≲ x' ≤ x   →   x ≅ y   (same 'block')
# i.e.   x ≲ y ≲ x' > x   →   x < y
#
# Cached under the 'cross_mrca' property.

def mrca_crosser(AB, equatee):

  # Cached
  def compute_cross_mrca(z):
    e = equatee(AB, z)
    if e: return e
    (x, in_left, in_right) = \
      AB.case(z,
              lambda x:(x, AB.in_left, AB.in_right),
              lambda y:(y, AB.in_right, AB.in_left))
    m = BOTTOM                  # identity for mrca
    for c in get_inferiors(x):
      q = get_cross_mrca(in_left(c))
      if q:
        r = AB.case(q, lambda x:x, lambda y:y)
        m = mrca(r, m)
    #clog("# cross-mrca of %s is %s" % (blurb(x), blurb(m)))
    return in_right(m) if m != BOTTOM else None
  get_cross_mrca = \
    prop.getter(prop.declare_property("cross_mrca",
                                      filler=compute_cross_mrca))
  return get_cross_mrca

def isinA(AB, z):
  return AB.case(z, lambda x: True, lambda x: False)

def isinB(AB, z):
  return AB.case(z, lambda x: False, lambda x: True)

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
  traverse(AB.B.top, AB.in_right)
  traverse(AB.A.top, AB.in_left)

def possible_mtrm_block(AB, z):
  rel = get_mtrm(z)             # a Relative
  return mtrm_block(z, rel.record) if rel else BOTTOM_BLOCK

(get_tipward, set_tipward) = prop.get_set(prop.declare_property("tipward"))

def get_accepted(z):  # in priority tree, if possible
  rel = get_superior(z, None)
  if rel and rel.relationship != ACCEPTED:
    return get_accepted(rel.record)
  return z

def lt_per_blocks(AB, x, y):
  b1 = get_block(x, BOTTOM_BLOCK)
  b2 = get_block(x, BOTTOM_BLOCK)
  if b1 == TOP_BLOCK: return False    # top not < anything
  elif b2 == TOP_BLOCK: return True   # non-top < top
  else: return b1 < b2

def relationship_per_blocks(AB, x, y):
  ship = block_relationship(get_block(x, BOTTOM_BLOCK),
                            get_block(y, BOTTOM_BLOCK))
  return COMPARABLE if ship == EQ else ship

# Implementation of blocks as Python sets of 'tipes.'

def block_relationship(e1, e2):   # can assume overlap
  if e1 == e2: return EQ          # same blocks
  if e1 == TOP_BLOCK: return GT
  if e2 == TOP_BLOCK: return LT
  elif e1.issubset(e2): return LT
  elif e2.issubset(e1): return GT
  elif e1.isdisjoint(e2): return DISJOINT
  else: return CONFLICT

BOTTOM_BLOCK = set()
TOP_BLOCK = True
def combine_blocks(e1, e2):
  if e1 == BOTTOM_BLOCK: return e2
  if e2 == BOTTOM_BLOCK: return e1
  return e1 | e2
def same_block(e1, e2): return e1 == e2
def is_empty_block(e): return e == BOTTOM_BLOCK
def mtrm_block(x, y):
  v1 = get_tipe(x, None)
  if v1 and v1 == get_tipe(y, None): return {v1}
  v1 = get_scientific(x, None)
  if v1 and v1 == get_scientific(y, None): return {v1}
  v1 = get_canonical(x, None)
  if v1 and v1 == get_canonical(y, None): return {v1}
  v1 = get_stemmed(x, None)     # 'Balaenoptera omurai and' in DH 1.1
  if v1 and v1 == get_stemmed(y, None): return {v1}
  v1 = get_managed_id(x, None)
  if v1 and v1 == get_managed_id(y, None): return {v1}
  # Phyllostoma bilabiatum / Phyllostomus bilabiatum  ... hmph
  # Why do these match? Dermanura anderseni, Artibeus anderseni*  ????
  # log("# Why do these match? %s, %s" % (blurb(x), blurb(y)))
  return {min(x.id, y.id)}

(get_block, set_block) = prop.get_set(prop.declare_property("block"))

# -----------------------------------------------------------------------------
# Alignment building...

# Propose that rs.record should be the parent (superior) of z

def propose_superior(AB, z, rs, ship, status, note):
  assert rs
  assert isinstance(z, prop.Record), blurb(z)
  assert isinstance(rs, Relative)
  assert ship == ACCEPTED or ship == SYNONYM
  # I don't feel good about these
  s = get_accepted(rs.record)
  rel = relation(ship,
                 s,
                 status or rs.status,    # accepted, etc
                 note)
  set_superior_carefully(z, rel)  # explanation of why < or <=

# Propose that x (in A) = y (in B).  x becomes an EQ synonym of y.

def propose_equation(AB, x, y, why_equiv):
  # Polarize
  assert isinA(AB, x) and isinB(AB, y)
  (x, y) = AB.case(x, lambda xx: (x, y), lambda yy: (y, x))

  # Treat the A record like a synonym of the B record
  set_superior_carefully(x, relation(EQ, y, "equivalent", why_equiv))
  # Set uppointer from B to A
  set_equated(y, relation(EQ, x, "equivalent", why_equiv))
  sci = get_scientific(x, None)
  if not get_scientific(y, None) and sci:
    if sci.startswith(get_canonical(y, None)):
      set_scientific(y, sci)
      #log("# transferring '%s' from A to B" % sci)

# x and y are candidate parents for node z, and they conflict with one
# another.  y has priority.

def propose_deprecation(AB, z, x, y):
  assert isinA(AB, x) and isinB(AB, y)
  # set_alt_parent(z, x) ...
  set_conflict(x, y)
  if False:
    log("# Deprecating %s because it conflicts with %s" %
        (blurb(x), blurb(y)))
    log("#  ... %s is now the 'accepted name' of %s" %
        (blurb(yy), blurb(x)))

(get_conflict, set_conflict) = prop.get_set(prop.declare_property("conflict"))

# -----------------------------------------------------------------------------
# Same-tree relations

#    x1 ? y1     'larger'
#   /       \
#  x         y   'smaller'

def simple_relationship(x, y):             # Within a single tree
  (x, synx) = nip_synonym(x)
  (y, syny) = nip_synonym(y)
  if x == y:
    if synx or syny:
      if synx and syny: return NOINFO         # blah
      else: return LE if synx else GE
    else:
      return EQ
  (x1, y1) = find_peers(x, y)    # Decrease levels as needed
  if x1 == y1:     # x <= x1 = y1 >= y
    if x1 == x:     # x = x1 = y1 > y, x > y
      return GT
    elif y1 == y:
      return LT
    else:
      assert False  # can't happen
  else:
    return DISJOINT

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

# This (?) is run pre-merge so should not chase A->B synonym links.
# If the parent is a synonym, that's a problem - shouldn't happen.

def get_left_superior(AB, z):
  x = get_left_persona(AB, z)
  if not x: return None
  rel_in_A = get_superior(x, None)
  # Return relation that's same as rel_in_A except in AB instead of A
  if rel_in_A:
    p = AB.in_left(rel_in_A.record)
    return relation(rel_in_A.relationship, p, rel_in_A.status, rel_in_A.note)
  return None

def get_right_superior(AB, z):
  y = get_right_persona(AB, z)
  if not y: return None
  rq = get_superior(y, None)
  if rq:
    q = AB.in_right(rq.record)
    return relation(rq.relationship, q, rq.status, rq.note)
  return None

# Returns a record in the left checklist, or None

def get_left_persona(AB, z):    # left = A = nonpriority
  assert isinstance(z, prop.Record)
  assert get_outject(z, None)
  return AB.case(z,
                 lambda x:x,       # if z in A, then z
                 lambda y:equal_persona(get_equated(z, None)))

# Returns a record in the right checklist, or None

def get_right_persona(AB, z):   # right = B = priority
  return AB.case(z,
                 lambda x:equal_persona(get_superior(z, None)),
                 lambda y:y)    # if z in B, then z

def equal_persona(rel):
  return (get_outject(rel.record)
          if rel and rel.relationship == EQ
          else None)

def equated(x, y):              # Are x and y equated?
  if x == y: return True
  z = get_equated(y, None)
  return z and z.record == x

# -----------------------------------------------------------------------------
# Find mutual tipward record matches (MTRMs)

def analyze_tipwards(AB):
  find_tipwards(AB.A, AB.in_left)
  find_tipwards(AB.B, AB.in_right)

# Postorder iteration over one of the two summands.

def find_tipwards(A, in_left, finish):
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

def get_mtrm(z):
  rel = get_tipward(z, None)
  if rel:
    q = rel.record
    rel2 = get_tipward(q, None)
    if rel2 and rel2.record == z:
      return q
  return None

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
      set_match(y, relation(reverse_relationship(ship), x, "nominal match",
                            reverse_note(note)))
    if x:
      set_match(x, relation(ship, y, "nominal match", note))
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

