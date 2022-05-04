#!/usr/bin/env python3

import types, functools
import property as prop, checklist, workspace

from util import log
from checklist import *
from workspace import *
from match_records import match_records

# ---- tools for working with sum taxonomy A+B

# This is run pre-merge so should not chase A->B synonym links.
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

def isinB(AB, z):
  return AB.case(z, lambda x: False, lambda x: True)

def isinA(AB, z):
  return AB.case(z, lambda x: True, lambda x: False)

# -----------------------------------------------------------------------------
# Decide about relations between two trees

# Returns a pair (overlap, outlier)
#   overlap is in y and x {< = > ><} y (but not !)
#   outlier is in y and x {< >< !} y (but not {> =}) x ⋧ y
#   either can be None
# x is in the A checklist, y is in the B checklist.

def analyze(AB, x, y):
  assert False   # Not in use yet!
  (over, out) = (None, None)
  for d in get_inferiors(y):    # in B
    (over2, out2) = analyze(x, d)
    over = over or over2; out = out or out2
    if over and out: return (over, out)
  # TBD: if over and out are both present, and one of the two is a synonym ...
  if over or out:      # Exclude peripherals
    return (over, out)
  m = AB.record_match(y)
  j = AB.join(m, y)
  return (j, None) if m else (None, j)

def outlier(AB, x, y):
  (_, out) = analyze(x, y)
  return out

# Cache this.
# Satisfies: if y = shadow(x) then x ~<= y, and if x' = shadow(y),
# then if x' <= x then x ~= y, or if x' > x then x > y.
#
# i.e.   x ≲ y ≲ x' ≤ x   →   x ≅ y   (same 'block')
# i.e.   x ≲ y ≲ x' > x   →   x < y

def mrca_crosser(AB):

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


# RCC5 decision procedure, where x and y are in different sources.

def cross_ge(AB, x, y):
  return not outlier(AB, x, y)

def cross_eq(AB, x, y):
  # see if they're peers ??
  return cross_ge(x, y) and cross_ge(y, x)

def gt(AB, x, y):
  return cross_ge(x, y) and not cross_ge(y, x)

# -----------------------------------------------------------------------------
# Precompute 'blocks' (tipe sets, implemented in one of various ways).
# A block is represented as either a set or TOP_BLOCK or BOTTOM_BLOCK

def compute_blocks(AB):
  def traverse(x, in_left):
    # x is in A or B
    e = BOTTOM_BLOCK
    for c in get_inferiors(x):  # inferiors in A/B
      e = join_blocks(e, traverse(c, in_left))
    if LIBERAL_TIPE_SETS or is_empty_block(e):
      e2 = possible_tipe(AB, in_left(x))
      if (not is_empty_block(e) and not is_empty_block(e2)
          and monitor(x)):
        log("# non-mutual TRM as tipe: %s" % blurb(x))
      e = join_blocks(e2, e)
    if is_top(x):
      e = TOP_BLOCK
    if e != BOTTOM_BLOCK:
      set_block(x, e)
    #log("# block(%s) = %s" % (blurb(x), e))
    return e
  traverse(AB.B.top, AB.in_right)
  traverse(AB.A.top, AB.in_left)

LIBERAL_TIPE_SETS = True

def possible_tipe(AB, z):       # tipe as in type specimen/series
  if LIBERAL_TIPE_SETS:
    w = get_cotipe(z, None)
  else:
    w = equatee(AB, z)
  return mtrm_block(z, w) if w else BOTTOM_BLOCK

(get_cotipe, set_cotipe) = prop.get_set(prop.declare_property("cotipe"))

def equated(x, y):              # Are x and y equated?
  if x == y: return True
  z = get_equated(y, None)
  return z and z.record == x

def equatee(AB, z):                 # in other tree
  if isinB(AB, z):
    ship = get_equated(z, None)  # Does z represent an MTRM ?
    return ship.record if ship else None
  else:
    ship = get_superior(z, None)
    if ship and ship.relationship == EQ:
      return ship.record    # z is in A
    return None

def get_accepted(z):  # in priority tree, if possible
  rel = get_superior(z, None)
  if rel and rel.relationship != ACCEPTED:
    return get_accepted(rel.record)
  return z

# Assuming x and y are in different trees, and that we are only
# looking at the 'blocks' (not topology otherwise), what is the
# relationship between x and y?

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
def join_blocks(e1, e2):
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
  # NO GOOD.
  if False:
    yy = AB.get_cross_mrca(x)
    set_superior_carefully(x, relation(SYNONYM, yy, "conflicting", note))
  else:
    # set_alt_parent(z, x) ...
    set_conflict(x, y)
  if False:
    log("# Deprecating %s because it conflicts with %s" %
        (blurb(x), blurb(y)))
    log("#  ... %s is now the 'accepted name' of %s" %
        (blurb(yy), blurb(x)))

(get_conflict, set_conflict) = prop.get_set(prop.declare_property("conflict"))

# -----------------------------------------------------------------------------
# Same-tree relationships (with a single tree)

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
# Find mutual tipward matches

def analyze_tipwards(AB):
  find_cotipes(AB.A, AB.in_left, lambda x,y:None)
  counter = [1]
  def finish(z, m):             # z in B
    # z is tipward.  If m is too, we have a MTRM.
    if get_cotipe(m, None) == z:
      #if monitor(z): log("# MTRM: %s :=: %s" % (blurb(z), blurb(m)))
      propose_equation(AB, m, z, "MTRM")
  find_cotipes(AB.B, AB.in_right, finish)    # of AB.flip()

def find_cotipes(A, in_left, finish):
  ensure_inferiors_indexed(A)
  def traverse(x):
    seen = False
    for c in get_inferiors(x):
      seen = traverse(c) or seen
    if not seen and not is_top(x):
      z = in_left(x)
      m = get_matched(z)
      if m:
        #if monitor(z): log("# cotipe: %s = %s" % (blurb(z), blurb(m)))
        finish(z, m)
        set_cotipe(z, m)
        set_cotipe(m, z)
        seen = z
    return seen
  traverse(A.top)

