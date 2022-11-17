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

def separated(x, y):
  AB = get_source(x)
  return isinA(AB, x) != isinA(AB, y)

def theorize(AB):
  AB.specimen_records = {}
  analyze_blocks(AB)                  # sets of 'tipes'
  compute_reflections(AB)
  find_equivalents(AB)

# -----------------------------------------------------------------------------
# Given a record v in the A checklist, return the least node in B
# that's definitely greater than v.

# AB.B could be priority, or not

def cross_superior(AB, v):
  # v = in_A(x)
  if monitor(v): log("# cross_superior %s:" % blurb(v))
  v_up = increase_until_overlap(AB, v)
  if not v_up:
    log("no overlap of %s with B" % blurb(v))
    return None

  if monitor(v): log("# xsup loop v = %s <= %s = v_up" % (blurb(v), blurb(v_up)))
  # increase w until v_up is less (< or <=) than it
  w = get_reflection(v_up, BOTTOM)     # candidate in AB
  assert w
  while True:
    if monitor(v): log("#  xsup iter w = %s" % (blurb(w), ))
    ship = cross_relation(AB, v_up, w).relationship
    if ship == LT or ship == LE or ship == PERI:
      break
    if ship == EQ and v_up != v:
      break
    (q, wsyn) = local_sup(AB, w)
    if not q:
      log("#  no node in B is bigger than %s" % blurb(v))
      return None
    w = q

  rel = cross_relation(AB, v, w)
  ship = rel.relationship
  assert ship == LT or ship == LE or ship == PERI, \
    (blurb(v), rcc5_symbol(ship), blurb(w))

  return relation(ship, w, rel.note)

#-----------------------------------------------------------------------------
# cross_relation: The implementation of the RCC-5 theory of AB (A+B).

# v and w are Records in AB, not in A or B
# Returns a Relative to w

# TBD: set status to "accepted" or "synonym" as appropriate

def cross_relation(AB, v, w):
  answer = None

  if isinA(AB, v) == isinA(AB, w):
    rel = simple_relationship(get_outject(v), get_outject(w))
    answer = (rel.relationship, rel.note)
  elif equivalent(v, w):
    answer = (EQ, "equivalent")
  else:
    (p,  vsyn) = local_sup(AB, v) # vsyn iff v is a synonym of p
    (q,  wsyn) = local_sup(AB, w)
    if vsyn and wsyn and equivalent(p, q):
      answer = (OVERLAP, "co-synonyms")
    elif vsyn and equivalent(p, w):
      if monitor(v) or monitor(w):
        log("# %s %s synonym of %s %s" % (blurb(v), blurb(p), blurb(w), blurb(ysup)))
      answer = (SYNONYM, "is a synonym of")
    elif wsyn and equivalent(v, q):
      answer = (SYNONYM, "has synonym")

    b1 = get_block(v, BOTTOM_BLOCK)
    b2 = get_block(w, BOTTOM_BLOCK)

    # Deal with peripheral nodes (empty block)
    v_up = increase_until_overlap(AB, v)
    w_up = increase_until_overlap(AB, w)
    if v_up != v and w_up != w:
      answer = (DISJOINT, "peripherals")
    elif v_up != v:
      answer = (PERI, "left peripheral")
    elif w_up != w:
      answer = (IREP, "right peripheral")

    else:
      # v and w have nonempty blocks.
      ship = block_relationship(b1, b2)
      if ship == EQ:

        # Look for a record match for v or w, and punt to simple case
        rel1 = rel2 = None
        v_eq = get_equivalent(v, None)  # v, v_eq, v_eq.record in AB (from B)
        if v_eq:
          # v = m ? w
          rel1 = simple_relationship(get_outject(v_eq.record), get_outject(w))
        w_eq = get_equivalent(w, None)   # w, w_eq, w_eq.record in AB
        if w_eq:
          # v ? n = w
          rel2 = simple_relationship(get_outject(v), get_outject(w_eq.record))
        # See whether the v/v_eq/w_eq/w diagram commutes
        if w_eq and v_eq:
          assert rel1.relationship == rel2.relationship
          ship = rel1.relationship & rel2.relationship
          if ship != 0:
            answer = (ship, rel1.note)
          else:
            log("shouldn't happen 1 %s %s %s / %s %s %s" %
                (blurb(v),
                 rcc5_symbol(rel1.relationship), 
                 blurb(w_eq.record),
                 blurb(v_eq.record),
                 rcc5_symbol(rel2.relationship), 
                 blurb(w)))
            answer = (OVERLAP, "shouldn't happen 1")
        elif v_eq:
          answer = (rel1.relationship, rel1.note)
          # log("# v %s %s %s" % (blurb(v), blurb(rel1), blurb(w)))
        elif w_eq:
          answer = (rel2.relationship, rel2.note)
          # log("# w %s %s %s" % (blurb(v), blurb(rel2), blurb(w)))

        else:
          assert b1 != BOTTOM_BLOCK
          # Neither v nor w has an equivalent, but they're in same block,
          # so OVERLAP.   (pathological monotype chain.)
          answer = (OVERLAP, "in parallel monotypic chains")

        assert answer[0] != EQ

      else:
        # ship is not EQ (i.e. specimen sets are different)
        answer = (ship, "via specimen set comparison")

  (ship, note) = answer
  if monitor(v) or monitor(w):
    log("# %s %s %s / %s" % (blurb(v), rcc5_symbol(ship), blurb(w), note))

  return relation(ship, w, "articulation", note)

# Returns <p, syn> where p in AB is superior in A (or B) of v,
#  and syn is true iff p is a synonym

def local_sup(AB, v):
  loc = get_superior(get_outject(v), None)
  if not loc: return (None, False)
  if isinA(AB, v):
    return (AB.in_left(loc.record), loc.relationship == SYNONYM)
  elif isinB(AB, v):
    return (AB.in_right(loc.record), loc.relationship == SYNONYM)
  else:
    assert False

# RCC-5 relationship across the two checklists
# x and y are in AB
# Could be made more efficient by skipping unused calculations

def cross_lt(AB, v, w):
  return cross_relation(AB, v, w).relationship == LT

# Find an ancestor of v (in A tree) that overlaps the B tree, i.e.
# has a nonempty block

def increase_until_overlap(AB, v):
  v0 = v
  if is_empty_block(get_block(v, BOTTOM_BLOCK)):
    x = get_outject(v)
    while True:
      if is_top(x): return None
      rel = get_superior(x, None)
      if not rel: return None   # Trees have no overlap !
      x = rel.record
      if isinA(AB, v):
        v = AB.in_left(x)
      else:
        v = AB.in_right(x)
      if not is_empty_block(get_block(v, BOTTOM_BLOCK)):
        # v0 < v = overlaps with B 
        break
  return v

#-----------------------------------------------------------------------------
# Store equivalents on nodes of AB...

def find_equivalents(AB):
  v = AB.in_left(AB.A.top)
  w = AB.in_right(AB.B.top)
  set_equivalent(v, relation(EQ, w, "top"))
  set_equivalent(w, relation(EQ, v, "top"))
  count = 1
  for x in checklist.postorder_records(AB.A):
    v = AB.in_left(x)
    rel = find_equivalent(v)    # -> a Relation
    if rel:
      set_equivalent(v, rel)
      set_equivalent(rel.record, relation(EQ, v, rel.status))
      count += 1
  log("# found %s B nodes with A equivalents" % count)

(get_equivalent, set_equivalent) = prop.get_set(prop.declare_property("equivalent"))

# Given x, returns y such that x = y.  x and y are in AB, not A or B.
# TBD: should cache this.

def find_equivalent(v):
  w = get_reflection(v, None)
  if w and equivalent(v, w):
    return relation(EQ, w, "equivalent")
  else:
    return None

# Are v and w equivalent?

def equivalent(v, w):
  e = really_equivalent(v, w)
  assert e == really_equivalent(w, v), \
    (blurb(v), blurb(w), blurb(get_matched(v)), blurb(get_matched(w)))
  return e

# This gets cached

def really_equivalent(v, w):
  # Same blocks, and
  # {same name or both at tops of monotype chains}
  b1 = get_block(v, BOTTOM_BLOCK)
  b2 = get_block(w, BOTTOM_BLOCK)
  # Is this right???
  if not same_block(b1, b2):
    # Happens all the time
    # log("# not same block %s %s" % (blurb(v), blurb(w)))
    return False
  rel1 = get_matched(v)
  if rel1:
    return rel1.record == w
  rel2 = get_matched(w)
  if rel2:
    return False

  if not is_empty_block(b1):
    # Neither is matched.
    # Propose equality if parents' blocks are both bigger.
    p_rel = get_superior(v, None)
    q_rel = get_superior(w, None)
    if not p_rel and not q_rel:
      log("# matched at top: %s = %s" % (blurb(v), blurb(w)))
      return True
    if not p_rel or not q_rel:
      return False

    if (get_block(p_rel.record, BOTTOM_BLOCK) != b1 and
        get_block(q_rel.record, BOTTOM_BLOCK) != b2):
      log("# matched: %s = %s" % (blurb(v), blurb(w)))
      return True

  # Happens all the time
  log("# not matched: %s != %s, %s" % (blurb(v), blurb(w), blurb(rel1.record)))
  return False


# plantain, pretzels

# -----------------------------------------------------------------------------
# ONLY USED WITHIN THIS FILE

# reflection of v (in AB) = w (in AB) satisfies: 
#   1. if y = reflection(x) then x ~<= y  ?
#   2. if furthermore x' = reflection(y), then
#      (a) if x' <= x then x ≅ y, 
#      (b) if x' > x then x > y.  [why?]
#
# i.e.   x ≲ y ≲ x' ≤ x   →   x ≅ y   (same 'block')
# i.e.   x ≲ y ≲ x' > x   →   x < y
#
# Cached in AB nodes under the 'reflection' property.
# Needed for equivalent and cosuperior calculations

def compute_reflections(AB):
  def do_reflections(AB):
    def traverse(x):            # arg in A, result in B
      v = AB.in_left(x)         # in AB
      probe = get_tipward(v, None)       # in AB
      if probe:
        m = get_outject(probe.record)  # in B
      else:
        m = BOTTOM                # identity for mrca
      for c in get_inferiors(x):
        q = traverse(c)       # in B
        if q:
          m = mrca(q, m)      # in B
      if m == BOTTOM:
        w = BOTTOM
      else:
        w = AB.in_right(m)
        # Elevate m to a name match, if there is one
        # Happens all the time - do not bother logging
        #  log("# promoting %s -> %s (for %s)" %
        #      (blurb(w), blurb(z), blurb(v)))
        nm = get_matched(x)     # returns relation in AB
        if nm:
          z = nm.record         # in AB
          n = get_outject(z)
          if (get_block(w, BOTTOM_BLOCK) == get_block(z, BOTTOM_BLOCK) and
              simple_lt(m, n)):
            w = z               # Override!
            m = n               # in B
      set_reflection(v, w)
      return m                  # m == get_outject(w)
    traverse(AB.A.top)
  do_reflections(AB)
  do_reflections(swap(AB))

def swap(AB):
  BA = AB.swap()
  BA.A = AB.B
  BA.B = AB.A
  BA.specimen_records = AB.specimen_records
  return BA

(get_reflection, set_reflection) = \
  prop.get_set(prop.declare_property("reflection"))

# -----------------------------------------------------------------------------
# Find tipward record matches (TRMs)

(get_tipward, set_tipward) = prop.get_set(prop.declare_property("tipward"))

def analyze_tipwards(AB):
  find_tipwards(AB.A, AB.in_left)
  find_tipwards(AB.B, AB.in_right)

# Postorder iteration over one of the two summands.
# Sets the 'tipward' property is we want to use the node to represent its
# quasi-type specimen.

def find_tipwards(A, in_left):
  def traverse(x):
    seen = seen2 = False
    for c in get_children(x, ()):
      seen = traverse(c) or seen
    for c in get_synonyms(x, ()):
      seen2 = traverse(c) or seen2
    if not seen and not is_top(x):
      # not seen means that some ancestor is tipward
      z = in_left(x)            # z is tipward...
      rel = get_matched(z)  # rel is a Relative...
      if monitor(x): log("tipwards %s %s" % (blurb(x), blurb(rel)))
      if rel:
        assert separated(z, rel.record)
        # z is tipward and has a match, not necessarily tipward ...
        set_tipward(z, rel)
        seen = x
    return seen or seen2
  traverse(A.top)

# -----------------------------------------------------------------------------
# Precompute 'blocks' (tipe sets, implemented in one of various ways).
# A block is represented as either a set or TOP_BLOCK.
# Blocks are stored on nodes in AB.

def analyze_blocks(AB):
  analyze_tipwards(AB)
  def traverse(x, in_left):
    # x is in A or B
    if monitor(x): log("computing block for %s" % (blurb(x),))
    # initial e = specimens from descendants
    e = BOTTOM_BLOCK
    for c in get_inferiors(x):  # inferiors in A/B
      e = combine_blocks(e, traverse(c, in_left))
      if monitor(c): log("got subblock %s -> %s" % (blurb(c), len(e),))
    v = in_left(x)
    if is_top(x):
      e = TOP_BLOCK
    else:
      specimen_id = get_specimen_id(AB, v)
      if specimen_id:
        e = combine_blocks(e, {specimen_id})
    if e != BOTTOM_BLOCK:
      set_block(v, e)
    #log("# block(%s) = %s" % (blurb(x), e))
    return e
  traverse(AB.A.top, AB.in_left)
  traverse(AB.B.top, AB.in_right)

def get_specimen_id(AB, z):
  rel = get_tipward(z, None)
  if rel:
    if isinB(AB, z):
      x = rel.record
      y = z
    else:
      x = z
      y = rel.record
    assert separated(x, y)
    id = get_primary_key(x)
    AB.specimen_records[id] = (x, y)
    return id
  else: return None

# For debugging

def show_specimens(z, tag, AB):
  def foo(vw):
    (v, w) = vw
    return blurb(w)
  log("# %s: {%s}" % (tag, ", ".join(map(foo, specimen_records(AB, z)))))

# The records in the B checklist corresponding to the type specimens
# in the block for z.

def specimen_records(AB, z):
  return map(lambda id: AB.specimen_records[id],
             get_block(z, BOTTOM_BLOCK))

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

def block_ge(e1, e2):
  return e1 >= e2

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

(get_block, set_block) = prop.get_set(prop.declare_property("block"))

# -----------------------------------------------------------------------------
# Same-tree relationships

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
