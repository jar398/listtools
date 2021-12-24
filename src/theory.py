#!/usr/bin/env python3

import types
import property as prop, checklist, workspace

from util import log
from checklist import *
from workspace import *
from match_records import match_records

# So how does this work.
# We determine an alignment heuristically.
# Concurrently, given a developing alignment, we determine a spanning
# tree.
# The alignment is represented as...
#   - Relateds values of the 'equated' property
#   - Relateds values of the 'superior' property
#   - A set of dropped nodes (conflicting or garbage)
# umm, that's pretty much it??

def analyze(AB):
  m = match_records(checklist_to_rows(AB.A), checklist_to_rows(AB.B))
  load_matches(m, AB)
  find_MTRMs(AB)
  compute_estimates(AB)
  spanning_tree(AB)

# -----------------------------------------------------------------------------
# Spanning tree computation.

def spanning_tree(AB):
  ensure_levels(AB.A)
  ensure_levels(AB.B)
  def traverse(v, in_leftright):
    z = in_leftright(v)
    if get_superior(z, None): return
    for c in get_inferiors(v):
      traverse(c, in_leftright)
    process_lineage(z,
                    get_left_superior(AB, z),
                    get_right_superior(AB, z),
                    AB) # Goes all the way up to the root
  traverse(AB.B.top, AB.in_right)
  traverse(AB.A.top, AB.in_left)

# x in A, y in B

def process_lineage(z, rx, ry, AB):

  state = [z, rx, ry]

  def propose_superior(z, rs, note, rx, ry):
    assert rs
    assert isinstance(rs, Related)
    if rx: assert isinstance(rx, Related)
    if ry: assert isinstance(ry, Related)
    set_superior(z, Related(rs.relation, rs.other, rs.status, note))
    state[0] = rs.other         # new z
    state[1] = rx
    state[2] = ry

  def propose_equation(x, y, note):
    set_equated(x, Related(EQ, y, "equated A", note))
    set_equated(y, Related(EQ, x, "equated B", note))    # or accepted

  while (rx or ry) and get_superior(z, None) == None:

    while rx and not block_lt(AB, z, rx.other):
      log("# %s not< x=%s, bump x to %s" % (blurb(z),
                                            blurb(rx),
                                            blurb(get_left_superior(AB, rx.other))))
      rx = get_left_superior(AB, rx.other)
    while ry and not block_lt(AB, z, ry.other):
      log("# %s not< x=%s, bump x to %s" % (blurb(z),
                                            blurb(ry),
                                            blurb(get_right_superior(AB, ry.other))))
      ry = get_right_superior(AB, ry.other)
    assert rx or ry

    rp = get_left_superior(AB, rx.other) if rx else None
    rq = get_right_superior(AB, ry.other) if ry else None

    if not ry:
      assert rx
      propose_superior(z, rx, "inherited A", rp, None)
    elif not rx:
      propose_superior(z, ry, "inherited B", rq, None)
    else:
      x = rx.other if rx else None
      y = ry.other if ry else None
      relation = block_relation(AB, x, y)
      p = rp.other if rp else None
      q = rq.other if rq else None
      if relation == LT:            # different blocks
        propose_superior(z, rx, "because MTRM sets A", rp, ry) # No choice
      elif relation == GT:          # different blocks
        propose_superior(z, ry, "because MTRM sets B", rx, rq)
      elif relation == CONFLICT:
        # x conflicts with y.  Delete x, take min(p, q) as parent
        # Parent of z is y, not x; skip x and go right to p
        propose_deprecation(x)
        propose_superior(z, ry, "conflict", rp, rq)    # skip bads in x-chain
      elif block_lt(AB, x, p):                             # COMPARABLE
        if block_lt(AB, y, q):
          propose_equation(x, y, "no reason not to")
          propose_superior(z, ry, "follows", rp, rq)
        else:
          propose_superior(z, ry, "triangle heuristic A", rx, rq)
      elif block_lt(AB, y, q):
        propose_superior(z, x, "triangle heuristic B", rp, ry)
      else:
        n = get_match(y)
        if n and n.relation == EQ:
          if n.other == x:
            propose_equation(x, y, "match + similar")
            propose_superior(z, ry, "follows", rp, rq)
          else:
            propose_superior(z, y, "priority B because hoping for match",
                             rx, rq)
        else:
          m = get_match(x)
          if m and m.relation == EQ:
            propose_superior(z, x, "priority A because hoping for match",
                             rp, ry)
          else:
            # Unmatched against unmatched = toss-up
            propose_superior(z, x, "priority A because toss-up",
                             rp, ry)

    (z, rx, ry) = state
    # end while loop

def get_left_persona(AB, z):
  x = AB.case(z,
              lambda x:z,
              lambda y:None)
  if not x:
    ship = get_equated(z, None)
    if ship: x = ship.other
  return x

def get_right_persona(AB, z):
  y = AB.case(z,
              lambda x:None,
              lambda y:z)
  if not y:
    ship = get_equated(z, None)
    if ship: y = ship.other
  return y

def get_left_superior(AB, z):
  if not z: return None
  x = get_left_persona(AB, z)
  if not x: return None
  ship = AB.case(x,
                 lambda x:get_superior(x, None),
                 lambda y:None)
  if not ship: return None
  return Related(ship.relation, AB.in_left(ship.other), ship.status, ship.note)

def get_right_superior(AB, z):
  if not z: return None
  y = get_right_persona(AB, z)
  if not y: return None
  ship = AB.case(z,
                 lambda x:None,
                 lambda y:get_superior(y, None))
  if not ship: return None
  return Related(ship.relation, AB.in_right(ship.other), ship.status, ship.note)

def equated(x, y):
  if x == y: return True
  z = get_equated(x, None)
  return z and z.other == y

def block_lt(AB, x, y):
  return block_relation(AB, x, y) == LT    # or LE for synonyms ...?

# -----------------------------------------------------------------------------
# Find estimates.  Assumes mtrm has run.

def block_relation(AB, x, y):
  assert x
  assert y
  if in_same_tree(AB, x, y):
    return simple_lt(AB.case(x, lambda x:x, lambda y:y),
                     AB.case(x, lambda x:x, lambda y:y))
  else:
    if is_toplike(x):
      if is_toplike(y): return EQ
      else: return GT
    elif is_toplike(y): return LT
    log("# Compare %s to %s, neither is top" % (blurb(x), blurb(y)))
    return estimate_relation(get_estimate(x), get_estimate(y))

(get_estimate, set_estimate) = prop.get_set(prop.get_property("estimate"))

def compute_estimates(AB):
  def traverse(x, in_left):
    e = null_estimate
    for c in get_inferiors(x):
      e = join_estimates(e, traverse(c, in_left))
    if is_empty(e):
      z = in_left(x)
      w = get_equated(z, None)
      if w:
        e = eta(z, w.other)
        set_estimate(x, e)
    else:
      set_estimate(x, e)
    log("# Estimate(%s) = %s" % (blurb(x), e))
    return e
  traverse(AB.B.top, AB.in_right)
  traverse(AB.A.top, AB.in_left)

# Estimates are sets in this instance, and aren't estimates

null_estimate = set()
def join_estimates(e1, e2): return e1 | e2
def eta(x, y): return {min(x.id, y.id)}
def is_empty(e): return len(e) == 0

def estimate_relation(e1, e2):  # ASSUMES OVERLAP
  if e1.issubset(e2):
    if e1 == e2: return COMPARABLE
    else: return LT
  elif e2.issubset(e1): return GT
  # elif e1.isdisjont(e2): return DISJOINT
  else: return CONFLICT

# -----------------------------------------------------------------------------
# Find mutual tipward matches

def find_MTRMs(AB):
  find_TRMs(AB.A, AB.in_left, set_trm)
  counter = [1]
  def set_if_mutual(z, m):
    set_trm(z, m)           # Nonmutual, might be of interest
    z2 = m.other
    m2 = get_trm(z2, None)      # tipward match to/in A
    if m2 and m2.relation == EQ and m2.other == z:
      if monitor(z): log("# MTRM: %s :=: %s" % (blurb(z), blurb(z2)))
      set_equated(z, m)
      set_equated(z2, m2)
  find_TRMs(AB.B, AB.in_right, set_if_mutual)    # of AB.flip()

(get_trm, set_trm) = prop.get_set(prop.get_property("TRM"))

def find_TRMs(A, in_left, set_trm):
  ensure_inferiors_indexed(A)
  def traverse(x):
    seen = False
    for c in get_inferiors(x):
      seen = traverse(c) or seen
    if not seen and not is_top(x):
      z = in_left(x)
      m = get_match(z)
      if m and m.relation == EQ:
        #if monitor(z): log("# TRM: %s = %s" % (blurb(z), blurb(m.other)))
        set_trm(z, m)
        seen = z
    return seen
  traverse(A.top)

def loser():
  if False: yield "lose"

# -----------------------------------------------------------------------------
# Load/dump a set of provisional matches (could be either record match
# or taxonomic matches... but basically, record matches)

# Fields of match records <A (matched), relation, B (taxon), remark>
get_match_relation = prop.getter(prop.get_property("relation"))
get_matched_key = prop.getter(prop.get_property("matched_id"))
get_match_note = prop.getter(prop.get_property("match_note"))

def load_matches(row_iterator, AB):

  header = next(row_iterator)
  plan = prop.make_plan_from_header(header)
  match_count = 0
  miss_count = 0
  for row in row_iterator:
    # row = [matchID, rel, taxonID, remark]
    match = prop.construct(plan, row)
    x = y = None
    xkey = get_matched_key(match, None)
    if xkey:
      x_in_A = look_up_record(AB.A, xkey)
      if x_in_A:
        x = AB.in_left(x_in_A)
    ykey = get_primary_key(match, None)
    if ykey:
      y_in_B = look_up_record(AB.B, ykey)
      if y_in_B:
        y = AB.in_right(y_in_B) 

    # The columns of the csv file
    rel = rcc5_relation(get_match_relation(match))
    note = get_match_note(match, MISSING)

    # x or y might be None with rel=NOINFO ... hope this is OK

    if x:
      set_match(x, Related(rel, y, "record match", note))
    if y:
      set_match(y, Related(reverse_relation(rel), x, "record match",
                           reverse_note(note)))
    if x and y: match_count += 1
    else: miss_count += 1

  log("# %s matches, %s misses" % (match_count, miss_count))

def reverse_note(note):
  return "↔ " + note            # tbd: deal with 'coambiguous'

def record_match(x):
  ship = get_match(x)
  if ship.relation == EQ: return ship.other
  return None

"""
TBD: filter out seniors
  seniors = 0
      # Filter out senior synonyms here
      if is_senior(item, accepted_item):
        seniors += 1
      else:
  if seniors > 0:     # Maybe interesting
    print("-- Suppressed %s senior synonyms" % seniors,
          file=sys.stderr)

"""

# -----------------------------------------------------------------------------
# Same-tree relations

(get_level, set_level) = prop.get_set(prop.get_property("level"))

def simple_le(x, y):
  y1 = x     # proceeds down to y
  # if {level(x) >= level(y1) >= level(y)}, then x <= y1 <= y or disjoint
  stop = get_level(y)
  while get_level(y1) > stop:
    y1 = get_superior(y1)    # y1 > previously, level(y1) < previously
  # level(y) = level(y1) <= level(x)
  return y1 == y    # Not > and not disjoint

def simple_lt(x, y): return simple_le(x, y) and x != y

def in_same_tree(AB, x, y):
  return (AB.case(x, lambda x:1, lambda x:2) ==
          AB.case(y, lambda x:1, lambda x:2))

def ensure_levels(S):
  def cache(x, n):
    set_level(x, n)
    for c in get_inferiors(x):
      cache(c, n+1)
  cache(S.top, 0)

# -----------------------------------------------------------------------------

if __name__ == '__main__':
  def testit(m, n):
    AB = workspace.workspace_from_newicks(m, n)
    analyze(AB)
    workspace.show_workspace(AB)
  testit("a", "a")
  testit("(c,d)a", "(c,e)b")
