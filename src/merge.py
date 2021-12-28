#!/usr/bin/env python3

import types, argparse
import property as prop, checklist, workspace, theory

from util import log
from checklist import *
from workspace import *
from theory import *
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

def merge(a_iter, b_iter, matches=None):
  A = rows_to_checklist(a_iter, {"name": "A"})  # meta
  B = rows_to_checklist(b_iter, {"name": "B"})  # meta
  AB = analyze(A, B, matches)
  return merge_preorder_rows(AB)

def analyze(A, B, m_iter=None):
  AB = make_workspace(A, B, {"name": "AB"})
  if not m_iter:
    m_iter = match_records(checklist_to_rows(AB.A), checklist_to_rows(AB.B))
  load_matches(m_iter, AB)
  AB.get_cross_mrca = mrca_crosser(AB)
  spanning_tree(AB)
  return AB

def merge_preorder_rows(AB, props=None):
  return checklist.preorder_rows(AB, props or usual_merge_props(AB))

def usual_merge_props(AB):

  def recover_left_id(z, default=None):
    x = get_left_persona(AB, z)
    return get_primary_key(get_outject(x)) if x else default

  def recover_right_id(z, default=None):
    y = get_right_persona(AB, z)
    return get_primary_key(get_outject(y)) if y else default

  return usual_workspace_props + \
     (prop.get_property("taxonID_A", getter=recover_left_id),
      prop.get_property("taxonID_B", getter=recover_right_id))

left_id_prop = prop.get_property("taxonID_A")
right_id_prop = prop.get_property("taxonID_B")
(get_left_id, set_left_id) = prop.get_set(left_id_prop)
get_right_id = prop.getter(right_id_prop)

# -----------------------------------------------------------------------------
# Spanning tree computation.

def spanning_tree(AB):
  find_MTRMs(AB)
  compute_blocks(AB)
  ensure_levels(AB.A)           # N.b. NOT levels in AB
  ensure_levels(AB.B)
  # Do this bottom up.  We'll sometimes hit nodes we've already processed.
  def traverse(v, in_leftright):
    for c in get_inferiors(v):
      traverse(c, in_leftright)
    z = in_leftright(v)
    process_lineage(AB, z) # Goes all the way up to the root
  log("#-- start B pass")
  traverse(AB.B.top, AB.in_right)
  log("#-- start A pass")
  traverse(AB.A.top, AB.in_left)
  ensure_inferiors_indexed(AB)

# x in A, y in B

def process_lineage(AB, z):

  rx = get_left_superior(AB, z)
  ry = get_right_superior(AB, z)

  # Consider rx and ry as possible parent relationships for z

  while z != AB.top:

    if get_superior(z, None):
      #clog("# 2nd(?) or so visit to %s" % (blurb(z),))
      break

    if rx and z == rx.record:
      clog("# Lift rx", rx)
      rx = get_left_superior(AB, rx.record)
    elif ry and z == ry.record:
      clog("# Lift ry", ry)
      ry = get_right_superior(AB, ry.record)

    if True:
      clog("# Candidates for superior of %s are %s, %s" %
           (blurb(z), blurb(rx), blurb(ry)))

    assert rx or ry
    if rx and ry:
      (ship, note) = relationship_per_heuristics(AB, rx.record, ry.record)
    elif rx:
      (ship, note) = (LT, "missing ry")
    elif ry:
      (ship, note) = (GT, "missing rx")
    assert isinstance(ship, int), ship

    rp = get_left_superior(AB, rx.record) if rx else None
    rq = get_right_superior(AB, ry.record) if ry else None

    if ship == EQ:
      propose_equation(AB, rx.record, ry.record, note)
      t = (z, ry, note, rp, rq)
    elif ship == LT or ship == SYNONYM:
      t = (z, rx, note, rp, ry)
    elif ship == CONFLICT:
      # x conflicts with y.  Delete x, take min(p, q) as parent
      # Parent of z is y, not x; skip x and go right to p
      note = "conflict over %s" % blurb(z)
      propose_deprecation(AB, rx.record, ry.record, note)
      t = (z, ry, note, rp, ry)    # z's parent could be y, or p (above the conflict)
    elif ship == GT:
      t = (z, ry, note, rx, rq)
    else:
      assert False, rcc5_symbol(ship)

    (z, rsup, note, rx, ry) = t
    sup = rsup.record
    (ship, status) = consider_unaccepted(AB, z)
    propose_superior(AB, z, rsup, ship, status, note)  # status, note
    z = sup

    # end while


def relationship_per_heuristics(AB, x, y):
  assert isinA(AB, x) and isinB(AB, y)

  rp = get_left_superior(AB, x)
  rq = get_right_superior(AB, y)

  ship = relationship_per_blocks(AB, x, y) # compare blocks
  assert ship != DISJOINT
  if ship == LT:            # different MTRM sets
    ans = (ship, "MTRMs(x) ⊂ MTRMs(y)")
  elif ship == GT:          # different MTRM sets
    ans = (GT, "MTRMs(y) ⊂ MTRMs(x)")
  elif ship == CONFLICT:
    ans = (CONFLICT, "MTRMs(x) >< MTRMs(y)")
  # *** ship is COMPARABLE at this point ***
  elif (rp and rq and
        get_block(x) != get_block(rp.record) and   # in different blocks?
        get_block(y) != get_block(rq.record)):
    ans = (EQ, "squeeze")    # z parent is x≡y
  else:
    n = get_match(y, None)
    if n and n.relationship == EQ:
      if n.record == x:
        ans = (EQ, "name match (%s) and compatible" % (n.note or 'match'))
      else:
        ans = (ship, "priority y because hoping for name match")
    else:
      ans = (GT, "no match, x priority arbitrary")

  return ans

# Tweak relationship to preserve accepted/unaccepted distinction.
# use this when set the superior of z.  not relevant for the x/y
# parent choice.  This is a policy decision.

def consider_unaccepted(AB, z):
  assert isinstance(z, prop.Record)
  if z == AB.top:
    # everything accepted?  or everything a synonym?
    log("# what is %s's relationship to top?" % blurb(z))
    # but there's nothing to relate to.
    return (LT, "at top")
  # If it was a synonym in B, make it one in AB
  rq = get_right_superior(AB, z)
  if rq:
    assert rq.relationship != None
    if rq.relationship != ACCEPTED:
      log("# making %s a synonym in A+B because it's one in B" % blurb(z))
    return (rq.relationship, rq.status)
  rp = get_left_superior(AB, z)
  if rp:
    assert rp.relationship != None
    if rp.relationship != ACCEPTED:
      log("# making %s a synonym in A+B because it's one in A" % blurb(z))
    return (rp.relationship, rp.status)
  else:
    assert False


# -----------------------------------------------------------------------------
# Find mutual tipward matches

def find_MTRMs(AB):
  find_TRMs(AB.A, AB.in_left, set_trm)
  counter = [1]
  def set_if_mutual(z, m):      # z in B
    set_trm(z, m)           # Nonmutual, might be of interest
    z2 = m.record           # z2 in A
    m2 = get_trm(z2, None)      # tipward match to/in A
    if m2 and m2.relationship == EQ and m2.record == z:
      #if monitor(z): log("# MTRM: %s :=: %s" % (blurb(z), blurb(z2)))
      propose_equation(AB, z2, z, "MTRM")
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
      if m and m.relationship == EQ:
        #if monitor(z): log("# TRM: %s = %s" % (blurb(z), blurb(m.record)))
        set_trm(z, m)
        seen = z
    return seen
  traverse(A.top)

# -----------------------------------------------------------------------------
# Load/dump a set of provisional matches (could be either record match
# or taxonomic matches... but basically, record matches)

def load_matches(row_iterator, AB):
  log("# Loading matches")

  header = next(row_iterator)
  plan = prop.make_plan_from_header(header)
  match_count = 0
  miss_count = 0
  for row in row_iterator:
    # row = [matchID, rel, taxonID, remark]
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

    # The columns of the csv file

    rel = rcc5_relationship(get_match_relationship(match)) # EQ, NOINFO
    note = get_match_note(match, MISSING)

    # x or y might be None with rel=NOINFO ... hope this is OK

    if y:
      set_match(y, relation(reverse_relationship(rel), x, "record match",
                           reverse_note(note)))
    if x:
      set_match(x, relation(rel, y, "record match", note))
    if x and y: match_count += 1
    else: miss_count += 1

  log("# %s matches, %s misses" % (match_count, miss_count))

get_match_relationship = prop.getter(prop.get_property("relation", inherit=False))

def reverse_note(note):
  if ' ' in note:
    return "↔ " + note            # tbd: deal with 'coambiguous'
  else:
    return note

def record_match(x):
  rel = get_match(x)
  if rel.relationship == EQ: return rel.record
  return None

# -----------------------------------------------------------------------------

# Load a merged checklist, recovering equated, match, and
# maybe... left/right id links? ... ...
# N.b. this will be a mere checklist, not a workspace.

def rows_to_merged(rows, meta):
  M = checklist.rows_to_checklist(rows, meta)  # sets superiors
  resolve_merge_links(M)
  return M

# This is used for reporting (not for merging).
# When we get here, superiors have already been set, so the ghosts
# will be bypassed in enumerations (preorder and all_records).

def resolve_merge_links(M):
  conjured = {}
  def register(rec):
    conjured[get_primary_key(rec)] = rec
  for record in all_records(M):

    note = get_match_note(record, None)
    key = get_match_key(record, None)
    if key != None:
      z = look_up_record(M, key, record) or conjured.get(key, None)
      if not z:
        z = conjure_record(M, key, record, register)
      if not get_match(z, None):
        set_match(z, relation(EQ, record, "record match",
                              reverse_note(note)))
      # z is record's record match regardless
      set_match(record, relation(EQ, z, "record match", note))
    elif note:
      set_match(record, relation(NOINFO, None, None, note))

    note = get_equated_note(record, None)
    key = get_equated_key(record, None)
    if key != None:
      # Only B records have equated's:
      #  have equation record EQ z, want synonym z EQ record
      z = look_up_record(M, key, record) or conjured.get(key, None)
      if not z:
        z = conjure_record(M, key, record, register)   # in A
      if not get_superior(z, None):
        set_superior(z, relation(EQ, record, "equivalent", note))
      set_equated(record, relation(EQ, z, "equivalent", note))  # B->A
      left_id = get_left_id(record, None)
      if left_id:
        z_left = get_left_id(z, None)
        if z_left: assert z_left == left_id
        else: set_left_id(z, left_id)
  return conjured

# Create new dummy records if needed so that we can link to them

def conjure_record(C, key, record, register):
  log("# Creating record %s on demand" % key)
  rec = checklist.make_record(key, C)
  set_left_id(rec, get_left_id(record))    # ???
  register(rec)
  return rec

#      x0 ---
#             \   match
#              \
#       x  <->  y    SYNONYM / EQ
#        \
#  match  \
#           --- y0

# -----------------------------------------------------------------------------

def test():
  def testit(m, n):
    log("\n-- Test: %s + %s --" % (m, n))

    B = rows_to_checklist(newick.parse_newick(n),
                          {"name": "B"})  # meta
    A = rows_to_checklist(newick.parse_newick(m),
                          {"name": "A"})  # meta

    AB = analyze(A, B)
    workspace.show_workspace(AB, props=usual_merge_props(AB))
  testit("a", "a")              # A + B
  testit("(c,d)a", "(c,e)b")
  testit("((a,b)e,c)d", "(a,(b,c)f)D")

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    TBD.  stdin = the A checklist
    """)
  parser.add_argument('--A', help="the A checklist (defaults to stdin)")
  parser.add_argument('--B', help="the B checklist")
  parser.add_argument('--matches', help="record matches")
  parser.add_argument('--test', action='store_true', help="run tests")
  args=parser.parse_args()

  if args.test:
    test()
  else:
    a_path = args.A
    b_path = args.B
    matches_path = args.matches

    with open(b_path) as b_file:
      def x(a_file):
        # 3.
        def y(matches_iter):
          rows = merge(csv.reader(a_file),
                       csv.reader(b_file),
                       matches=matches_iter)
          writer = csv.writer(sys.stdout)
          for row in rows:
            writer.writerow(row)

        # 2.
        if matches_path:
          with open(matches_path) as matches_file:
            y(csv.reader(matches_file))
        else:
          y(None)

      # 1.
      if a_path:
        with open(a_path) as a_file:
          x(a_file)
      else:
        x(sys.stdin)

