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

def merge(a_iter, b_iter, m=None):
  A = rows_to_checklist(a_iter, {"name": "A"})  # meta
  B = rows_to_checklist(b_iter, {"name": "B"})  # meta
  AB = analyze(A, B, m)
  log("here")
  return merge_preorder_rows(AB)

def analyze(A, B, m=None):
  AB = make_workspace(A, B, {"name": "AB"})
  if not m:
    m = match_records(checklist_to_rows(AB.A), checklist_to_rows(AB.B))
  load_matches(m, AB)
  spanning_tree(AB)
  return AB

def merge_preorder_rows(AB, props=None):
  return workspace.preorder_rows(AB, props or usual_merge_props(AB))

def usual_merge_props(AB):

  def recover_id_a(z, default=None):
    x = get_left_persona(AB, z)
    return get_primary_key(get_outject(x)) if x else default

  def recover_id_b(z, default=None):
    y = get_right_persona(AB, z)
    return get_primary_key(get_outject(y)) if y else default

  return usual_workspace_props + \
     (prop.get_property("taxonID_A", getter=recover_id_a),
      prop.get_property("taxonID_B", getter=recover_id_b))

# -----------------------------------------------------------------------------
# Spanning tree computation.

def spanning_tree(AB):
  find_MTRMs(AB)
  compute_mtrmsets(AB)
  ensure_levels(AB.A)
  ensure_levels(AB.B)
  def traverse(v, in_leftright):
    z = in_leftright(v)
    if get_superior(z, None): return
    for c in get_inferiors(v):
      traverse(c, in_leftright)
    assert isinstance(z, prop.Record)
    process_lineage(AB, z) # Goes all the way up to the root
  traverse(AB.B.top, AB.in_right)
  traverse(AB.A.top, AB.in_left)

# x in A, y in B

def process_lineage(AB, z):
  rx = get_left_superior(AB, z)
  ry = get_right_superior(AB, z)
  if not (rx or ry): return
  state = [z, rx, ry]

  # Links should only go from most dominant to most dominant
  # b->b, a->b (a not = anything), b->a (similarly), a->a (similarly)

  def propose_sup(z, rs, note, rx, ry):
    propose_superior(AB, z, rs, None, note)
    # TBD: if there's an equation, and z is in A, make z a synonym
    if rx: assert isinstance(rx, Relative)
    if ry: assert isinstance(ry, Relative)
    state[0:3] = (rs.record, rx, ry)

  def propose_deprecation(x, z, rp, ry):
    clog("# Need to deprecate", x)
    state[0:3] = (z, rp, ry)

  def propose_eq(rx, ry, note):
    propose_equation(AB, rx.record, ry.record, note)
    return ry

  def consider(z, rx, ry):
    clog("# Candidates for superior %s are %s, %s" %
         (blurb(z), blurb(rx), blurb(ry)))

    if False:
      while rx and not mtrmset_lt(AB, z, rx.record):
        log("# %s not< x=%s, bump x to %s" %
            (blurb(z), blurb(rx), blurb(get_left_superior(AB, rx.record))))
        rx = get_left_superior(AB, rx.record)
      while ry and not mtrmset_lt(AB, z, ry.record):
        log("# %s not< x=%s, bump x to %s" %
            (blurb(z), blurb(ry), blurb(get_right_superior(AB, ry.record))))
        ry = get_right_superior(AB, ry.record)
    elif False:
      assert (not rx) or mtrmset_lt(AB, z, rx.record)
      assert (not ry) or mtrmset_lt(AB, z, ry.record)

    assert rx or ry

    rp = get_left_superior(AB, rx.record) if rx else None
    rq = get_right_superior(AB, ry.record) if ry else None

    if not ry:
      propose_sup(z, rx, None, rp, None)
    elif not rx:
      propose_sup(z, ry, None, rq, None)
    else:
      x = rx.record if rx else None
      y = ry.record if ry else None
      ship = mtrmset_relationship(AB, x, y)
      assert ship != DISJOINT
      if ship == LT:            # different mtrmsets
        propose_sup(z, rx, "MTRMs(x) ⊂ MTRMs(y)", rp, ry) # No choice
      elif ship == GT:          # different mtrmsets
        propose_sup(z, ry, "MTRMs(y) ⊂ MTRMs(x)", rx, rq)
      elif ship == CONFLICT:
        # x conflicts with y.  Delete x, take min(p, q) as parent
        # Parent of z is y, not x; skip x and go right to p
        propose_deprecation(x, z, rp, ry)
      elif not rp:
        if rq:
          propose_sup(z, ry, None, rx, rq)
        else:
          rw = propose_eq(rx, ry, "top")
          propose_sup(z, rw, None, rp, rq)
      elif not rq:
        propose_sup(z, x, None, rp, ry)
      else:
        assert ship == COMPARABLE
        n = get_match(y)
        if n and n.relationship == EQ:
          if n.record == x:
            rw = propose_eq(rx, ry, "similar|%s" % (n.note or 'match'))
            propose_sup(z, rw, "parents match", rp, rq)
          else:
            propose_sup(z, y, "priority B because hoping for match",
                        rx, rq)
        else:
          m = get_match(x)
          if m and m.relationship == EQ:
            propose_sup(z, rx, "priority A because hoping for match", 
                        rp, ry)
          elif True:
            # Unmatched against unmatched = toss-up
            propose_sup(z, rx, "tossup A", rp, ry)
          else:
            # Doesn't seem to work so well
            rw = propose_eq(rx, ry, "shot in the dark")
            propose_sup(z, rw, "shot in the dark", rp, rq)
    return state

  while (rx or ry) and get_superior(z, None) == None:
    (z, rx, ry) = consider(z, rx, ry)

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

def test():
  def testit(m, n):
    log("\n-- Test: %s + %s --" % (m, n))

    B = rows_to_checklist(newick.parse_newick(n),
                          {"name": "B"})  # meta
    A = rows_to_checklist(newick.parse_newick(m),
                          {"name": "A"})  # meta

    AB = analyze(A, B)
    workspace.show_workspace(AB, props=usual_merge_props(AB))
  testit("a", "a")
  testit("(c,d)a", "(c,e)b")
  testit("((a,b)e,c)d1", "(a,(b,c)f)d2")

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    TBD.  stdin = the A checklist
    """)
  parser.add_argument('--target', help="the B checklist")
  parser.add_argument('--matches', help="record matches")
  parser.add_argument('--test', action='store_true', help="run tests")
  args=parser.parse_args()

  if args.test:
    test()
  else:
    a_file = sys.stdin
    b_path = args.target
    rm_sum_path = args.matches

    a_iter = csv.reader(a_file)
    with open(b_path) as b_file:
      b_iter = csv.reader(b_file)
      # TBD: Compute sum if not provided ?
      if rm_sum_path:
        with open(rm_sum_path) as rm_sum_file:
          rm_sum_iter = csv.reader(rm_sum_file)
          rows = merge(a_iter, b_iter, rm_sum_iter)
      else:
        rows = merge(a_iter, b_iter)

    writer = csv.writer(sys.stdout)
    for row in rows:
      writer.writerow(row)
