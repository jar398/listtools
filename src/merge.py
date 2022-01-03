#!/usr/bin/env python3

import types, argparse, os
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

def merge(a_iter, b_iter, A_name='A', B_name='B', matches=None):
  A = rows_to_checklist(a_iter, {"name": A_name})  # meta
  B = rows_to_checklist(b_iter, {"name": B_name})  # meta
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
  def filtered_records():
    for z in checklist.preorder_records(AB):  #checklist.all_records(AB):
      if keep_record(AB, z):
        yield z
  return records_to_rows(AB, filtered_records(),
                         props or usual_merge_props(AB))

def usual_merge_props(AB):

  def get_difference(z, default=None):
    return find_difference(AB, z, default)

  return usual_workspace_props + \
     (prop.declare_property("difference", filler=get_difference),)

difference_prop = prop.declare_property("difference")

# -----------------------------------------------------------------------------
# Spanning tree computation.

def spanning_tree(AB):
  analyze_tipwards(AB)                # also find tipes
  compute_blocks(AB)
  ensure_levels(AB.A)           # N.b. NOT levels in AB
  ensure_levels(AB.B)
  # Do this bottom up.  We'll sometimes hit nodes we've already processed.
  def traverse(v, in_leftright):
    u = get_conflict(v, None)
    if u:
      log("# conflict: %s >< %s" % (blurb(u), blurb(v)))
    for c in get_inferiors(v):
      traverse(c, in_leftright)
    z = in_leftright(v)
    process_lineage(AB, z) # Goes all the way up to the root
  log("-- merge: start B pass")
  traverse(AB.B.top, AB.in_right)
  log("-- merge: start A pass")
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
      #clog("# Lift rx", rx)
      rx = get_left_superior(AB, rx.record)
    elif ry and z == ry.record:
      #clog("# Lift ry", ry)
      ry = get_right_superior(AB, ry.record)

    if False:
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
  if ship == LT:            # different tipe sets
    ans = (LT, "types(x) ⊂ types(y)")
  elif ship == GT:          # different tipe sets
    ans = (GT, "types(y) ⊂ types(x)")
  elif ship == CONFLICT:
    ans = (CONFLICT, "types(x) >< types(y)")
  # *** ship is COMPARABLE at this point ***
  elif (rp and rq and
        get_block(x, BOTTOM_BLOCK) != get_block(rp.record, BOTTOM_BLOCK) and   # in different blocks?
        get_block(y, BOTTOM_BLOCK) != get_block(rq.record, BOTTOM_BLOCK)):
    ans = (EQ, "squeeze")    # z parent is x≡y
  else:
    n = get_matched(y)
    if n:
      if n == x:
        ans = (EQ, "name match (%s) and compatible" % (get_basis_of_match(y, None) or '...'))
      else:
        ans = (LT, "priority y because hoping for name match")
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
    if False and rq.relationship != ACCEPTED:
      log("# making %s a synonym in A+B because it's one in B" % blurb(z))
    return (rq.relationship, rq.status)
  rp = get_left_superior(AB, z)
  if rp:
    assert rp.relationship != None
    if False and rp.relationship != ACCEPTED:
      log("# making %s a synonym in A+B because it's one in A" % blurb(z))
    return (rp.relationship, rp.status)
  else:
    assert False


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

# -----------------------------------------------------------------------------
# Load/dump a set of provisional matches (could be either nominal match
# or taxonomic matches... but basically, nominal matches)

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
    note = get_basis_of_match(match, MISSING)
    # x or y might be None with rel=NOINFO ... hope this is OK
    if y:
      set_match(y, relation(reverse_relationship(rel), x, "nominal match",
                            reverse_note(note)))
    if x:
      set_match(x, relation(rel, y, "nominal match", note))
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

# -----------------------------------------------------------------------------

# Returns a string summarizing the main way in which a B node differs from
# its corresponding A node

"""

extension
   not in B x in both x not in A
   synonym x accepted
   parent=y

name
   named(x) rcc5 named(x)


"""


def find_difference(AB, z, default=None):
  diffs = []
  if isinB(AB, z):
    y = z
    x = None  # overridden
    rq = get_superior(z, None)
    if rq:               # Not top?
      rx = get_equated(z, None)
      if rx:                      # In A too?
        x = rx.record
        def check_field(prope):
          get = prop.getter(prope)
          # Field value match?
          name1 = get(y, None)         # in B
          name2 = get(x, None) # in A
          if name1 and name1 != name2:
            # Similarly, keep if canonical differs.
            # Could do the same for scientificName and/or other fields.
            diffs.append(prope.label)
            return True
          return False
        check_field(canonical_prop) or \
          check_field(scientific_prop)
        check_field(rank_prop)

        # What about their parents?
        rp = get_left_superior(AB, z)  # y -> x -> p
        if False:
          p = dequate(rp.record) if rp else None
          q = rq.record if rq else None
          if p != q:
            diffs.append('superior')
        if rp.relationship != rq.relationship:
          if rp.relationship == ACCEPTED:
            diffs.append('synonym/accepted')
          else:
            diffs.append('accepted/synonym')
      else:
        # Present in B but not in A
        diffs.append("not in A")   # String matters, see report.py

  else:                   # Not in B
    x = z
    y = None
    diffs.append("not in B")    # filter out ??

  if x:
    # Say something about mismatched equivalent nodes.
    m = get_matched(x)       # m in B (ry is also in B)
    if m:
      ship = post_hoc_relationship(AB, x, m)
      if ship != EQ and ship != NOINFO:
        diffs.append("use of A name (in A %s in B)" %
                     rcc5_symbol(ship))

  if y:
    m = get_matched(y)       # n in A (rx is also in A)
    if m:
      ship = post_hoc_relationship(AB, m, y)
      if ship != EQ and ship != NOINFO:
        diffs.append("use of B name (in A %s in B)" %
                     rcc5_symbol(ship))
  if len(diffs) == 0:             # In both, no differences
    return None
  else:
    return ';'.join(diffs)

def post_hoc_relationship(AB, x, y):

  (x, synx) = nip_synonym(x)
  (y, syny) = nip_synonym(y)
  if x == y:
    if synx or syny:
      if synx and syny: return NOINFO         # blah
      else: return LE if synx else GE
    else:
      return EQ

  rel = relationship_per_blocks(AB, x, y)
  if rel == COMPARABLE:
    return simple_relationship(x, y) # ordering within block
  else:
    return rel

def keep_record(AB, x):
  sup = get_superior(x, None)
  if sup and sup.relationship == EQ:
    y = sup.record
    assert isinB(AB, y)
    # This is an A record that's = to some B record.
    # Get rid of it unless there is some important difference.
    diffs = find_difference(AB, y, None)
    #if diffs: log("# Keeping A record %s (= %s) because %s" %
    #              (blurb(x), blurb(sup.record), diffs))
    return diffs
  return True


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

class Stdin:
  def __enter__(self): return sys.stdin
  def __exit__(self, exc_type, exc_val, exc_tb): return

def stdopen(x):
  if x == '-':
    return Stdin()
  else:
    return open(x)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Only one of A or B can come from standard input.
    Checklist names if not provided come from removing the extension
    (usually '.csv') from the basename of the path name.
    """)
  parser.add_argument('--A', help="the A checklist, as path name or -",
                      default='-')
  parser.add_argument('--B', help="the B checklist, as path name or -",
                      default='-')
  parser.add_argument('--Aname',
                      help="short durable name of the A checklist")
  parser.add_argument('--Bname',
                      help="short durable name of the B checklist")
  parser.add_argument('--matches', help="file containing match-records output")
  parser.add_argument('--test', action='store_true', help="run smoke tests")
  args=parser.parse_args()

  if args.test:
    test()
  else:
    a_path = args.A
    b_path = args.B
    a_name = args.Aname or os.path.splitext(os.path.basename(a_path))[0]
    b_name = args.Aname or os.path.splitext(os.path.basename(b_path))[0]
    assert a_path != b_path
    matches_path = args.matches

    with stdopen(a_path) as a_file:
      with stdopen(b_path) as b_file:

        def y(matches_iter):
          rows = merge(csv.reader(a_file),
                       csv.reader(b_file),
                       A_name = a_name,
                       B_name = b_name,
                       matches=matches_iter)
          writer = csv.writer(sys.stdout)
          for row in rows:
            writer.writerow(row)

        if matches_path:
          with open(matches_path) as matches_file:
            y(csv.reader(matches_file))
        else:
          y(None)
