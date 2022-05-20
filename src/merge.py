#!/usr/bin/env python3

import types, argparse, os
import property as prop, checklist, workspace, theory

from util import log, stdopen
from checklist import *
from workspace import *
from theory import *
from match_records import match_records
from rcc5 import rcc5_symbol

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
  theory.load_matches(m_iter, AB)
  AB.get_cross_mrca = theory.mrca_crosser(AB)
  spanning_tree(AB)
  return AB

def merge_preorder_rows(AB, props=None):
  def filtered_records():
    for z in checklist.preorder_records(AB):
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
  ensure_inferiors_indexed(AB)  # children and synonyms properties

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
      propose_deprecation(AB, z, rx.record, ry.record)
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

# relationship_per_blocks(AB, get_superior(x).record, get_cosuperior(x)) ...

def relationship_per_blocks(AB, x, y):
  ship = block_relationship(get_block(x, BOTTOM_BLOCK),
                            get_block(y, BOTTOM_BLOCK))
  if ship == EQ and get_matched(x) == y:
    return ship
  elif ship == EQ:
    return COMPARABLE
  else:
    return ship

"""
# What on earth is going on here !!??

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
"""

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
    rx = get_equated(z, None)
    if rx:                      # In A too?
      diffs.append(NO_CHANGE)
      rq = get_superior(z, None)
      if rq:               # Not top?
        x = rx.record
        def field_differs(prope):
          get = prop.getter(prope)
          # Field value match?
          name1 = get(y, None)         # in B
          name2 = get(x, None) # in A
          if name1 and name2 and name1 != name2:
            # Similarly, keep if canonical differs.
            # Could do the same for scientificName and/or other fields.
            diffs.append(prope.label)
            return True
          return False
        field_differs(canonical_prop) or \
          field_differs(tipe_prop)
        field_differs(rank_prop)

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
      diffs.append('not in A')   # String matters, see report.py

  else:                   # Not in B
    x = z
    y = None

    a = get_superior(z, None)
    if a and a.relationship == EQ:
      diffs.append(REDUNDANT)  # cf. keep_record
    else:
      diffs.append('not in B')    # filter out ??

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

NO_CHANGE = 'in both'
REDUNDANT = 'merged'

def keep_record(AB, x):
  sup = get_superior(x, None)
  if sup and sup.relationship == EQ:
    y = sup.record
    assert isinB(AB, y)
    # This is an A record that's = to some B record.
    # Get rid of it unless there is some important difference.
    # diffs = find_difference(AB, y, None) ...

    # Only keep it if we get a new canonical name out of it
    can1 = get_canonical(x, None)
    can2 = get_canonical(y, None)

    return can1 and can2 and can1 != can2
  return True

# Returns the RCC5 relationship between x and y.  This is to be called
# only after the trees have been merged, i.e. after all the
# within-block relationships have been established by process_lineage.

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

def equated(x, y):              # Are x and y equated?
  if x == y: return True
  z = get_equated(y, None)
  return z and z.record == x

# -----------------------------------------------------------------------------
# Match building...

def get_accepted(z):  # in priority tree, if possible
  rel = get_superior(z, None)
  if rel and rel.relationship != ACCEPTED:
    return get_accepted(rel.record)
  return z

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
