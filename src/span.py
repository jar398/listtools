#!/usr/bin/env python3

# Create a spanning tree for A+B based on B priority

# For each equivalent pair, only the B record will have a superior and inferiors.

# Would it be worthwhile to have a way to write a spanning tree as a
# checklist, so that it can be re-read without loss? ...
# Actually, would it be worthwhile to be able to do that with a workspace ...

import argparse, sys
import theory
from checklist import rows_to_checklist, checklist_to_rows, get_superior, \
  set_superior, relation, get_children, set_children, get_synonyms, set_synonyms, \
  get_source, get_outject, blurb
import newick
import checklist, match_records, workspace
from rcc5 import *
from util import log

# Assign a parent to every node of AB except for the top.

# B checklist has priority
# Each node z has two possible parents p and q, one in each checklist.
# (If z is equivalent to z' in the other checklist then by
# construction the possible parents of the two are the same.)
# The 'correct' parent r would be the least one that is greater than both:
# r = p if r < p < q, r = q if r < q < p.
# If p = q then pick the one in the priority checklist.
# If p and q are incomparable (>< or !) then the parent is the one
# from the same checklist (r = p if z in A, r = q if z in B).

def span(AB):

  AB.top = None                 # Updated further down

  def traverse(C):              # C is AB or swap(AB)

    for y in checklist.preorder_records(C.B):    # starts with top

      # y could be A.top or B.top
      w = C.in_right(y)
      assert get_superior(w, None) == None

      # A.top is AB.B.top and has no superior
      if y == AB.B.top:
        assert not AB.top
        AB.top = w
        continue

      sup = None

      # Normalize non-dominant w to dominant equivalent
      if theory.isinA(AB, w):
        v_rel = theory.get_equivalent(w, None)
        # N.b. w and v could be top
        if v_rel:    # in A
          # If w has an equivalent in B, synonymize the two
          sup = relation(SYNONYM, v_rel.record, "equivalent",
                         "treat equivalent as synonym")
      if not sup:
        # p is one possible parent of w
        p_rel = theory.cross_superior(C, w)

        # q is another possible parent of w
        qy_rel = get_superior(y, None)

        # If there are two possibilities, pick the smallest
        if qy_rel and p_rel:
          p = p_rel.record
          q = C.in_right(qy_rel.record)
          pq_rel = theory.cross_relation(C, p, q)
          ship = pq_rel.relationship
          blot = "%s %s %s" % (blurb(p), rcc5_symbol(ship), blurb(q))
          if ship == GT or ship == GE:
            sup = relation(qy_rel.relationship, q,
                           "%s, choosing right" % blot)
          elif ship == LT or ship == LE:
            sup = relation(p_rel.relationship, p,
                           "%s, choosing left" % blot)
          elif ship == EQ:
            if theory.isinB(AB, q):
              sup = relation(qy_rel.relationship, q,
                             "%s, choosing left" % blot)
            else:
              sup = relation(p_rel.relationship, p,
                             "%s, choosing right" % blot)
          else:
            log("incomparable parent candidates: %s" % blot)
            sup = relation(qy_rel.relationship, q,
                           "%s, choosing right" % blot)

        elif qy_rel:
          sup = relation(qy_rel.relationship,
                         C.in_right(qy_rel.record),
                         "root 1")
        elif p_rel:
          sup = relation(p_rel.relationship,
                         p_rel.record,
                         "root 2")
        else:
          log("no parent candidates!? %s" % blurb(w))

        # Normalize sup to dominant equivalent (transfer children and
        # synonyms over to dominant node)
        r = sup.record
        s_rel = theory.get_equivalent(r, None)
        if s_rel and theory.isinB(AB, s_rel.record):
          sup = relation(sup.relationship, s_rel.record, sup.note + " normalized")

      log("sup %s = %s" % (blurb(w), blurb(sup)))

      assert isinstance(sup, checklist.Relative)
      checklist.link_superior(w, sup) # adds w to children/synonyms list

  traverse(AB)
  traverse(theory.swap(AB))
  assert AB.top
  AB.indexed = True
  return AB

# -----------------------------------------------------------------------------

def test():
  def testit(m, n):
    log("\n-- Test: %s + %s --" % (m, n))
    B = rows_to_checklist(newick.parse_newick(n),
                          {"name": "B"})  # meta
    A = rows_to_checklist(newick.parse_newick(m),
                          {"name": "A"})  # meta
    matches_iter = match_records.match_records(checklist_to_rows(A), checklist_to_rows(B))
    AB = workspace.make_workspace(A, B, {"name": "AB"})
    theory.load_matches(matches_iter, AB)
    theory.analyze_tipwards(AB)                # also find tipes
    theory.compute_blocks(AB)
    theory.find_equivalents(AB)
    theory.compute_cross_mrcas(AB)
    S = span(AB)
    print(newick.compose_newick(checklist.preorder_rows(S)))
  testit("a", "a")              # A + B
  testit("(c,d)a", "(c,e)b")
  testit("((a,b)e,c)d", "(a,(b,c)f)D")

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    """)
  parser.add_argument('--A', help="the A checklist, as path name or -",
                      default='-')
  parser.add_argument('--B', help="the B checklist, as path name or -",
                      default='-')
  parser.add_argument('--test', action='store_true', help="run smoke tests")
  args=parser.parse_args()
  if args.test:
    test()
  else:
    a_path = args.A
    b_path = args.B
    assert a_path != b_path
    with util.stdopen(a_path) as a_file:
      with util.stdopen(b_path) as b_file:
        report = span(csv.reader(a_file),
                      csv.reader(b_file))
        util.write_rows(report, sys.stdout)