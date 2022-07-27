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
  def traverse(AB, B_priority):
    AB.top = AB.in_right(AB.B.top)

    for y in checklist.preorder_records(AB.B):    # starts with top
      z = AB.in_right(y)
      assert get_superior(z, None) == None

      # Bury A records that have B equivalents
      equ = theory.get_equivalent(z, None)
      if equ and not B_priority:
        w = equ.record
        log("burying %s under %s" % (blurb(z), blurb(w)))
        sup = relation(SYNONYM, w, "equivalent")

      else:
        # One possible parent of z
        qy_rel = get_superior(y, None)

        if qy_rel:
          q = AB.in_right(qy_rel.record)

          # Another possible parent of z
          p_rel = theory.cross_superior(AB, z)
          assert p_rel
          p = p_rel.record

          # See which one is best
          pq_rel = theory.cross_relationship(AB, p, q)
          if pq_rel.relationship == EQ:
            # Best superior is p = q
            sup = relation(p_rel.relationship, p, "equivalent parents")
          elif pq_rel.relationship == LT:
            # Override in B, jump to A.
            # Best superior is p < q
            # Wait, shouldn't the relationship come from q's checklist ??
            sup = relation(p_rel.relationship, p, "graft")
          else:
            # Best superior is q > p
            sup = relation(qy_rel.relationship,
                           q,
                           qy_rel.note)
        elif not B_priority:
          # this may be redundant ... ?
          sup = relation(SYNONYM, AB.top, 'spanning tree top rule')
        else:
          continue

      checklist.link_superior(z, sup) # adds z to children/synonyms list

  traverse(AB, True)
  traverse(theory.swap(AB), False)
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
    span(AB)
    print(newick.compose_newick(checklist.preorder_rows(AB)))
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
