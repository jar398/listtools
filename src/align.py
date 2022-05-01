#!/usr/bin/env python3

import argparse, sys, os
import util, csv, newick
import checklist, workspace, theory

from rcc5 import *

from checklist import rows_to_checklist, checklist_to_rows, load_matches, \
  get_superior
from newick import parse_newick
from util import log, MISSING
from match_records import match_records


def align(a_iter, b_iter, A_name='A', B_name='B', matches_iter=None):
  A = rows_to_checklist(a_iter, {"name": A_name})  # meta
  B = rows_to_checklist(b_iter, {"name": B_name})  # meta
  if not matches_iter:
    matches_iter = match_records(checklist_to_rows(A), checklist_to_rows(B))
  AB = workspace.make_workspace(A, B, {"name": "AB"})
  load_matches(matches_iter, AB)
  AB.get_cross_mrca = theory.mrca_crosser(AB)
  return align_checklists(AB)

# Traverse either A or B to find equivalences
#   If record match but not topo match: show true articulation
#   If not topo match: (should be A/B symmetric)
#     find witness to >< if any
#     show < if not 'obvious' (commutative diagram)

def align_checklists(AB):
  A = AB.A
  B = AB.B
  theory.analyze_tipwards(AB)                # also find tipes
  theory.compute_blocks(AB)
  theory.ensure_levels(A)           # N.b. NOT levels in AB
  theory.ensure_levels(B)
  alignment = []     # [idA nameA ship nameB idB comment]
  alignment.append(["id in A", "name in B",
                    "rel",
                    "name in B", "id in B",
                    "status", "note"])
  for x in checklist.preorder_records(A): # in A
    xx = AB.in_left(x)                  # in A+B
    xsup = get_superior(x, None)       # in A
    rel = checklist.get_match(xx, None) # Relative
    if rel:
      if rel.record:
        yy = rel.record         # in A+B
        y = theory.get_right_persona(AB, yy) # in B
        ysup = get_superior(y, None) # in B
        #if ysup and ysup.relationship == SYNONYM:
        #  y = ysup.record       # Get accepted
        # Suppress synonym/synonym matches - not interesting
        if not (xsup and xsup.relationship == SYNONYM and
                ysup and ysup.relationship == SYNONYM):
          ship = theory.relationship_per_blocks(AB, xx, yy)
          if ship == COMPARABLE: ship = EQ
          alignment.append([checklist.get_primary_key(x),
                            checklist.blurb(x),
                            rcc5_symbol(ship),
                            checklist.blurb(y),
                            checklist.get_primary_key(y),
                            rel.status,
                            rel.note])
      elif xsup and xsup.relationship != SYNONYM:
        alignment.append([checklist.get_primary_key(x),
                          checklist.blurb(x),
                          MISSING,
                          "ambiguous",
                          MISSING,
                          rel.status,
                          rel.note])
  log("%s articulations" % (len(alignment) - 1))
  return alignment

def test():
  print("ok")
  def testit(m, n):
    log("\n-- Test: %s + %s --" % (m, n))
    return align(newick.parse_newick(m), newick.parse_newick(n))
  testit("a", "a")              # A + B

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Only one of A or B can come from standard input.
    Checklist short names if not provided come from removing the extension
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

    with util.stdopen(a_path) as a_file:
      with util.stdopen(b_path) as b_file:

        def y(matches_iter):
          rows = align(csv.reader(a_file),
                       csv.reader(b_file),
                       A_name = a_name,
                       B_name = b_name,
                       matches_iter = matches_iter)
          writer = csv.writer(sys.stdout)
          for row in rows:
            writer.writerow(row)

        if matches_path:
          with open(matches_path) as matches_file:
            y(csv.reader(matches_file))
        else:
          y(None)
