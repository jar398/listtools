#!/usr/bin/env python3

import sys, csv, argparse
import util, property as prop
import checklist, workspace, match_records
import theory, span

from util import windex, MISSING
from property import mep_get, mep_set
from rcc5 import *
from checklist import *
from workspace import *
from theory import is_accepted, get_accepted

def demo(A_iter, B_iter):
  A = rows_to_checklist(A_iter, {"name": "A"})  # meta
  B = rows_to_checklist(B_iter, {"name": "B"})  # meta
  AB = workspace.make_workspace(A, B, {"name": "AB"})
  matches_iter = match_records.match_records(checklist_to_rows(A),
                                             checklist_to_rows(B))
  theory.load_matches(matches_iter, AB)
  theory.analyze_tipwards(AB)                # also find tipes
  theory.compute_blocks(AB)
  theory.find_equivalents(AB)
  theory.compute_cross_mrcas(AB)
  span.span(AB)
  return generate_report(AB)

def generate_report(AB):
  # Too difficult to do a generator
  report = []

  report.append(("A id", "A name", "rcc5", "B name", "B id"))

  def do_row(x, y):
    report.append((get_primary_key(x),
                   get_canonical(x, None) or get_scientific(x),
                   rcc5_symbol(theory.cross_relationship(AB, x, y).relationship),
                   get_canonical(y, None) or get_scientific(y),
                   get_primary_key(y)))

  for z in checklist.preorder_records(AB):

    assert len(get_source(z).meta['name']) == 2

    # Eight cases, of which 3 are to be ignored.

    #   yy    y     x     xx
    # 1       acc = acc
    # 2       acc = syn < acc
    # 4 acc > syn = acc
    # 5 acc > syn = syn < acc   IGNORE

    # 3       acc.      < acc
    # 6 acc > syn.      < acc   IGNORE

    # 7 acc >       acc
    # 8 acc >       syn.        IGNORE
    # 9     >                   CAN'T HAPPEN

    def when_A(x):
      # cases 7-8
      if is_accepted(x):
        # case 7, acc >       acc.
        do_row(x, theory.cross_superior(AB, x))

    def when_B(y):
      # cases 1-6
      rel = theory.get_equivalent(z, None)
      if rel:
        #       acc = acc
        #       acc = syn < acc
        # acc > syn = acc
        x = get_outject(rel.record)
        if is_accepted(y):
          # cases 1-2
          if is_accepted(x):
            do_row(x, y)
          else:
            do_row(get_accepted(x), y)
        elif is_accepted(x):
          # case 4, acc.  ... < acc
          do_row(x, get_accepted(y))

      elif is_accepted(y):
        # case 3, acc. ... < acc
        do_row(theory.cross_superior(AB, y), y)

    AB.case(z, when_A, when_B)

  return report

# -----------------------------------------------------------------------------

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    standard input = the merged checklist we're reporting on
    """)
  parser.add_argument('--A', help="the A checklist, as path name or -",
                      default='-')
  parser.add_argument('--B', help="the B checklist, as path name or -",
                      default='-')
  args=parser.parse_args()
  a_path = args.A
  b_path = args.B
  assert a_path != b_path
  with util.stdopen(a_path) as a_file:
    with util.stdopen(b_path) as b_file:
      report = demo(csv.reader(a_file),
                    csv.reader(b_file))
      util.write_rows(report, sys.stdout)
