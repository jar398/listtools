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

# Returns an Iterable

def generate_report(AB):
  # Too difficult to do a generator
  report = []

  report.append(("A id", "A name", "rcc5", "B name", "B id", "comment"))

  def do_row(v, w, comment):
    x_info = ('', '')
    y_info = ('', '')
    sym = ''
    if v:
      x = get_outject(v)
      x_info = (get_primary_key(x),
                get_canonical(x, None) or get_scientific(x))
    if w:
      y = get_outject(w)
      y_info = (get_primary_key(y),
                get_canonical(y, None) or get_scientific(y))
    if v and w:
      rel = theory.cross_relationship(AB, v, w)
      sym = rcc5_symbol(rel.relationship)

    report.append((x_info[0],
                   x_info[1],
                   sym,
                   y_info[1],
                   y_info[0],
                   comment))

  def traverse(z):

    assert len(get_source(z).meta['name']) == 2
    if monitor(z): log("traverse")

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
      if get_rank(x, 'species') == 'species':
        # cases 7-8
        if is_acceptable(x):
          # case 7, acc >       acc.
          if monitor(z): log("when_A")
          comment = 'A only'
          rel = get_match(z, None) # in B
          if rel:
            if rel.relationship == EQ:
              comment = 'match has different extension'
            elif rel.record:
              comment = 'ambiguous in A; match in B is %s' % blurb(rel.record)
            else:
              comment = 'ambiguous in B'
          do_row(z, theory.cross_superior(AB.swap(), z), comment)

    def when_B(y):
      if get_rank(y, 'species') == 'species':
        # cases 1-6
        rel = theory.get_equivalent(z, None)
        if rel:
          #       acc = acc
          #       acc = syn < acc
          # acc > syn = acc
          v = rel.record
          x = get_outject(rel.record)
          if is_acceptable(y):
            # cases 1-2
            if is_acceptable(x):
              if monitor(z): log("when_B case 1")
              do_row(v, z, '')
            else:
              if monitor(z): log("when_B case 2")
              do_row(AB.in_left(get_accepted(x)), z,
                     'via A synonym')
          elif is_acceptable(x):
            # case 4, acc.  ... < acc
            if monitor(z): log("when_B case 4")
            do_row(v, AB.in_right(get_accepted(y)),
                   'via B synonym')

        elif is_acceptable(y):
          # case 3, acc. ... < acc
          if monitor(z): log("when_B case 3")
          comment = 'B only'
          rel = get_match(z, None)
          if rel:
            if rel.relationship == EQ:
              comment = 'match has different extension'
            elif rel.record:
              comment = "ambiguous in B; match in A is %s" % blurb(rel.record)
            else:
              comment = 'ambiguous in A'
          do_row(theory.cross_superior(AB, z), z, comment)

    if z != AB.top:
      AB.case(z, when_A, when_B)
    for c in sorted(get_inferiors(z), key=sort_key):
      traverse(c)

  traverse(AB.top)

  return report

def is_acceptable(x):
  return is_accepted(x) and not get_rank(x, None) == 'subspecies'

def sort_key(c):
  return (get_canonical(c, ''),
          get_scientific(c, ''),
          get_primary_key(c, ''))

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
