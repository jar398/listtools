#!/usr/bin/env python3

import sys, csv, argparse
import util, property as prop
import checklist, workspace, match_records
import theory, span, align

from util import windex, MISSING
from property import mep_get, mep_set
from rcc5 import *
from checklist import *
from workspace import *
from theory import is_accepted, get_accepted

def demo(A_iter, B_iter):
  A = rows_to_checklist(A_iter, {"name": "A"})  # meta
  B = rows_to_checklist(B_iter, {"name": "B"})  # meta

  al = list(align.generate_alignment(A, B))
  return (generate_report(A, B, al),
          generate_eulerx(A, B, al))

# Returns an Iterable

def generate_report(A, B, al):
  yield ("kind",
         "A id", "A name", "rcc5", "B name", "B id",
         "action", "comment")
  for (kind, A_id, rcc5, B_id, action, comment) in al:
    r = look_up_record(A, A_id)
    s = look_up_record(B, B_id)
    if r and s:
      yield (kind,
             A_id, get_canonical(r), rcc5, get_canonical(s), B_id,
             action, comment)

# Returns generator of lines (strings)

def generate_eulerx(A, B, al):
  yield from generate_eulerx_checklist(A)
  yield from generate_eulerx_checklist(B)
  yield from eulerx_alignment(A, B, al)

def eulerx_alignment(A, B, al):
  yield ("alignment %s-%s %s-%s" %
         (checklist_tag(A), checklist_tag(B),
          checklist_description(A), checklist_description(B)))
  for (kind, A_id, rcc5, B_id, action, comment) in al:
    r = look_up_record(A, A_id)
    s = look_up_record(B, B_id)
    if r and s:
      yield eulerx_articulation(r, rcc5, s)

def eulerx_articulation(r, rcc5, s):
  re = eulerx_relationship(rcc5)
  if re:
    return "[%s %s %s]" % (get_eulerx_qualified_name(r), re, get_eulerx_qualified_name(s))
  else:
    return "#[%s %s %s]" % (get_eulerx_qualified_name(r), rcc5, get_eulerx_qualified_name(s))


def eulerx_relationship(rcc5):
  if rcc5 == '=': return 'equals'
  elif rcc5 == '<': return 'is_included_in'
  elif rcc5 == '>': return 'includes'
  elif rcc5 == '><': return 'overlaps'
  elif rcc5 == '!': return 'disjoint'
  else: return False

# -----------------------------------------------------------------------------

def test():
  def testit(m, n):
    log("\n-- Test: %s + %s --" % (m, n))
    (report, eulerx) = demo(newick.parse_newick(n),
                            newick.parse_newick(m))
    util.write_rows(report, sys.stdout)
    for line in eulerx:
      print(line, file=sys.stdout)
  testit("a", "a")              # A + B
  testit("(c,d)a", "(c,e)b")
  testit("((a,b)e,c)d", "(a,(b,c)f)D")
  testit("(a,b*)c", "(a)c")
  testit("((b*)a)c", "(a,b)c")     # WRONG WRONG

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    """)
  parser.add_argument('--A', help="the A checklist, as path name or -",
                      default='-')
  parser.add_argument('--B', help="the B checklist, as path name or -",
                      default='-')
  parser.add_argument('--eulerx', help="where to put the Euler/X version of the alignment",
                      default=None)
  parser.add_argument('--test', action='store_true', help="run smoke tests")
  args=parser.parse_args()
  if args.test:
    test()
  else:
    a_path = args.A
    b_path = args.B
    e_path = args.eulerx
    assert a_path != b_path
    with util.stdopen(a_path) as a_file:
      with util.stdopen(b_path) as b_file:
        (report, eulerx) = demo(csv.reader(a_file),
                                csv.reader(b_file))
        util.write_rows(report, sys.stdout)
        if e_path:
          with open(e_path, "w") as e_file:
            for line in eulerx:
              print(line, file=e_file)
