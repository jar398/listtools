#!/usr/bin/env python3

import sys, csv, argparse
import util
import workspace
import align

from checklist import *
from workspace import *
from eulerx import generate_eulerx

def test():
  def testit(m, n):
    log("\n-- Test: %s + %s --" % (m, n))
    (report, eulerx, tipwards) = demo(newick.parse_newick(n), "A",
                                      newick.parse_newick(m), "B")
    util.write_rows(report, sys.stdout)
    for line in eulerx:
      print(line, file=sys.stdout)
  testit("a", "a")              # A + B
  #testit("(c,d)a", "(c,e)b")
  #testit("((a,b)e,c)d", "(a,(b,c)f)D")
  #testit("(a,b*)c", "(a)c")
  testit("((b*)a)c", "(a,b)c")     # WRONG WRONG

def demo(A_iter, A_name, B_iter, B_name):
  (al, AB) = align.ingest(A_iter, A_name, B_iter, B_name)
  return (align.generate_alignment_report(al, AB),
          align.generate_short_report(al, AB),
          generate_eulerx(AB, al))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    """)
  parser.add_argument('--A', help="the A checklist, as path name or -",
                      default='-')
  parser.add_argument('--Aname', help="short name of the A checklist",
                      default=None)
  parser.add_argument('--B', help="the B checklist, as path name or -",
                      default='-')
  parser.add_argument('--Bname', help="short name of the B checklist",
                      default=None)
  parser.add_argument('--eulerx', help="where to put the Euler/X version of the alignment",
                      default=None)
  parser.add_argument('--short', help="where to put the diff report",
                      default=None)
  parser.add_argument('--test', action='store_true', help="run smoke tests")
  args=parser.parse_args()
  if args.test:
    test()
  else:
    aname = args.Aname
    bname = args.Bname
    a_path = args.A
    b_path = args.B
    e_path = args.eulerx
    d_path = args.short
    assert a_path != b_path
    with util.stdopen(a_path) as a_file:
      with util.stdopen(b_path) as b_file:
        (al, AB) = align.ingest(csv.reader(a_file),
                                aname,
                                csv.reader(b_file),
                                bname)
        al = list(al)
        report = align.generate_alignment_report(al, AB)
        util.write_rows(report, sys.stdout)
        if d_path:
          short = align.generate_short_report(al, AB)
          with open(d_path, "w") as d_file:
            util.write_rows(short, d_file)
        if e_path:
          eulerx = generate_eulerx(AB, al)
          with open(e_path, "w") as e_file:
            for line in eulerx:
              print(line, file=e_file)

