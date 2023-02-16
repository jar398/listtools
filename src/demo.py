#!/usr/bin/env python3

import sys, csv, argparse
import util
import workspace
import align

from rcc5 import *
from checklist import *
from workspace import *

# Returns an Iterable of rows

# Returns generator of lines (strings)

def generate_eulerx(AB, al):
  yield from generate_eulerx_checklist(AB.A)
  yield from generate_eulerx_checklist(AB.B)
  yield from eulerx_alignment(AB, al)

def eulerx_alignment(AB, al):
  A = AB.A; B = AB.B
  yield ("articulation %s-%s %s-%s" %
         (get_tag(A), get_tag(B),
          checklist_description(A), checklist_description(B)))
  for (v, ship, w, note, comment, forwardp) in al:
    assert note
    x = get_outject(v); y = get_outject(w)
    if not is_top(x) and not is_top(y):
      yield eulerx_articulation(x, ship, y, note)

def eulerx_articulation(x, ship, y, note):
  sym = rcc5_symbol(ship)
  if ok_for_eulerx(ship):
    return "[%s %s %s]" % (get_eulerx_qualified_name(x),
                           sym,
                           get_eulerx_qualified_name(y))
  else:
    return ("#[%s %s %s] %s" %
            (get_eulerx_qualified_name(x),
             sym,
             get_eulerx_qualified_name(y),
             note or ''))

def ok_for_eulerx(ship):
  if False:
    if ship == EQ: return '='
    elif ship == LT: return '<'
    elif ship == GT: return '>'
    elif ship == CONFLICT: return '><'
    elif ship == DISJOINT: return '!'
    else: return False
  return True

# -----------------------------------------------------------------------------

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

