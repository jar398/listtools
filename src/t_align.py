#!/usr/bin/env python3

from align import *
import sys, argparse
import newick, match_records

def test(n1, n2):
  import newick, match_records
  a_iterator = newick.parse_newick(n1)
  b_iterator = newick.parse_newick(n2)
  rm_sum_iterator = match_records.match_records(a_iterator, b_iterator)

  a_iterator = newick.parse_newick(n1)
  b_iterator = newick.parse_newick(n2)
  sum_iterator = align(a_iterator, b_iterator, rm_sum_iterator)
  rows = [row for row in sum_iterator]
  write_generated((row for row in rows), sys.stdout)
  return newick.generate_newick((row for row in rows))

def run_test(A, B, expect):
  n3 = test(A, B)
  if expect and n3 != expect:
    print("\n**** A+B = %s != %s\n" % (n3, expect))
  return n3

def run_tests():
  run_test("a", "a", "a")
  run_test("(b)a", "(c)a", "(b,c)a")
  run_test("(b)a", "(b)c", "((b)a)c")
  run_test("(b,c)a", "(c,d)a", "(b,c,d)a")
  run_test("((b,c)g,d)a", "(b,c,d)a", "((b,c)g,d)a")
  run_test("(b,c,d)a", "((b,c)g,d)a", "((b,c)g,d)a")
  run_test("((b,c)g,d,e)a", "(b,c,(d,e)h)a", "((b,c)g,(d,e)h)a")
  # src/t_align.py "((b,c)g,d)a" "(b,(c,d)h)a" ...   LOSE

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('A', nargs='?', default=None)
  parser.add_argument('B', nargs='?', default=None)
  parser.add_argument('--expect', default=None)
  args=parser.parse_args()

  if args.A and args.B:
    n3 = run_test(args.A, args.B, args.expect)
    print("A+B = %s" % n3, file=sys.stderr)
  else:
    run_tests()
