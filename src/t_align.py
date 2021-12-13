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
  rows = align(a_iterator, b_iterator, rm_sum_iterator)
  #rows = [row for row in sum_iterator]
  #util.write_generated((row for row in rows), sys.stdout)
  return newick.generate_newick((row for row in rows))

def run_test(A, B, expect, noisy):
  print("\nsrc/t_align.py '%s' '%s'" % (A, B))
  print("  or %s + %s" % (newick.generate_newick(newick.parse_newick(A)),
                          newick.generate_newick(newick.parse_newick(B))))
  n3 = test(A, B)
  if expect and n3 != expect:
    print("*** TEST FAILED *** A+B = %s != %s\n" % (n3, expect))
    return 1
  elif noisy:
    print(n3)
  return 0

def run_tests():
  tests = [
    ("a", "a", "a"),
    ("(b)a", "(c)a", "(b,c)a"),
    #("(b)a", "(b)c", "((b)a)c"),
    ("(b,c)a", "(c,d)a", "(b,c,d)a"),
    ("((b,c)g,d)a", "(b,c,d)a", "((b,c)g,d)a"),
    ("(b,c,d)a", "((b,c)g,d)a", "((b,c)g,d)a"),
    ("((b,c)g,d,e)a", "(b,c,(d,e)h)a", "((b,c)g,(d,e)h)a"),
    ("((a,b)c,d)r", "((a,b)e,d)r", "((a,b)e,d)r"),
    # src/t_align.py "((b,c)g,d)a" "(b,(c,d)h)a" ...   LOSE
  ]
  failures = 0
  for (n1, n2, want) in tests:
    failures += run_test(n1, n2, want, False)
  if failures > 0:
    print("\n*** %s tests failed\n" % failures, file=sys.stderr)
  return failures

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('A', nargs='?', default=None)
  parser.add_argument('B', nargs='?', default=None)
  parser.add_argument('--expect', default=None)
  args=parser.parse_args()

  if args.A and args.B:
    n3 = run_test(args.A, args.B, args.expect, True)
    print("A+B = %s" % n3, file=sys.stderr)
  else:
    failures = run_tests()
    sys.exit(0 if failures == 0 else 1)
