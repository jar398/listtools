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

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('A')
  parser.add_argument('B')
  parser.add_argument('--expect', default=None)
  args=parser.parse_args()
  n3 = test(args.A, args.B)
  if args.expect:
    if n3 != args.expect:
      print("A+B = %s != %s" ^ (n3, args.expect))
      assert False
  else:
    print("A+B = %s" % n3)
