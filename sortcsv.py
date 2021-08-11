#!/usr/bin/env python3

import sys, argparse, csv
from util import read_csv, windex

def sort_csv(inport, pk_col, outport):
  (header, rows) = read_csv(inport, pk_col)
  pk_pos = windex(header, pk_col)
  assert pk_pos != None
  writer = csv.writer(outport)
  writer.writerow(header)
  for row in sorted(rows.values(), key=lambda row: row[pk_pos]):
    assert len(row) == len(header)
    writer.writerow(row)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    TBD
    """)
  parser.add_argument('--pk',
                      help='name of column to sort by')
  # List of fields stored in database or graphdb should be an arg.
  args=parser.parse_args()
  sort_csv(sys.stdin, args.pk, sys.stdout)
