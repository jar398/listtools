#!/usr/bin/env python3

import sys, argparse, csv
from util import read_rows, windex

def sort_csv(inport, key_columns, outport):
  reader = csv.reader(inport)
  header = next(reader)
  key_positions = [windex(header, pk_col) for pk_col in key_columns.split(",")]
  print("## Sort key positions: %s" % (key_positions,), file=sys.stderr)

  def sort_key(row):
    return tuple(row[pk_pos] for pk_pos in key_positions)

  rows = read_rows(reader, key=sort_key)

  writer = csv.writer(outport)
  writer.writerow(header)
  for row in sorted(rows.values(), key=sort_key):
    assert len(row) == len(header)
    writer.writerow(row)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    TBD
    """)
  parser.add_argument('--key',
                      help='comma-separated names of columns to sort by')
  # List of fields stored in database or graphdb should be an arg.
  args=parser.parse_args()
  sort_csv(sys.stdin, args.key, sys.stdout)
