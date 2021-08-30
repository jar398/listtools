#!/usr/bin/env python3

import sys, argparse, csv
from util import windex

def sort_csv(inport, key_columns, outport):
  reader = csv.reader(inport)
  header = next(reader)
  key_positions = [windex(header, pk_col) for pk_col in key_columns.split(",")]
  print("# sort key positions: %s" % (key_positions,), file=sys.stderr)

  def sort_key(row):
    return (tuple(row[pk_pos] for pk_pos in key_positions), row)

  rows = read_rows(reader)
  writer = csv.writer(outport)
  print("# sorting %s rows" % len(rows), file=sys.stderr)
  rows.sort(key=sort_key)
  writer.writerow(header)
  for row in rows:
    assert len(row) == len(header)
    writer.writerow(row)

def read_rows(reader):
  all_rows = []
  for row in reader:
    all_rows.append(row)
  return all_rows

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    TBD
    """)
  parser.add_argument('--key',
                      help='comma-separated names of columns to sort by')
  # List of fields stored in database or graphdb should be an arg.
  args=parser.parse_args()
  sort_csv(sys.stdin, args.key, sys.stdout)
