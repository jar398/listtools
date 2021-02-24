#!/bin/env python3

import sys, csv

def project(specified_columns, csvp, inport, outport):
  if specified_columns.startswith("+"):
    specified_columns = specified_columns[1:]
    keeping = True
  elif specified_columns.startswith("-"):
    specified_columns = specified_columns[1:]
    keeping = False
  else:
    print("** project: List of columns should start with + or -: %s" %
          specified_columns, file=sys.stderr)
    assert False
  columns = specified_columns.split(",")
  (d, q, g) = csv_parameters(csvp)
  reader = csv.reader(inport, delimiter=d, quotechar=q, quoting=g)
  header = next(reader)
  if keeping:
    print("# project: Keeping %s" % (columns,), file=sys.stderr)
    keepers = columns
  else:
    print("# project: Flushing %s" % (columns,), file=sys.stderr)
    keepers = [column for column in header if not (column in columns)]
  positions = [header.index(keeper) for keeper in keepers]
  print("# project: Keeping %s" % (positions,), file=sys.stderr)
  writer = csv.writer(outport, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL)
  writer.writerow(keepers)
  for row in reader:
    assert len(row) == len(header)
    writer.writerow([row[position] for position in positions])

def csv_parameters(path):
  if ".csv" in path:
    return (",", '"', csv.QUOTE_MINIMAL)
  else:
    return ("\t", "\a", csv.QUOTE_NONE)

if __name__ == '__main__':
  specified_columns = sys.argv[1]
  csvp = sys.argv[2] if len(sys.argv) > 2 else "stdin.csv"
  project(specified_columns, csvp, sys.stdin, sys.stdout)
