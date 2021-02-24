#!/bin/env python3

# Remove some columns

import sys, csv, argparse

def project(keep, drop, inport, outport):
  reader = csv.reader(inport)
  header = next(reader)
  if keep:
    keepers = keep.split(",")
    print("# project: Keeping %s" % (keep,), file=sys.stderr)
  elif drop:
    droppers = drop.split(",")
    print("# project: Flushing %s" % (drop,), file=sys.stderr)
    keepers = [column for column in header if not (column in droppers)]
  positions = [header.index(keeper) for keeper in keepers]
  print("# project: Keeping %s" % (positions,), file=sys.stderr)
  writer = csv.writer(outport, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL)
  writer.writerow(keepers)
  for row in reader:
    assert len(row) == len(header)
    writer.writerow([row[position] for position in positions])

if __name__ == '__main__':
  specified_columns = sys.argv[1]
  parser = argparse.ArgumentParser(description="""
    CSV rows are read from standard input and written to standard output.
    """)
  parser.add_argument('--keep',
                      help="a,b,c where a,b,c are columns to keep (removing all others)")
  parser.add_argument('--drop',
                      help="a,b,c where a,b,c are columns to drop (keeping all others)")
  args=parser.parse_args()
  project(args.keep, args.drop, sys.stdin, sys.stdout)
