#!/usr/bin/env python3

# Remove some columns

import sys, csv, argparse

def project(keep, drop, inport, outport):
  reader = csv.reader(inport)
  header = next(reader)
  keepers = header
  if keep:
    keepers = keep.split(",")
  if drop:
    droppers = drop.split(",")
    print("# project: Flushing %s" % (droppers,), file=sys.stderr)
    keepers = [column for column in keepers if not (column in droppers)]
  losers = []
  for keeper in keepers:
    if not (keeper in header):
      print("-- column %s not in header %s" % (keeper, header,),
            file=sys.stderr)
      losers.append(keepers.index(keeper))
  losers.reverse()
  for loser in losers:
    del keepers[loser]
  keep_positions = [header.index(keeper) for keeper in keepers]
  print("# project: Keeping %s" % (keepers,), file=sys.stderr)
  print("# project: Keeping %s" % (keep_positions,), file=sys.stderr)
  writer = csv.writer(outport, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL)
  writer.writerow(keepers)
  for row in reader:
    assert len(row) == len(header)
    writer.writerow([row[position] for position in keep_positions])

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
