#!/usr/bin/env python3

import sys, os, csv, argparse, pathlib
from util import windex

mode_column = "mode"

def scatter(inport, dest):

  import pathlib
  pathlib.Path(dest).mkdir(parents=True, exist_ok=True)

  reader = csv.reader(inport)
  header = next(reader)

  files = {}
  writers = {}
  counts = {}

  # open three files for output: remove, update, add

  mode_pos = windex(header, mode_column)
  subheader = [x for x in header]
  del subheader[mode_pos]

  for row in reader:
    mode = row[mode_pos]
    del row[mode_pos]
    if not mode in writers:
      fname = os.path.join(dest, mode + ".csv")
      file = open(fname + ".new", "w")
      files[mode] = file
      writers[mode] = csv.writer(file)
      writers[mode].writerow(subheader)
      counts[mode] = 0
    writer = writers[mode]
    writer.writerow(row)
    counts[mode] += 1
  for file in files.values():
    file.close()
  for mode in counts.keys():
    fname = os.path.join(dest, mode + ".csv")
    n = fname + ".new"
    print("%s -> %s (%s rows)" % (n, fname, counts[mode]))
    os.replace(n, fname)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    TBD
    """)
  parser.add_argument('--dest',
                      help='name of directory to hold the various outputs')
  args=parser.parse_args()
  scatter(sys.stdin, args.dest)

