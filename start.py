#!/bin/env python3

import sys, csv, argparse

def start_csv(filename, outport):
  with open(filename, "r") as inport:
    (d, q, g) = csv_parameters(filename)
    reader = csv.reader(inport, delimiter=d, quotechar=q, quoting=g)
    header = next(reader)
    if len(header) == 1:
      if "," in header[0] or "\t" in header[0]:
        print("** start: Suspicious header", file=sys.stderr)
        print("** start: Header is %s" % (row,), file=sys.stderr)
    writer = csv.writer(outport) # CSV not TSV
    writer.writerow(header)
    count = 0
    trimmed = 0
    for row in reader:
      if len(row) > len(header):
        row = row[0:len(header)]
        trimmed += 1
      if len(row) != len(header):
        print(("** start: Unexpected number of columns: have %s want %s" %
               (len(row), len(header))),
              file=sys.stderr)
        print(("** start: Row is %s" % (row,)), file=sys.stderr)
        assert False
      writer.writerow(row)
      count += 1
  print("start: %s rows, %s columns" % (count, len(header)),
        file=sys.stderr)
  if trimmed > 0:
    # Ignoring extra values is appropriate behavior for DH 0.9.  But
    # elsewhere we might want ragged input to be treated as an error.
    print("start: trimmed extra values from %s rows" % (trimmed,),
          file=sys.stderr)
    
def csv_parameters(path):
  if ".csv" in path:
    return (",", '"', csv.QUOTE_MINIMAL)
  else:
    return ("\t", "\a", csv.QUOTE_NONE)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Normalize to csv and check that number of columns is same in every row.
    CSV rows are written to standard output.
    """)
  parser.add_argument('input', default=None,
                      help='name of input file.  TSV assumed unless name contains ".csv"')
  args=parser.parse_args()
  start_csv(args.input, sys.stdout)
