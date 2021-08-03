
# Things used in more than one plotter python file

import sys, io, argparse, csv

MISSING = ''

def read_csv(inport, pk_col):
  reader = csv.reader(inport)
  header = next(reader)
  pk_pos = windex(header, pk_col)
  all_rows = {}
  for row in reader:
    pk = row[pk_pos]
    assert pk != MISSING
    all_rows[pk] = row
  print("Read %s rows" % len(all_rows), file=sys.stderr)
  return (header, all_rows)

def windex(header, fieldname):
  if fieldname in header:
    return header.index(fieldname)
  else:
    return None

# Returns row with same length as correspondence

def map_row(correspondence, rowa):
  def m(j):
    i = correspondence[j]
    return rowa[i] if i != None else MISSING
  return [m(j) for j in range(0, len(correspondence))]

