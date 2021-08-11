
# Things used in more than one plotter python file

import sys, io, argparse, csv

MISSING = ''

def read_csv(inport, pk_col):
  reader = csv.reader(inport)
  header = next(reader)
  pk_pos = windex(header, pk_col)
  if pk_pos == None:
    print("Column %s not found in %s" % (pk_col, header),
          file=sys.stderr)
    assert False
  all_rows = {}
  for row in reader:
    assert len(row) == len(header)
    pk = row[pk_pos]
    if pk == MISSING:
      print("Column %s not set in %s" % (pk_col, row),
            file=sys.stderr)
      assert False
    all_rows[pk] = row
  print("Read %s rows" % len(all_rows), file=sys.stderr)
  return (header, all_rows)

def windex(header, fieldname):
  if fieldname in header:
    return header.index(fieldname)
  else:
    return None

# Returns row with same length as correspondence

# correspondence[j] is position of column that maps to jth column
def correspondence(headera, headerb):
  return (len(headera), [windex(headera, col) for col in headerb])

def precolumn(corr, j):
  (n, v) = corr
  return v[j]

def apply_correspondence(corr, rowa):
  (n, v) = corr
  if len(rowa) != n:
    print("Incorrect length input: apply %s\n to %s" % (corr, rowa,),
          file=sys.stderr)
    assert False
  def m(j):
    i = v[j]
    return rowa[i] if i != None else MISSING
  return [m(j) for j in range(0, len(v))]

# test
_corr = correspondence([1,2,3], [3,1])
assert apply_correspondence(_corr, [10,20,30]) == [30,10]

def csv_parameters(path):
  if not path or ".csv" in path:
    return (",", '"', csv.QUOTE_MINIMAL)
  else:
    return ("\t", "\a", csv.QUOTE_NONE)
