
# Things used in more than one plotter python file

import sys, io, argparse, csv, hashlib

MISSING = ''

def read_csv(inport, pk_col=None, key=None):
  reader = csv.reader(inport)
  header = next(reader)
  assert pk_col == None or key == None
  assert pk_col != None or key != None
  if pk_col:
    pk_pos = windex(header, pk_col)
    if pk_pos == None:
      print("Column %s not found in %s" % (pk_col, header),
            file=sys.stderr)
      assert False
    keyfun = lambda row: row[pk_pos]
  else:
    keyfun = key
  return (header, read_rows(reader, key=keyfun))

def read_rows(reader, key=lambda row: row):
  all_rows = {}
  for row in reader:
    ky = key(row)
    if ky in all_rows:
      print("Multiple records for key %s, e.g. %s" % (ky, row),
            file=sys.stderr)
      assert False
    all_rows[ky] = row
  print("# read_rows: read %s rows" % len(all_rows), file=sys.stderr)
  return all_rows

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

def stable_hash(lst, without=None):
  if without != None:
    lst = lst + []
    del lst[without]
  return hashlib.sha1("^".join(lst).encode('utf-8')).hexdigest()[0:8]
