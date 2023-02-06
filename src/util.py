
# Miscellaneous things used by more than one listtools module

import sys, io, argparse, csv, hashlib

MISSING = ''  # property.MISSING

def read_csv(inport, pk_col=None, key=None):
  return ingest_csv(csv.reader(inport), pk_col, key)

# Iterator -> (header, dict)

def ingest_csv(reader, pk_col=None, key=None):
  header = next(reader)
  assert isinstance(header, list) or isinstance(header, tuple)
  assert pk_col == None or key == None
  assert pk_col != None or key != None
  if key:
    keyfun = key
  else:
    if pk_col == None: pk_col = "taxonID"
    pk_pos = windex(header, pk_col)
    if pk_pos == None:
      print("Column %s not found in %s" % (pk_col, header),
            file=sys.stderr)
      assert False
    keyfun = lambda row: row[pk_pos]
  return (header, read_rows(reader, key=keyfun))

# Read rows and index by specified primary key

def read_rows(reader, key=lambda row: row[0]):
  all_rows = {}
  for row in reader:
    assert isinstance(row[0], str)
    ky = key(row)
    if ky in all_rows:
      print("Multiple records for key %s, e.g. %s" % (ky, row),
            file=sys.stderr)
      assert False
    all_rows[ky] = row
  # print("# read_rows: read %s non-header rows" % len(all_rows), file=sys.stderr)
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

def write_rows(gen, outfile):
  writer = csv.writer(outfile)
  for row in gen:
    writer.writerow(row)

# Random.  Simplifies processing of command line args that are
# filename or - (for stdin); lets you use 'with' without having to
# discriminate between the two argument forms.

class Stdin:
  def __enter__(self): return sys.stdin
  def __exit__(self, exc_type, exc_val, exc_tb): return

def stdopen(x):
  if x == '-' or not x:
    return Stdin()
  else:
    return open(x)

# 

log_allowance = 50

def log(mess):
  global log_allowance
  if log_allowance > 0:
    print(mess, file=sys.stderr)
  elif log_allowance == 0:
    print("*** logging truncated ***", file=sys.stderr)
  log_allowance -= 1
