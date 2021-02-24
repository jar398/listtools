#!/bin/env python3

# Detect duplicates and sort according to some primary key.

import sys, csv

def prepare(pk_spec, csvp, inport, outport):
  pk_fields = pk_spec.split(",")
  print("# prepare: Primary key fields = %s" % pk_fields, file=sys.stderr)
  (d, q, g) = csv_parameters(csvp)
  reader = csv.reader(inport, delimiter=d, quotechar=q, quoting=g)
  header = next(reader)
  for field in pk_fields:
    if not windex(header, field):
      print("# prepare: Primary key field %s not found in header" % field,
            file=sys.stderr)
      print("# prepare: header = %s" % (header,),
            file=sys.stderr)
  pk_positions = [windex(header, field) for field in pk_fields]
  merged = {}
  conflicts = {}
  merges = 0
  count = 0
  writer = csv.writer(outport)
  writer.writerow(header)
  for row in reader:
    pk = primary_key(row, pk_positions)
    if not any(pk):
      print("# prepare: Bad key, row = %s" % (row,),
            file=sys.stderr)
      assert False
    have_row = merged.get(pk)
    if have_row:
      merged[pk] = merge(have_row, row)
      conflicts[pk] = True
      merges += 1
      if merges < 10:
        print("# prepare: Merge <- %s" % have_row, file=sys.stderr)
    else:
      merged[pk] = row
    if count % 500000 == 0:
      print("# prepare: %s" % count, file=sys.stderr)
    count += 1
  s = sorted(merged.keys())
  print("# prepare: sorted %s rows" % len(s), file=sys.stderr)
  for key in s:
    writer.writerow(merged[key])
  print("prepare: %s rows resulting from merges" % len(conflicts), file=sys.stderr)

def primary_key(row1, pk_positions1):
  # not so good I think.
  key = tuple(mint(row1[pos]) for pos in pk_positions1)
  return key

def mint(val):
  if not val:
    return val
  elif val.isdigit():
    return int(val)
  else:
    return val

def merge(row1, row2):
  return [merge_values(x, y) for (x, y) in zip(row1, row2)]

def merge_values(x, y):
  if x == y: return x
  else: return "?"

def windex(header, fieldname):
  if fieldname in header:
    return header.index(fieldname)
  else:
    return None

def csv_parameters(path):
  if ".csv" in path:
    return (",", '"', csv.QUOTE_MINIMAL)
  else:
    return ("\t", "\a", csv.QUOTE_NONE)

if __name__ == '__main__':
  pk = sys.argv[1]
  csvp = sys.argv[2] if len(sys.argv) > 2 else "stdin.csv"
  prepare(pk, csvp, sys.stdin, sys.stdout)
