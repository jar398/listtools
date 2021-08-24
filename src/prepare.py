#!/bin/env python3

# Detect duplicates and sort according to some primary key.

import sys, csv, argparse

def prepare(pk_spec, inport, outport):
  pk_fields = pk_spec.split(",")
  print("# prepare: Primary key fields = %s" % pk_fields, file=sys.stderr)
  reader = csv.reader(inport)
  header = next(reader)
  for field in pk_fields:
    if windex(header, field) == None:
      print("# prepare: Primary key field %s not found in header" % field,
            file=sys.stderr)
      print("# prepare: header = %s" % (header,),
            file=sys.stderr)
  pk_positions = [windex(header, field) for field in pk_fields]
  can_pos = windex(header, "canonicalName")
  sci_pos = windex(header, "scientificName")
  merged = {}
  conflicts = {}
  scinames = {}
  ambiguous_scinames = {}
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
    if sci_pos != None:
      name = row[sci_pos]
      if not name:
        name = row[can_pos]
    else:
      name = row[can_pos]
    if name in scinames:
      print("# ambiguous scientific name: %s" % name)
    scinames[name] = row
    if count % 500000 == 0:
      print("# prepare: %s" % count, file=sys.stderr)
    count += 1
  s = sorted(merged.keys())
  print("# prepare: sorted %s rows" % len(s), file=sys.stderr)
  for key in s:
    writer.writerow(merged[key])
  print("prepare: %s rows resulting from merges" % len(conflicts), file=sys.stderr)
  print("prepare: %s ambiguous names" % len(ambiguous_scinames), file=sys.stderr)

# Compare diff.py
def primary_key(row1, pk_positions1):
  # not so good I think.
  return tuple(mint(row1[pos]) for pos in pk_positions1)

def mint(val):
  if val and val.isdigit():
    return (int(val), val)
  else:
    return (1000000000, val)

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
  parser = argparse.ArgumentParser(description="""
    Merge duplicate records, and sort all records.
    Rows merge when they have the same key.  Sort is according to key value.
    CSV rows are read from standard input.  Merged and sorted are written
    to standard output.
    """)
  parser.add_argument('--key',
                      help="names 'a,b,c' for columns that together form the sort key")
  args=parser.parse_args()
  prepare(args.key, sys.stdin, sys.stdout)
