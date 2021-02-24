#!/bin/env python3

# Find the differences between two source tables, and generate a list
# of patch directives that would transform one into the other.

# Records are matched between the sources via their primary keys.

import sys, csv

def prepare_diff_report(pk_spec, path2, inport1, outport):
  pk_fields = pk_spec.split(",")
  print("# diff: Primary key fields = %s" % pk_fields, file=sys.stderr)
  (reader1, header1, pk_positions1) = start(inport1, pk_fields, ".csv")    #foooooo
  print("# diff: Header 1 = %s" % header1, file=sys.stderr)
  with open(path2, "r") as inport2:
    (reader2, header2, pk_positions2) = start(inport2, pk_fields, path2)
    print("# diff: Header 2 = %s" % header2, file=sys.stderr)
    if sorted(header1) != sorted(header2):
      print("# diff: Tables have different columns", file=sys.stderr)
    column_map = [(windex(header1, header2[j]), j) for j in range(len(header2))]
    print("# diff: Column mapping is %s" % column_map, file=sys.stderr)
    writer = csv.writer(outport)
    writer.writerow(["op"] + header2)
    row1 = None
    row2 = None
    added = 0
    removed = 0
    changed = 0
    continued = 0
    count = 0
    while True:
      if row1 == None:
        try:
          row1 = next(reader1)
          if len(row1) != len(header1):
            print("** Row %s of stdin is ragged" % (count,))
          pk1 = primary_key(row1, pk_positions1)
        except StopIteration:
          row1 = False
      if row2 == None:
        try:
          row2 = next(reader2)
          if len(row2) != len(header2):
            print("** Row %s of %s is ragged" % (count, path2,))
          pk2 = primary_key(row2, pk_positions2)
        except StopIteration:
          row2 = False
      assert row1 != None and row2 != None    # should be obvious
      if not row1 and not row2:
        break
      if count % 500000 == 0:
        print("# diff: %s" % count, file=sys.stderr)
      count += 1
      if row1 and (not row2 or pk1 < pk2):
        writer.writerow(["remove"] + row1)
        row1 = None
        removed += 1
      elif row2 and (not row1 or pk2 < pk1):
        writer.writerow(["add"] + row2)
        row2 = None
        added += 1
      elif pk1 == pk2:
        # Which columns?
        d = row_diff(row1, row2, column_map)
        if d:
          writer.writerow([d] + row2)
          changed += 1 
        else:
          continued += 1          
        row1 = row2 = None
      else:
        print("** Fail. Keys are %s %s" % (pk1, pk2), file=sys.stderr)
        assert False

    print("# diff: Added:     %s" % added, file=sys.stderr)
    print("# diff: Removed:   %s" % removed, file=sys.stderr)
    print("# diff: Changed:   %s" % changed, file=sys.stderr)
    print("# diff: Continued: %s" % continued, file=sys.stderr)

def primary_key(row1, pk_positions1):
  return tuple(mint(row1[pos]) for pos in pk_positions1)

def mint(val):
  if val and val.isdigit():
    return int(val)
  else:
    return val

def row_diff(row1, row2, column_map):
  if not isinstance(row1, list):
    print("!!! %s" % row1, file=sys.stderr)
  if not isinstance(row2, list):
    print("!!!! %s" % row2, file=sys.stderr)
  return "".join([value_diff(row1[i], row2[j], j) for (i, j) in column_map])

def value_diff(x1, x2, j):
  if x1 == x2:
    return ""
  elif not x1:
    return "a%s " % j
  elif not x2:
    return "r%s " % j
  else:
    return "c%s " % j  

def continues(row1, row2, column_map):
  # True if values in row2 all come from row1
  for (i, j) in column_map:
    if row2[j] != row1[i]:
      return False

def start(inport, pk_fields, extension):
  (d, q, g) = csv_parameters(extension)
  reader = csv.reader(inport, delimiter=d, quotechar=q, quoting=g)
  header = next(reader)
  for field in pk_fields:
    if windex(header, field) == None:
      print("# diff: Primary key field %s not found in header\n  %s" %
            (field, header),
            file=sys.stderr)
  pk_positions = [windex(header, field) for field in pk_fields]
  return (reader, header, pk_positions)

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
  pk_spec = sys.argv[1]
  path2 = sys.argv[2] if len(sys.argv) > 2 else "stdin.csv"
  prepare_diff_report(pk_spec, path2, sys.stdin, sys.stdout)
