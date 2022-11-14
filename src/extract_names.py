#!/usr/bin/env python3

import sys, csv
from util import windex, MISSING

reader = csv.reader(sys.stdin)
header = next(reader)
sci_col = windex(header, "scientificName")
can_col = windex(header, "canonicalName")
status_col = windex(header, "nomenclaturalStatus")
print("# Columns: sci %s can %s status %s" % (sci_col, can_col, status_col),
      file=sys.stderr)
sci_count = 0
for row in reader:
  assert len(row) == len(header), (len(row), len(header), row, header)
  name = MISSING
  if sci_col != None:
    name = row[sci_col]
    if name != MISSING:
      sci_count += 1
  if name == MISSING and can_col != None:
    name = row[can_col]
  if status_col != None and 'common' in row[status_col]:
    name = MISSING  #name.lower()   # force Quality = 0
  print(name.strip(), file=sys.stdout)
print("# Found %s scientific names" % sci_count, file=sys.stderr)
