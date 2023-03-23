#!/usr/bin/env python3

import sys, csv
from util import windex, MISSING

reader = csv.reader(sys.stdin)
header = next(reader)
sci_pos = windex(header, "scientificName")
can_pos = windex(header, "canonicalName")
status_pos = windex(header, "nomenclaturalStatus")
print("# Columns: sci %s can %s status %s" % (sci_pos, can_pos, status_pos),
      file=sys.stderr)
sci_count = 0
for row in reader:
  assert len(row) == len(header), (len(row), len(header), row, header)
  name = MISSING
  if sci_pos != None:
    name = row[sci_pos]
    if name != MISSING:
      sci_count += 1
  if name == MISSING and can_pos != None:
    name = row[can_pos]
  if status_pos != None and 'common' in row[status_pos]:
    name = MISSING  #name.lower()   # force Quality = 0
  if name.startswith('? ') or  name == '?':
    name = 'Wildcard' + name[1:]     # Undone in use_gnparse.py
  print(name, file=sys.stdout)
print("# Found %s scientific names" % sci_count, file=sys.stderr)
