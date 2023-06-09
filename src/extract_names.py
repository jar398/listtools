#!/usr/bin/env python3

import sys, csv
from util import windex, MISSING

reader = csv.reader(sys.stdin)
header = next(reader)
sci_pos = windex(header, "scientificName")
can_pos = windex(header, "canonicalName")
status_pos = windex(header, "nomenclaturalStatus")
gen_pos = windex(header, "genericName")
spec_pos = windex(header, "specificEpithet")
infra_pos = windex(header, "infraspecificEpithet")
rank_pos = windex(header, "taxonRank")
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
  if name == MISSING:
    if can_pos != None and row[can_pos] != MISSING:
      name = row[can_pos]
    elif rank_pos != None and gen_pos != None and spec_pos != None:
      # genericName,specificEpithet  if rank is species
      # genericName,specificEpithet,infraspecificEpithet  if rank is subspecies
      rank = row[rank_pos]
      if rank == "species":
        name = "%s %s" % (row[gen_pos], row[spec_pos])
      elif rank == "subspecies":
        name = "%s %s %s" % (row[gen_pos], row[spec_pos], row[infra_pos])
  if status_pos != None and 'common' in row[status_pos]:
    name = MISSING  #name.lower()   # force Quality = 0
  if name.startswith('? ') or  name == '?':
    name = 'Wildcard' + name[1:]     # Undone in use_gnparse.py
  print(name, file=sys.stdout)
print("# Found %s scientific names" % sci_count, file=sys.stderr)
