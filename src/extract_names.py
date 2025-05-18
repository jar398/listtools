#!/usr/bin/env python3

import sys, csv
from util import windex, MISSING, log
import regex

year_re_string = '\\b([12][0-9]{3})\\b'
year_re = regex.compile(year_re_string)

reader = csv.reader(sys.stdin)
header = next(reader)
sci_pos = windex(header, "scientificName")
auth_pos = windex(header, "scientificNameAuthorship") # COL
can_pos = windex(header, "canonicalName")
status_pos = windex(header, "nomenclaturalStatus")
gen_pos = windex(header, "genericName")
spec_pos = windex(header, "specificEpithet")
infra_pos = windex(header, "infraspecificEpithet")
rank_pos = windex(header, "taxonRank")
log("# Columns: sci %s can %s status %s" % (sci_pos, can_pos, status_pos))
poly_count = 0
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
        poly_count += 1
      elif rank == "subspecies" and infra_pos != None:
        name = "%s %s %s" % (row[gen_pos], row[spec_pos], row[infra_pos])
        poly_count += 1
  # Do we really want to do this?
  if name != MISSING and auth_pos != None:
    y_match = year_re.search(row[sci_pos])
    if not y_match:
      auth = row[auth_pos]
      a_match = year_re.search(auth)
      if a_match and auth and (auth[0].isupper() or auth[1].isupper()):
        name = "%s %s" % (name, auth)
  if status_pos != None and 'common' in row[status_pos]:
    name = MISSING  #name.lower()   # force Quality = 0 for vernaculars
  name = name.replace(' [sic]', '')
  # E.g. Pecten medius Lamarck, 1819 sensu Daniel, 1884
  if ' sensu ' in name:
    name = name[:name.index(' sensu ')]

  # We want gnparse to treat ? as if it were alphabetic
  if name.startswith('? ') or  name == '?':
    name = 'Xyzzy' + name[1:]     # Undone in use_gnparse.py
  name = name.replace('?', 'xyzzy')

  print(name, file=sys.stdout)
if poly_count > 0:
  log("# Synthesized %s polynomials" % poly_count)
log("# Found %s scientific names" % sci_count)
