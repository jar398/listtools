#!/bin/env python3

# Filter that adds an page id (EOLid) column
# Assumes CSV input and output

import sys, csv

page_id_label = "EOLid"
parent_page_id_label = "parentEOLid"
accepted_page_id_label = "acceptedEOLid"

def apply_mappings(mappings, csvp, inport, outport):
  (d, q, g) = csv_parameters(csvp)
  reader = csv.reader(inport, delimiter=d, quotechar=q, quoting=g)
  header = next(reader)
  taxon_id_pos = windex(header, "taxonID")
  parent_taxon_id_pos = windex(header, "parentNameUsageID")
  accepted_taxon_id_pos = windex(header, "acceptedNameUsageID")
  page_id_pos = windex(header, page_id_label)

  if taxon_id_pos == None:
    print("** map: No taxonID in header: %s" % (header,),
          file=sys.stderr)
    assert False
  overwrite = (page_id_pos != None)
  if overwrite:
    print("map: Overwriting %s column with mapped page ids" % page_id_label,
          file=sys.stderr)
  writer = csv.writer(outport, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL)
  out_header = [field for field in header]  # copy
  if overwrite:
    pass
  else:
    out_header.append(page_id_label)
  if parent_taxon_id_pos != None:
    out_header.append(parent_page_id_label)
  if accepted_taxon_id_pos != None:
    out_header.append(accepted_page_id_label)
  writer.writerow(out_header)
  count = 0
  self_parent = 0
  for row in reader:
    did_map = False
    if(count % 500000 == 0): print("# map: Row %s" % count, file=sys.stderr)
    taxon_id = row[taxon_id_pos]
    if not taxon_id:
      print("** map: No taxon id in column %s: %s" % (taxon_id_pos, row,),
            file=sys.stderr)
      assert False
    page_id = mappings.get(taxon_id)
    if page_id:
      did_map = True
    # Deal with page id
    if overwrite:
      if page_id:
        have_page_id = row[page_id_pos]
        if have_page_id:
          have_page_id = int(have_page_id)
          if have_page_id != page_id:
            print("map: Page id conflict for taxon %s; replacing %s by %s" %
                  (taxon_id, have_page_id, page_id),
                  file=sys.stderr)
      row[page_id_pos] = page_id
    else:
      row.append(page_id)
    if parent_taxon_id_pos != None:
      parent_taxon_id = row[parent_taxon_id_pos]
      parent_page_id = mappings.get(parent_taxon_id) if parent_taxon_id else None
      if parent_page_id and parent_page_id == page_id:
        # Flush self-parent nodes!
        self_parent += 1
        continue
      row.append(parent_page_id)
    if accepted_taxon_id_pos != None:
      accepted_taxon_id = row[accepted_taxon_id_pos]
      if accepted_taxon_id and accepted_taxon_id == taxon_id:
        accepted_page_id = None
      else:
        # synonym
        accepted_page_id = mappings.get(accepted_taxon_id)
        if accepted_page_id: did_map = True
      row.append(accepted_page_id)
    assert len(row) == len(out_header)
    if not did_map:
      print("** map: No mapping in row %s" % (row,),
            file=sys.stderr)
      continue
    writer.writerow(row)
    count += 1
  print("map: Emitted %s rows" % count, file=sys.stderr)
  if self_parent > 0:
    print("map: Suppressed %s self-parent rows" % self_parent, file=sys.stderr)

def read_mappings(mapfile):
  if not mapfile: return None
  mappings = {}
  with open(mapfile, "r") as infile:
    reader = csv.reader(infile)
    next(reader)
    for [node_id, page_id] in reader:
      mappings[node_id] = int(page_id)
  print("map: %s mappings" % len(mappings),
        file=sys.stderr)
  return mappings

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
  mapfile = sys.argv[1]
  csvp = sys.argv[2] if len(sys.argv) > 2 else "stdin.csv"
  apply_mappings(read_mappings(mapfile), csvp, sys.stdin, sys.stdout)
