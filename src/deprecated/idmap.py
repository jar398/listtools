#!/usr/bin/env python3

# Filter that joins a hierarchy item id column to a table.
# Assumes CSV input and output.
# This is pretty specific to EOL, where items are called 'pages'.

import sys, csv, argparse

item_id_col = "EOLid"
parent_item_id_col = "parentEOLid"

def apply_mappings(mappings, inport, outport):
  reader = csv.reader(inport)
  header = next(reader)
  usage_id_pos = windex(header, "taxonID")
  parent_usage_id_pos = windex(header, "parentNameUsageID")
  accepted_usage_id_pos = windex(header, "acceptedNameUsageID")
  item_id_pos = windex(header, item_id_col)

  if usage_id_pos == None:
    print("** map: No taxonID in header: %s" % (header,),
          file=sys.stderr)
    assert False
  writer = csv.writer(outport, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL)
  out_header = [field for field in header]  # copy
  if item_id_pos == None:
    out_header.append(item_id_col)
  if parent_usage_id_pos != None:
    out_header.append(parent_item_id_col)
  writer.writerow(out_header)
  map_count = 0
  map_parent_count = 0
  map_accepted_count = 0
  self_parent = 0
  for row in reader:
    did_map = False
    if map_count % 500000 == 0:
      print("# map: loading map, row %s" % map_count, file=sys.stderr)
    usage_id = row[usage_id_pos]
    if not usage_id:
      print("** map: No taxon id in column %s: %s" % (usage_id_pos, row,),
            file=sys.stderr)
      assert False
    item_id = mappings.get(usage_id)
    if item_id:
      did_map = True
      map_count += 1
    # Deal with item id
    if item_id_pos != None:
      if item_id:
        have_item_id = row[item_id_pos]
        if have_item_id:
          if have_item_id != item_id:
            print("map: Item id conflict for usage %s; replacing %s with %s" %
                  (usage_id, have_item_id, item_id),
                  file=sys.stderr)
      elif row[item_id_pos]:
        did_map = True
      row[item_id_pos] = item_id
    else:
      row.append(item_id)
    if parent_usage_id_pos != None:
      parent_usage_id = row[parent_usage_id_pos]
      parent_item_id = mappings.get(parent_usage_id) if parent_usage_id else None
      if parent_item_id and parent_item_id == item_id:
        # Flush self-parent links!
        self_parent += 1
        continue
      if parent_item_id:
        map_parent_count += 1
      row.append(parent_item_id)
    if accepted_usage_id_pos != None and not item_id:
      accepted_usage_id = row[accepted_usage_id_pos]
      if accepted_usage_id and accepted_usage_id != usage_id:
        accepted_item_id = mappings.get(accepted_usage_id)
        if accepted_item_id:
          did_map = True
          map_accepted_count += 1
          row[item_id_pos] = accepted_item_id
    assert len(row) == len(out_header)
    if not did_map:
      print("** map: No mapping in row %s" % (row,),
            file=sys.stderr)
      continue
    writer.writerow(row)
  print("map: Mapped %s taxon ids, %s parents, %s accepteds" %
        (map_count, map_parent_count, map_accepted_count),
        file=sys.stderr)
  if self_parent > 0:
    print("map: Suppressed %s self-parent rows" % self_parent, file=sys.stderr)

def read_mappings(mapfile):
  if not mapfile: return None
  mappings = {}
  with open(mapfile, "r") as infile:
    reader = csv.reader(infile)
    next(reader)
    for [usage_id, item_id] in reader:
      mappings[usage_id] = item_id
  print("map: %s mappings" % len(mappings),
        file=sys.stderr)
  return mappings

def windex(header, fieldname):
  if fieldname in header:
    return header.index(fieldname)
  else:
    return None

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Add column for item id with contents determined by mapping each taxonID.
    Similarly for parentNameUsageID and acceptedNameUsageID, if present.
    """)
  parser.add_argument('--mapping',
                      help='name of file where taxonID to item id mapping is stored')
  args=parser.parse_args()
  apply_mappings(read_mappings(args.mapping), sys.stdin, sys.stdout)
