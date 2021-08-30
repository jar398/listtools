#!/usr/bin/env python3

# This is mainly an EOL tool but could be useful as a general
# hierarchy checking tool.

# Input: a usages table, with taxonID as primary key
#   also acceptedNameUsageId, taxonomicStatus, etc.

# Output: a hierarchical items table, with one row per item

import sys, csv, argparse
from util import apply_correspondence, correspondence, windex, MISSING

def hierarchy(keep, infile, outfile, usage_to_item):

  usage_id_col = "taxonID"
  # This choice should be controlled independently of presence of map...
  if usage_to_item:
    item_id_col = "EOLid"
    parent_id_col = "parentEOLid"
  else:
    item_id_col = "taxonID"
    parent_id_col = "parentNameUsageId"

  unmapped = []
  def itemize(usage_id):
    if usage_to_item:
      item_id = usage_to_item.get(usage_id)
      if item_id: return item_id
      item_id = "[%s]" % usage_id
      usage_to_item[usage_id] = item_id
      unmapped.append(usage_id)
      return item_id
    else:
      assert usage_id != MISSING
      item_id = usage_id

  item_rows = {}
  reader = csv.reader(infile)
  header = next(reader)

  usage_pos = windex(header, usage_id_col)
  assert usage_pos != None
  item_id_pos = windex(header, item_id_col)
  parent_usage_pos = windex(header, "parentNameUsageID")
  accepted_usage_pos = windex(header, "acceptedNameUsageID")
  status_pos = windex(header, "taxonomicStatus")
  if keep == "":
    keep_cols = []
  else:
    keep_cols = keep.split(",")
    print("Extra item columns: %s" % (keep_cols,), file=sys.stderr)

  out_header = [item_id_col, parent_id_col, usage_id_col] + keep_cols
  corr = correspondence(header, out_header)
  print("Correspondence: %s" % (corr,), file=sys.stderr)

  item_rows = {} # usage id to output row
  seen_item_ids = {}
  parent = {}    # usage id to parent usage id
  roots = []
  children = {}
  synonyms = 0
  discards = []

  for row in reader:
    usage_id = row[usage_pos]

    # Two ways to test whether a usage is accepted
    au = MISSING
    if accepted_usage_pos != None:
      au = row[accepted_usage_pos]
    indication_1 = (au == MISSING or au == usage_id)
    stat = row[status_pos] if status_pos != None else "accepted"
    indication_2 = (stat == "accepted" or stat == "valid")
    if indication_1 != indication_2:
      print("Conflicting evidence concerning acceptedness: %s %s %s" %
             (usage_id, au, stat),
            file=sys.stderr)

    # Ignore any non-accepted rows
    if indication_1:
      item_id = itemize(usage_id)
      if item_id in seen_item_ids:
        discards.append((usage_id, item_id, seen_item_ids[item_id]))
      else:

        # Work out parents and children
        # This could be done in terms of items not usages
        parent_usage_id = row[parent_usage_pos]
        if parent_usage_id == MISSING:
          roots.append(parent_usage_id)
        else:
          parent[usage_id] = parent_usage_id
          ch = children.get(parent_usage_id)
          if ch:
            ch.append(usage_id)
          else:
            children[parent_usage_id] = [usage_id]

        # Form an output row
        seen_item_ids[item_id] = usage_id
        item_row = apply_correspondence(corr, row)
        if item_row[0] != MISSING and item_row[0] != item_id:
          print("For usage %s, mapping %s will override input file %s" %
                (usage_id, item_id, item_row[0]),
                file=sys.stderr)
        item_row[0] = item_id
        item_rows[usage_id] = item_row    # don't need all of it
    else:
      synonyms += 1

  writer = csv.writer(outfile)
  writer.writerow(out_header)

  # Now fill in parent pointers, generate the output, and get the topology
  for (usage_id, item_row) in item_rows.items():

    parent_item_id = MISSING
    parent_usage_id = parent.get(usage_id)
    if parent_usage_id:
      if parent_usage_id in item_rows:
        parent_item_id = itemize(parent_usage_id)
      else:
        roots.append(usage_id)
    item_row[1] = parent_item_id

    writer.writerow(item_row)

  print("%s items, %s roots, %s items with children, %s non-items, %s unmapped, %s discards" %
        (len(item_rows), len(roots), len(children), synonyms, len(unmapped), len(discards)),
        file=sys.stderr)

  # As a diagnostic service, check that the hierarchy is well-formed.
  seen = {}

  def descend(usage_id):
    seen[usage_id] = True
    for child_id in children.get(usage_id, ()):
      descend(child_id)
  for root in roots:
    descend(root)
  if len(seen) != len(item_rows):
    print("Reached only %s items out of %s by recursive descent" %
          (len(seen), len(item_rows)),
          file=sys.stderr)
    throttle = 0
    for usage_id in item_rows:
      if not usage_id in seen:
        throttle += 1
        if throttle <= 10:
          print("Missed: %s = %s" % (usage_id, itemize(usage_id)),
                file=sys.stderr)

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

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    CSV rows are read from standard input and written to standard output.
    """)
  parser.add_argument('--keep',
                      default="",
                      help="a,b,c where a,b,c are columns to keep")
  parser.add_argument('--mapping',
                      help='name of file where usage id to item id mapping is stored')
  args=parser.parse_args()
  hierarchy(args.keep, sys.stdin, sys.stdout, read_mappings(args.mapping))
