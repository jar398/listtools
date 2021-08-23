#!/usr/bin/env python3

# Input: a usages table, with public key taxonID, EOLid, parentEOLid,
# and lots of other columns

# Output: a pages table, with public key EOLid, parentEOLid, and
# accepted taxonID
# and maybe landmark_status

import sys, csv, argparse
from util import apply_correspondence, correspondence, windex, MISSING

def pages(keep, infile, outfile):
  pages = {}
  reader = csv.reader(infile)
  header = next(reader)
  page_pos = windex(header, "EOLid")
  usage_pos = windex(header, "taxonID")
  assert usage_pos != None
  parent_usage_pos = windex(header, "parentNameUsageID")
  accepted_usage_pos = windex(header, "acceptedNameUsageID")
  status_pos = windex(header, "taxonomicStatus")
  if keep == "":
    keep_cols = []
  else:
    keep_cols = keep.split(",")
    print("Extra page columns: %s" % (keep_cols,), file=sys.stderr)
  out_header = ["EOLid", "parentEOLid", "taxonID"] + keep_cols
  corr = correspondence(header, out_header)
  print("Correspondence: %s" % (corr,), file=sys.stderr)

  pages = {}     # page id to output row
  mapping = {}   # usage id to page id
  parent = {}    # usage id to parent usage id
  accepted = {}  # usage id to True status "accepted" or "valid"
  unpaged = 0

  for row in reader:
    usage_id = row[usage_pos]
    page_id = row[page_pos]

    mapping[usage_id] = page_id
    parent_usage_id = row[parent_usage_pos]
    if parent_usage_id == usage_id:
      print("Self-parent %s" % parent_usage_id, file=sys.stderr)
      parent_usage_id = MISSING
    elif parent_usage_id != MISSING:
      parent[usage_id] = parent_usage_id
    stat = row[status_pos] if status_pos != None else "accepted"
    acc = (stat == "accepted" or stat == "valid")
    if acc: accepted[usage_id] = acc

    # Decide which usage (input row) is the canonical one for this page
    if page_id == MISSING:
      if acc: unpaged += 1
    else:
      # Compare this usage with previous usage for same page, if any
      prev = pages.get(page_id)
      x = apply_correspondence(corr, row)
      if prev:
        prev_usage_id = prev[2]
        if prev_usage_id in parent and parent_usage_id == MISSING:
          # If we'd lose a parent by replacing prev with x, don't
          pass
        elif prev_usage_id in accepted and not acc:
          # If we'd pass from accepted to synonym, don't
          pass
        else:
          pages[page_id] = x
      else:
        # First usage for this page id
        pages[page_id] = x

  print("%s accepted usages lack a page id" % unpaged, file=sys.stderr)

  # Fill in parent page id for each output row
  for (page_id, x) in pages.items():
    usage_id = x[2]
    parent_usage_id = parent.get(usage_id)
    if parent_usage_id:
      x[1] = mapping.get(parent_usage_id, MISSING)

  (roots, children) = get_topology(pages)

  writer = csv.writer(outfile)
  writer.writerow(out_header)
  count = [0]
  def descend(page_id):
    writer.writerow(pages[page_id])
    count[0] += 1
    for child_id in children.get(page_id, ()):
      descend(child_id)
  for root in roots:
    descend(root)
  if count[0] != len(pages):
    print("Reached only %s pages out of %s by recursive descent" %
          (count[0], len(pages)),
          file=sys.stderr)

def is_accepted_status(status):
  return status == "accepted" or status == "valid"

def get_topology(pages):
  roots = []
  children = {}
  for (page_id, page) in pages.items():
    parent_id = page[1]
    if parent_id == MISSING:
      roots.append(page_id)
    elif not parent_id in pages:
      print("Missing parent page %s for %s" % (parent_id, page_id),
            file=sys.stderr)
      page[1] = MISSING
      roots.append(page_id)
    else:
      ch = children.get(parent_id)
      if ch:
        ch.append(page_id)
      else:
        children[parent_id] = [page_id]
  print("%s pages, %s roots, %s pages with children" %
        (len(pages), len(roots), len(children)),
        file=sys.stderr)
  return (roots, children)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    CSV rows are read from standard input and written to standard output.
    """)
  parser.add_argument('--keep',
                      default="",
                      help="a,b,c where a,b,c are columns to keep")
  args=parser.parse_args()
  pages(args.keep, sys.stdin, sys.stdout)
