#!/usr/bin/env python3

# Processes output of plugin.py

import sys, csv
import util

def mdd_diff(in_rows):
  header = next(in_rows)
  a_pos = util.windex(header, "A concept name")
  b_pos = util.windex(header, "B concept name")
  op_pos = util.windex(header, "operation")
  keepers = []
  for row in in_rows:
    a_name = row[a_pos]
    b_name = row[b_pos]
    op = row[op_pos]
    auth = "year" in op or "author" in op
    if a_name == b_name and not auth:
      continue
    cat = '?'
    if "side" in op:
      if "split" in op: cat = "split"
      elif "overlaps" in op: cat = "overlap"      # ???
    elif "lumped" in op: cat = "lump"
    elif "added" in op: cat = "de novo"
    elif "removed" in op: cat = "removed"
    elif "moved" in op: cat = "genus change"
    elif a_name != b_name:
      cat = "name change" if name_change(a_name, b_name) \
                          else "spelling change"
    elif auth:
      cat = "authorship correction"
    else:
      assert False
    keepers.append((a_name or "NA", b_name or "NA", cat))  #op, 
  def get_sort_key(row):
    return row[1] if row[0] == "NA" else row[0]
  keepers.sort(key=get_sort_key)
  yield ("Before_Name", "After_Name", "Category")  #"Comment", 
  for out_row in keepers:
    yield out_row

def name_change(a_name, b_name):
  m = len(a_name)
  n = len(b_name)
  if abs(m-n) > 3: return True
  z = min(m, n) - 2
  if z <= 2:
    return True
  else:
    return a_name[0:z] != b_name[0:z]

writer = csv.writer(sys.stdout)
for out_row in mdd_diff(csv.reader(sys.stdin)):
  writer.writerow(out_row)
