#!/usr/bin/env python3

# Processes output of plugin.py

import sys, csv
import util

def clean_diff(in_rows):
  header = next(in_rows)
  a_pos = 0
  b_pos = 1
  con_pos = util.windex(header, "Comment")
  cat_pos = util.windex(header, "Category")
  assert con_pos != None
  assert cat_pos != None
  keepers = []
  for row in in_rows:
    a_name = row[a_pos].replace("_", " ")
    b_name = row[b_pos].replace("_", " ")
    cat = row[cat_pos]
    keepers.append((a_name, b_name, cat))
  def get_sort_key(row):
    return row[1] if row[0] == "NA" else row[0]
  keepers.sort(key=get_sort_key)
  yield ("Before_Name", "After_Name", "Category")
  for out_row in keepers:
    yield out_row

writer = csv.writer(sys.stdout)
for out_row in clean_diff(csv.reader(sys.stdin)):
  writer.writerow(out_row)
