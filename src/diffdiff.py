#!/usr/bin/env python3

import sys
import rows
import property as prop

# mdd1.10-mdd1.11-short.csv
# A name,rcc5,B name,category,note,kind,witnesses

# Diff_v1.10-v1.11.csv
# MDDv1.10_Name,MDDv1.11_Name,Comment,Category,Reference

d1_path = 'Diff_v1.10-v1.11.csv'
d2_path = 'work/mdd1.10-mdd1.11-short.csv'

get_A_name = prop.getter(prop.declare_property("A name"))
get_B_name = prop.getter(prop.declare_property("B name"))
get_rcc5_symbol = prop.getter(prop.declare_property("rcc5"))
get_category = prop.getter(prop.declare_property("category"))

# MDDv1.10_Name,MDDv1.11_Name,Comment,Category,Reference
# A name,rcc5,B name,category,note,kind,witnesses

def load_diff(diff_path):
  with rows.open(diff_path) as a_rows:
    row_iterator = iter(a_rows.rows())
    header = next(row_iterator)
    for i in range(0, len(header)):
      if header[i] == "MDDv1.10_Name":
        header[i] = "A name"
      elif header[i] == "MDDv1.11_Name":
        header[i] = "B name"
      elif header[i] == "Category":
        header[i] = "category"
    print(header, file=sys.stderr)
    plan = prop.make_plan_from_header(header)
    table = {}
    for row in row_iterator:
      rec = prop.construct(plan, row)
      b_name = get_B_name(rec).replace('_', ' ')
      if b_name == 'NA': b_name = ''
      elif get_rcc5_symbol(rec, None) == '<!': b_name = ''
      a_name = get_A_name(rec).replace('_', ' ')
      if a_name == 'NA': a_name = ''
      elif get_rcc5_symbol(rec, None) == '>!': a_name = ''
      key = (b_name, a_name)
      table[key] = rec
  return table

d1 = load_diff(d1_path)
d2 = load_diff(d2_path)
print("In: mdd %s listtools %s" % (len(d1), len(d2)), file=sys.stderr)

joined = {}
def enjoin(d, which):
  for (key, rec) in d.items():
    if key in joined:
      j = joined[key]
    else:
      j = [None, None]
      joined[key] = j
    j[which] = rec
enjoin(d1, 0)
enjoin(d2, 1)

k = 0
for (r1, r2) in joined.values():
  if r1 and r2: k += 1

print("Out: in either %s, in both %s" % (len(joined), k), file=sys.stderr)

def normalize_category(cat):
  if cat == 'add': cat = 'de novo'
  elif cat == 'rename': cat = 'spelling change'
  elif cat == 'split (promotion)': cat = 'split'
  elif cat == 'lump (demotion)': cat = 'lump'
  return cat

# with rows.open("foo.csv", "w") as w:
with rows.open("-", "w") as w:
  def generate():
    yield ('A name', 'B name', 'mdd category', 'listtools category')
    for ((b_name, a_name), (r1, r2)) in joined.items():
      cat1 = normalize_category(get_category(r1)) if r1 else '-'
      cat2 = normalize_category(get_category(r2)) if r2 else '-'
      if cat1 != cat2:
        yield (a_name, b_name, cat1, cat2)
  w.write_rows(generate())
