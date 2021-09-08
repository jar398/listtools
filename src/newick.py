#!/usr/bin/env python3

import sys, re, csv, argparse
from util import MISSING, windex

tokenize = re.compile("[^,()]+")

# Returns a generator yielding rows

def parse_newick(newick):
  rows = []

  rows.append(("taxonID", "canonicalName", "parentNameUsageID",))
  counter = [0]

  def gen():
    counter[0] += 1
    return counter[0]

  def parse(n, parent):
    id = gen()

    if n >= len(newick):
      print("END")
      return (None, n)

    children = []

    if newick[n] == '(':
      n += 1
      while True:
        n = parse(n, id)
        if n >= len(newick):
          break
        elif newick[n] == ')':
          n += 1
          break
        elif newick[n] == ',':
          n += 1
        else:
          assert False

    name = MISSING
    m = tokenize.match(newick, n)
    if m:
      name = m[0]
      n += len(name)
    rows.append((str(id), name.strip(), str(parent),))

    return n
  parse(0, MISSING)
  return (row for row in rows)

def generate_newick(rows):
  header = next(rows)
  key_pos = windex(header, "taxonID")
  parent_pos = windex(header, "parentNameUsageID")
  accepted_pos = windex(header, "acceptedNameUsageID")
  name_pos = windex(header, "canonicalName")
  assert key_pos != None
  assert parent_pos != None
  rows_dict = {row[key_pos] : row for row in rows}
  children = {}
  roots = []
  for (key, row) in rows_dict.items():
    parent_id = row[parent_pos]
    if parent_id != MISSING:
      if key in children:
        children[parent_id].append(key)
      else:
        children[parent_id] = [key]
    else:
      roots.append(key)
  def traverse(key):
    row = rows_dict[key]
    name = row[name_pos]
    if key in children:
      newicks = (traverse(child_key) for child_key in children[key])
      childs = ",".join(sorted(newicks))
      return "(%s)%s" % (childs, name)
    else:
      return name
  return ".".join(map(lambda key:traverse(key), roots))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    TBD
    """)
  parser.add_argument('--newick',
                      default=None,
                      help="tree given in proto-newick syntax")
  args=parser.parse_args()

  if args.newick != None:
    gen = parse_newick(args.newick)
    writer = csv.writer(sys.stdout)
    for row in gen: writer.writerow(row)
  else:
    print(generate_newick(csv.reader(sys.stdin)))
