#!/usr/bin/env python3

import sys, re, csv, argparse
from util import MISSING, windex

tokenize = re.compile("[^,()]+")

# Given a Newick string, returns an iterable of rows.

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
      return n

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
  n = parse(0, MISSING)
  if n < len(newick):
    print("! extra stuff after end of newick: %s" % newick[n:],
          file=sys.stderr)
  return rows

# Consumes an iterable of rows, and returns a Newick string.

def compose_newick(rows):
  rows = iter(rows)
  header = next(rows)
  key_pos = windex(header, "taxonID")
  parent_pos = windex(header, "parentNameUsageID")
  name_pos = windex(header, "canonicalName")
  assert key_pos != None
  assert parent_pos != None
  rows_dict = {row[key_pos] : row for row in rows}
  children = {}
  roots = []
  # print("newick: %s rows" % len(rows_dict), file=sys.stderr)
  for (key, row) in rows_dict.items():
    parent_key = row[parent_pos]
    if parent_key != MISSING:
      if parent_key in children:
        children[parent_key].append(key)
      else:
        children[parent_key] = [key]
    else:
      roots.append(key)
  # print("newick: %s roots, %s superiors" % (len(roots), len(children)),
  #      file=sys.stderr)
  def traverse(key):
    row = rows_dict[key]
    name = row[name_pos]
    if key in children:
      newicks = (traverse(child_key) for child_key in children[key])
      childs = ",".join(sorted(newicks))
      return "(%s)%s" % (childs, name)
    else:
      return name
  return ",".join(map(traverse, roots))

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
    print(compose_newick(csv.reader(sys.stdin)))
