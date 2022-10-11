#!/usr/bin/env python3

import sys, re, csv, argparse
from util import MISSING, windex, log

tokenize = re.compile("[^,()]+")

# Given a Newick string, returns an iterable of rows.

def parse_newick(newick):
  rows = []

  rows.append(("taxonID", "canonicalName",
               "parentNameUsageID", "acceptedNameUsageID",
               "taxonomicStatus"))
  counter = [0]

  def gen():
    counter[0] += 1
    return counter[0]

  def parse(n, sup):
    id = gen()

    if n >= len(newick):
      return n

    inferiors = []

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
      if name.endswith('*'):
        can = name[0:-1]
        parent = MISSING
        accepted = sup
        status = 'synonym'
      else:
        can = name
        parent = sup
        accepted = MISSING
        status = 'accepted'
      n += len(name)
    # taxonID, canonicalName,
    #           parentNameUsageID, acceptedNameUsageID,
    #           taxonomicStatus
    rows.append((str(id), can.strip(), str(parent), str(accepted), status))

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
  accepted_pos = windex(header, "acceptedNameUsageID")
  taxstat_pos = windex(header, "taxonomicStatus")
  name_pos = windex(header, "canonicalName")
  assert key_pos != None
  assert parent_pos != None
  rows_dict = {row[key_pos] : row for row in rows}
  inferiors = {}
  roots = []
  # print("newick: %s rows" % len(rows_dict), file=sys.stderr)
  for (key, row) in rows_dict.items():
    accepted_key = row[accepted_pos] if accepted_pos != None else MISSING
    parent_key = row[parent_pos]
    ship = (not accepted_key, key)    # True for children, False for synonyms
    sup_key = accepted_key or parent_key    # '' is falsish
    if sup_key:
      if sup_key in inferiors:
        inferiors[sup_key].append(ship)
      else:
        inferiors[sup_key] = [ship]
    else:
      roots.append(ship)
  print("newick: %s roots, %s superiors" % (len(roots), len(inferiors)),
        file=sys.stderr)
  def suppressp(ship):
    (_, key) = ship
    row = rows_dict[key]
    return (row[taxstat_pos] == "equivalent"
            if taxstat_pos != None
            else False)
  def traverse(ship):
    (accepted, key) = ship
    row = rows_dict[key]
    name = row[name_pos]
    if not accepted: name = name + '*'
    if key in inferiors:
      childs = traverse_seq(inferiors[key])
      if childs == '':
        return name
      else:
        return "(%s)%s" % (childs, name)
    else:
      return name
  def traverse_seq(seq):
    newicks = (traverse(child_ship)
               for child_ship in seq
               if not suppressp(child_ship))
    return ",".join(sorted(newicks))
  return traverse_seq(roots)

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
