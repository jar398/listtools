#!/usr/bin/env python3

import sys, re, csv, argparse, util
from util import windex, log, MISSING

tokenize = re.compile("[^,()]+")

# Given a Newick string, returns an iterable of rows.

def parse_newick(newick):
  yield ("taxonID", "canonicalName",
         "parentNameUsageID", "acceptedNameUsageID",
         "taxonomicStatus", "taxonRank")

  id_counter = [0]
  def gen_id():
    id_counter[0] += 1
    return id_counter[0]

  cursor = [0]

  def parse(superior):
    if cursor[0] >= len(newick):
      return

    id = gen_id()

    if newick[cursor[0]] == '(':
      cursor[0] += 1
      while True:
        yield from parse(id)
        if cursor[0] >= len(newick):
          break
        elif newick[cursor[0]] == ')':
          cursor[0] += 1
          break
        elif newick[cursor[0]] == ',':
          cursor[0] += 1
        else:
          assert False

    name = MISSING
    m = tokenize.match(newick, cursor[0])
    if m: name = m[0]
    cursor[0] += len(name)

    if name.endswith('*'):
      can = name[0:-1]
      parent = MISSING
      accepted = superior
      status = 'synonym'
    else:
      can = name
      parent = superior
      accepted = MISSING
      status = 'accepted'
    # taxonID, canonicalName,
    #           parentNameUsageID, acceptedNameUsageID,
    #           taxonomicStatus
    can = can.strip()
    rank = MISSING
    if can[0].isupper():
      parts = can.split(' ')
      if len(parts) == 2 and parts[1].islower():
        rank = 'species'
    yield (str(id), can, str(parent), str(accepted), status, rank)

  yield from parse(MISSING)
  if False:
    if cursor[0] < len(newick):
      print("! extra stuff after end of newick: %s" % newick[cursor[0]:],
            file=sys.stderr)

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
    With the --newick argument, converts the command line newick
    notation string into CSV, and writes the CSV to standard output.

    Without it, converts CSV read from standard input to newick
    notation, which is written to standard output.
    """)
  parser.add_argument('--newick',
                      default=None,
                      help="tree given in proto-newick syntax")
  args=parser.parse_args()

  if args.newick != None:
    writer = csv.writer(sys.stdout)
    for row in parse_newick(args.newick):
      writer.writerow(row)
  else:
    print(compose_newick(csv.reader(sys.stdin)))
