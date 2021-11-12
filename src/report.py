#!/usr/bin/env python3

import sys, csv, argparse
from util import windex, MISSING
import property as prop, align
from property import mep_get, mep_set

from align import key_prop, get_key, \
  out_a_prop, out_a, \
  out_b_prop, out_b, \
  canonical_prop, get_canonical, \
  scientific_prop, get_scientific, \
  remark_prop, get_remark, \
  get_parent, set_parent, \
  get_accepted, set_accepted, \
  get_children, get_synonyms, \
  get_canonical, get_rank, \
  EQ, LT, GT, DISJOINT, UNCLEAR

previous_prop = prop.Property("previous")
get_previous = prop.getter(previous_prop)
set_previous = prop.setter(previous_prop)

change_prop = prop.Property("change")
get_change = prop.getter(change_prop)
set_change = prop.setter(change_prop)

def report(a_iter, b_iter, sum_iter):
  (a_usage_dict, a_roots) = align.load_usages(a_iter)
  (b_usage_dict, b_roots) = align.load_usages(b_iter)
  print("%s A usages, %s B usages" % (len(a_usage_dict), len(b_usage_dict)),
        file=sys.stderr)
  al = load_alignment(sum_iter, a_usage_dict, b_usage_dict)
  return generate_report(al)

verbose = False

def generate_report(al):
  (_, roots) = al
  yield ("A name", "B name", "rank", "comment")
  def traverse(u):
    x = out_a(u, None)
    y = out_b(u, None)
    change = get_change(u, None)
    if not change:
      if y:
        comment = "new/split/renamed"
      else:
        comment = "deprecated/lumped/renamed"
    elif change == '=':
      if get_blurb(x) == get_blurb(y):
        comment = None
        if verbose: comment = " "
      else:
        comment = "renamed"
    elif change == '<': comment = "widened"
    elif change == '>': comment = "narrowed"
    elif change == '!': comment = "homonym/ambiguous/uncertain"
    elif change == '><': comment = "conflict"
    elif change == '?': comment = "synonym choice"
    else: assert False

    if False:
      m = get_previous(u, None)
      if m == u:
        if get_blurb(x) == get_blurb(y):
          comment = None
          if verbose: comment = " "
        else:
          comment = "renamed"
      elif m:
        how = align.how_related(m, u)
        if how == LT: comment = "widened"
        elif how == GT: comment = "narrowed"
        elif how == EQ: comment = "shouldn't happen"
        elif how == UNCLEAR: comment = "synonym shenanigans"
        else: comment = "reconstituted"
      elif y:
        comment = "new/split/renamed"
      else:
        assert m == None
        comment = "deprecated/lumped/renamed"

    if comment:
      yield [get_blurb(x),
             get_blurb(y),
             get_rank(y or x, MISSING),
             comment]
    for c in get_children(u, []):
      for row in traverse(c): yield row
    if verbose:
      for s in get_synonyms(u, []):
        for row in traverse(s): yield row
  for root in roots:
    for row in traverse(root): yield row

def get_blurb(z):
  if z:
    return get_canonical(z, None) or get_scientific(z, None) or get_key(z)
  else:
    return MISSING

make_union = prop.constructor(key_prop, out_a_prop, out_b_prop)

# Compare load_usages in align.py

def load_alignment(iterator, a_usage_dict, b_usage_dict):
  header = next(iterator)
  # taxonID,taxonID_A,taxonID_B,parentNameUsageID,acceptedNameUsageID,
  #   canonicalName,remark
  key_pos = windex(header, "taxonID")
  remark_pos = windex(header, "remark")
  usage_a_pos = windex(header, "taxonID_A")
  usage_b_pos = windex(header, "taxonID_B")
  previous_pos = windex(header, "previousID")
  change_pos = windex(header, "change")
  parent_pos = windex(header, "parentNameUsageID")
  accepted_pos = windex(header, "acceptedNameUsageID")

  key_to_union = {}   # taxonID -> union
  fixup = []
  accepted_ids = {}

  for row in iterator:
    key = row[key_pos]
    x = a_usage_dict.get(row[usage_a_pos])
    y = b_usage_dict.get(row[usage_b_pos])
    # name, remark
    union = make_union(key, x, y)
    key_to_union[key] = union

    accepted_key = row[accepted_pos] if accepted_pos != None else MISSING
    if accepted_key == key: accepted_key = MISSING
    parent_key = row[parent_pos]
    if accepted_key != MISSING: parent_key = MISSING
    previous_key = row[previous_pos] if previous_pos != None else None
    fixup.append((union, parent_key, accepted_key, previous_key))
    if change_pos != None:
      change = row[change_pos]
      if change != MISSING:
        set_change(union, change)

  for (union, parent_key, accepted_key, previous_key) in fixup:
    if accepted_key != MISSING:
      probe = key_to_union.get(accepted_key)
      if probe: set_accepted(union, probe)
    elif parent_key != MISSING:
      probe = key_to_union.get(parent_key)
      if probe: set_parent(union, probe)
    if previous_key != MISSING:
      probe = key_to_union.get(previous_key)
      if probe: set_previous(union, probe)

  # Collect children so we can say children[x]
  roots = align.collect_children(key_to_union.values())

  print("%s rows in alignment" % len(key_to_union),
        file=sys.stderr)

  align.cache_levels(roots)

  return (key_to_union, roots)

def write_generated(gen, outfile):
  writer = csv.writer(outfile)
  for row in gen:
    writer.writerow(row)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    TBD.  stdin = B hierarchy
    """)
  parser.add_argument('--source', help="A hierarchy")
  parser.add_argument('--alignment', help="alignment")
  args=parser.parse_args()

  b_file = sys.stdin
  a_path = args.source
  sum_path = args.alignment

  with open(a_path) as a_file:
    with open(sum_path) as sum_file:
      rep = report(csv.reader(a_file),
                   csv.reader(b_file),
                   csv.reader(sum_file))
      write_generated(rep, sys.stdout)
