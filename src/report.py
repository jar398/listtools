#!/usr/bin/env python3

import sys, csv, argparse
import util, property as prop, align
from util import windex, MISSING
from property import mep_get, mep_set

from align import key_prop, get_key, \
  out_a_prop, out_a, \
  out_b_prop, out_b, \
  canonical_prop, get_canonical, \
  scientific_prop, get_scientific, \
  get_remark, set_remark, \
  get_parent, set_parent, \
  get_accepted, set_accepted, \
  get_children, get_synonyms, \
  get_canonical, get_rank, get_year

previous_prop = prop.Property("previous")
get_previous = prop.getter(previous_prop)
set_previous = prop.setter(previous_prop)

change_prop = prop.Property("change")
get_change = prop.getter(change_prop)
set_change = prop.setter(change_prop)

def report(a_iter, b_iter, sum_iter, diff_mode):
  (a_usage_dict, a_roots) = align.load_usages(a_iter)
  (b_usage_dict, b_roots) = align.load_usages(b_iter)
  print("%s A usages, %s B usages" % (len(a_usage_dict), len(b_usage_dict)),
        file=sys.stderr)
  al = load_alignment(sum_iter, a_usage_dict, b_usage_dict)
  return generate_report(al, diff_mode)

# Full mode shows every concept in the sum.
# Diff mode only shows changed/new/removed concepts.

def generate_report(al, diff_mode):
  (_, roots) = al
  yield ("A name", "B name", "rank", "year", "comment", "remark")
  stats = {}
  species_stats = {}
  def traverse(u):
    def tick(stat):
      if stat in stats:
        s = stats[stat]
      else:
        s = [0, 0]
        stats[stat] = s
      s[0] += 1
      rank = ((y and get_rank(y, None)) or (x and get_rank(x, None)))
      if rank == "species":
        s[1] += 1
    x = out_a(u, None)
    y = out_b(u, None)
    change = get_change(u, None)   # relative to record match
    name_changed = (get_blurb(x) != get_blurb(y))

    if x and y:
      if get_accepted(x, None):
        if get_accepted(y, None):
          if name_changed:
            # This doesn't happen
            comment = "renamed synonym"
            tick(comment)
          else:
            comment = "kept synonym"
            tick(comment)
          if diff_mode: comment = None
        else:
          # synonym -> accepted
          comment = "promoted to accepted"
          tick(comment)
      elif get_accepted(y, None):
        # accepted -> synonym
        comment = "lumped into %s" % get_blurb(get_accepted(y))
        tick("demoted to synonym")
      else:
        # accepted -> accepted
        if name_changed:
          comment = "renamed accepted"
          tick(comment)
        else:
          comment = None if diff_mode else "kept"
          tick("kept accepted")
    else:
      if x and get_accepted(x, None):
        comment = None
        tick("dropped synonym")
      elif y and get_accepted(y, None):
        comment = None
        tick("added synonym")
      elif not change:            # who knows what happened
        if y:
          comment = "new/split/renamed"
          tick("added accepted")
        else:
          comment = "deprecated/lumped/renamed"
          tick("dropped accepted")
      else:
        if change == '<':
          comment = "widened"
        elif change == '>':
          comment = "narrowed"
        elif change == '><':
          comment = "changed incomparably"
        elif change == '?':
          # Not completely sure about this
          comment = "promoted to accepted"
        else:
          # = or !
          comment = "alignment failure"
        if x: comment = "(%s)" % comment
        tick(comment)

    remark = get_remark(u, None)
    if comment or not diff_mode:
      year = get_year(y or x, MISSING)
      rank = get_rank(y or x, MISSING)
      noise = noises.get(rank, ". . . . .")
      bx = get_blurb(x)
      if x and get_accepted(x, None): bx = bx + "*"
      by = get_blurb(y)
      if y and get_accepted(y, None): by = by + "*"
      if not diff_mode:
        bx = bx + " " + noise
        by = noise + " " + by
      yield [bx, by, rank, year, comment, remark]
    for c in get_children(u, []):
      for row in traverse(c): yield row
    for s in get_synonyms(u, []):
      for row in traverse(s): yield row
  for root in roots:
    for row in traverse(root): yield row
  for key in sorted(list(stats.keys())):
    s = stats[key]
    print("%7d %7s %s" % (s[0], s[1], key), file=sys.stderr)


noises = {"subspecies": "",
          "species": "_",
          "subgenus": " _",
          "genus": "_ _",
          "subfamily": " _ _",
          "family": "_ _ _",
          "order": "_ _ _ _",
          "class": "_ _ _ _ _",
          }

def get_blurb(z):
  if z:
    blurb = get_canonical(z, None) or get_scientific(z, None) or get_key(z)
    return blurb
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
  remark_pos = windex(header, "remark")
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
      if change != MISSING: set_change(union, change)
    if remark_pos != None:
      remark = row[remark_pos]
      if remark != MISSING: set_remark(union, remark)

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
  roots = align.collect_inferiors(key_to_union.values())

  print("%s taxon concepts in alignment" % len(key_to_union),
        file=sys.stderr)

  align.cache_levels(roots)

  return (key_to_union, roots)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    TBD.  stdin = B hierarchy
    """)
  parser.add_argument('--source', help="A hierarchy")
  parser.add_argument('--alignment', help="alignment")
  parser.add_argument('--diff', action='store_true',
                      help="difference mode")
  parser.add_argument('--full', action='store_true',
                      help="full report")
  args=parser.parse_args()

  b_file = sys.stdin
  a_path = args.source
  sum_path = args.alignment
  diff_mode = not args.full

  with open(a_path) as a_file:
    with open(sum_path) as sum_file:
      rep = report(csv.reader(a_file),
                   csv.reader(b_file),
                   csv.reader(sum_file),
                   diff_mode)
      util.write_generated(rep, sys.stdout)
