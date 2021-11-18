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
  get_canonical, get_rank

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

# Full style shows every concept in the sum.
# Diff style only shows changed/new/removed concepts.
use_diff_style = True

def generate_report(al):
  (_, roots) = al
  yield ("A name", "B name", "rank", "comment", "remark")
  def traverse(u):
    x = out_a(u, None)
    y = out_b(u, None)
    change = get_change(u, None)
    name_changed = (get_blurb(x) != get_blurb(y))
    if not y:
      comment = "deprecated/lumped/renamed"
      # Never report added synonyms
      if get_accepted(x, None): comment = None
    elif not x:
      comment = "new/split/renamed"
      # Never report lost synonyms
      if get_accepted(y, None): comment = None
    elif not change:
      # Not sure how this can happen!
      comment = "uncertain/ambiguous"
    elif change == '=':
      comment = "="
      # Never report on pure synonyms even if renamed
      if ((not x or get_accepted(x, None)) and
          (not y or get_accepted(y, None))):
        comment = None
      # In diff style, do not report on = pairs with unchanged name
      elif not name_changed and use_diff_style:
        comment = None
      else:
        comment = "renamed"
    elif change == '<':
      if (y and x and
          get_accepted(y, None)):
        if get_accepted(x, None):
          comment = None      # synonym alignments are uninteresting
        else:
          comment = "synonymized under %s" % get_blurb(get_accepted(y))
      elif name_changed:
        comment = "lumped"
      else:
        comment = "widened"
    elif change == '>': comment = "split" if name_changed else "narrowed"
    elif change == '!': comment = "ambiguous/uncertain" if name_changed else "homonym??"
    elif change == '><': comment = "conflict"
    elif change == '?': comment = "synonym choice"
    else: assert False

    remark = get_remark(u, None)
    if comment:
      rank = get_rank(y or x, MISSING)
      noise = noises.get(rank, ". . . . .")
      bx = get_blurb(x)
      if x and get_accepted(x, None): bx = "*" + bx
      by = get_blurb(y)
      if y and get_accepted(y, None): by = "*" + by
      yield [bx + " " + noise,
             noise + " " + by,
             rank, comment, remark]
    for c in get_children(u, []):
      for row in traverse(c): yield row
    for s in get_synonyms(u, []):
      for row in traverse(s): yield row
  for root in roots:
    for row in traverse(root): yield row

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
      if change != MISSING:
        set_change(union, change)
    if remark_pos != None:
      remark = row[remark_pos]
      if remark != MISSING:
        set_remark(union, remark)

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
      util.write_generated(rep, sys.stdout)
