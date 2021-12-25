#!/usr/bin/env python3

import sys, csv, argparse
import util, property as prop
import workspace, merge
from util import windex, MISSING
from property import mep_get, mep_set

from rcc5 import *

from checklist import primary_key_prop, get_primary_key, \
  canonical_prop, get_canonical, \
  scientific_prop, get_scientific, \
  get_superior, \
  get_children, \
  get_canonical, get_rank

id_a_prop = prop.getter(prop.get_property("taxonID_A"))
id_b_prop = prop.getter(prop.get_property("taxonID_B"))

previous_prop = prop.get_property("previous")
get_record_match = prop.getter(previous_prop)
set_record_match = prop.setter(previous_prop)

change_prop = prop.get_property("change")
get_change = prop.getter(change_prop)

def report(a_iter, b_iter, merged_iter, diff_mode):
  #B = workspace.rows_to_checklist(b_iter, {'name': 'B'})
  #A = workspace.rows_to_checklist(a_iter, {'name': 'A'})
  AB = workspace.rows_to_checklist(merged_iter, {'name': 'AB'})
  return generate_report(AB, diff_mode)

# Full mode shows every concept in the sum.
# Diff mode only shows changed/new/removed concepts.

def generate_report(AB, diff_mode):
  yield ("A name", "B name", "rank", "comment", "remarks")
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
    m = get_record_match(u, None)

    # WORK IN PROGRESS
    if False and y and m:
      change2 = rcc5_symbols[related_how(m, y)[0]]
    else:
      change2 = None

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
        # New taxon concept with prior record match
        if change == '<':
          comment = "widened"
        elif change == '>':
          comment = "narrowed"
        elif change == '><':
          comment = "changed incomparably"
        elif change == '?':
          # Accepted with sibling synonym or within-cluster peer.
          # No way to distinguish.
          comment = "promoted to accepted"
        else:
          # = or !
          comment = "merge failure"
        if x: comment = "(%s)" % comment
        tick(comment)

    if comment or not diff_mode:
      rank = get_rank(y or x, MISSING)
      noise = noises.get(rank, ". . . . .")
      bx = get_blurb(x)
      if x and get_accepted(x, None): bx = bx + "*"
      by = get_blurb(y)
      if y and get_accepted(y, None): by = by + "*"
      if not diff_mode:
        bx = bx + " " + noise
        by = noise + " " + by
      add_remark(u, comment)
      yield [bx, by, rank, comment, get_remarks(u, MISSING)]
    for c in get_children(u, []):
      for row in traverse(c): yield row
    for s in get_synonyms(u, []):
      for row in traverse(s): yield row
  for row in traverse(AB.top): yield row
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
    blurb = get_canonical(z, None) or get_scientific(z, None) or get_primary_key(z)
    return blurb
  else:
    return MISSING

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    B hierarchy is read from stdin
    """)
  parser.add_argument('--source', help="A hierarchy")
  parser.add_argument('--merged', help="merged checklist")
  parser.add_argument('--diff', action='store_true',
                      help="difference mode")
  parser.add_argument('--full', action='store_true',
                      help="full report")
  args=parser.parse_args()

  b_file = sys.stdin
  a_path = args.source
  merged_path = args.merged
  diff_mode = not args.full

  with open(a_path) as a_file:
    with open(merged_path) as merged_file:
      rep = report(csv.reader(a_file),
                   csv.reader(b_file),
                   csv.reader(merged_file),
                   diff_mode)
      util.write_generated(rep, sys.stdout)
