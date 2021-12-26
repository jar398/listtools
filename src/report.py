#!/usr/bin/env python3

import sys, csv, argparse
import util, property as prop
import checklist, workspace, merge

from util import windex, MISSING
from property import mep_get, mep_set
from rcc5 import *
from checklist import *
from workspace import *

from checklist import primary_key_prop, get_primary_key, \
  canonical_prop, get_canonical, \
  scientific_prop, get_scientific, \
  get_superior, \
  get_children, \
  get_canonical, get_rank

change_prop = prop.get_property("change")
get_change = prop.getter(change_prop)

#

def report(merged_iter, full_report):
  AB = merge.rows_to_merged(merged_iter, {'name': 'AB'})
  return generate_report(AB, full_report)

# Full mode shows every concept in the sum.
# Diff mode only shows changed/new/removed concepts.

def generate_report(AB, full_report):
  yield ("A name", "B name", "rank", "comment", "remarks")
  stats = {}

  def tick(stat):
    if stat in stats:
      s = stats[stat]
    else:
      s = [0, 0]      # [all, species only]
      stats[stat] = s
    s[0] += 1
    rank = ((y and get_rank(y, None)) or (x and get_rank(x, None)))
    if rank == "species":
      s[1] += 1

  for z in all_records(AB):
    def foo():
      m = get_match(z)
      pass
    e = get_equivalent(z)    # z is a B node; get the A node
    foo(z)
    foo(e)

  def traverse(u):
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
          if full_report: comment = None
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
          comment = None if full_report else "kept"
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

    if comment or not full_report:
      rank = get_rank(y or x, MISSING)
      noise = noises.get(rank, ". . . . .")
      bx = get_blurb(x)
      if x and get_accepted(x, None): bx = bx + "*"
      by = get_blurb(y)
      if y and get_accepted(y, None): by = by + "*"
      if not full_report:
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
    log("%7d %7s %s" % (s[0], s[1], key))


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

# -----------------------------------------------------------------------------

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    standard input = the merged checklist we're reporting on
    """)
  parser.add_argument('--mode',
                      help="diff: difference mode, full: full report")
  args=parser.parse_args()

  full_report = (args.mode != 'full')

  rep = report(csv.reader(sys.stdin), full_report)
  util.write_rows(rep, sys.stdout)
