#!/usr/bin/env python3

import sys, csv, argparse
import util, property as prop
import checklist, workspace, merge
import theory

from util import windex, MISSING
from property import mep_get, mep_set
from rcc5 import *
from checklist import *
from workspace import *

from merge import get_left_id, get_right_id


#

def report(merged_iter, full_report):
  AB = merge.rows_to_merged(merged_iter, {'name': 'AB'})
  merge.resolve_merge_links(AB)
  ensure_inferiors_indexed(AB)
  theory.ensure_levels(AB)
  return generate_report(AB, full_report)

# Full mode shows every concept in the sum.
# Diff mode only shows changed/new/removed concepts.

readable = {'++': 'kept, renamed',
            '==': 'kept',
            'o+': 'new',
            '+o': 'deprecated',
            }

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

    x_id = get_left_id(z, None)
    y_id = get_right_id(z, None)
    if y_id:
      y = z
      if x_id:
        rx = get_equated(y, None)  # how can this be missing??
        if rx: x = rx.record
    else:
      y = None
      x = z
    n1 = status_of_name(x)
    n2 = status_of_name(y)
    log("# status of %s (in A) is %s" % (blurb(x), n1))
    log("# status of %s (in B) is %s" % (blurb(y), n2))
    status_of_name(y)
    blob = readable.get(n1+n2) or "%s,%s" % (n1, n2)
    tick(blob)

  for key in sorted(list(stats.keys())):
    s = stats[key]
    log("%7d %7d %s" % (s[0], s[1], key))

# y = possessor of name

def status_of_name(y):
  if not y:
    return 'o'   # extension not present
  else:
    m = get_match(y, None)
    if m == None:
      return '+'   # no match (deprecated / new) ⏚
    else:
      r = theory.simple_relationship(y, m.record)
      if r == EQ:
        can1 = get_canonical(y, None)
        can2 = get_canonical(m.record, None)
        if can1 and can2 and can1 != can2:
          return '~'
      return rcc5_symbol(r)


def fooooo():
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
    name_changed = (blurb(x) != blurb(y))

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
        comment = "lumped into %s" % blurb(get_accepted(y))
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
      bx = blurb(x)
      if x and get_accepted(x, None): bx = bx + "*"
      by = blurb(y)
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


noises = {"subspecies": "",
          "species": "_",
          "subgenus": " _",
          "genus": "_ _",
          "subfamily": " _ _",
          "family": "_ _ _",
          "order": "_ _ _ _",
          "class": "_ _ _ _ _",
          }

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
