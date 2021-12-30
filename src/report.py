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

get_difference = prop.getter(merge.difference_prop)

#

def report(merged_iter, mode):
  AB = merge.rows_to_merged(merged_iter, {'name': 'AB'})
  merge.resolve_merge_links(AB)
  for r in all_records(AB): assert get_superior(r)  # sanity
  ensure_inferiors_indexed(AB)
  #theory.ensure_levels(AB)
  return generate_report(AB, mode)

readable_template = "%s"

# Full mode shows every concept in the sum.
# MDD mode only shows changed/new/removed rows.

def generate_report(AB, mode):
  yield ("A name", "B name", "rank", "differences", "match_note")
  stats = {}

  def tick(stat, y, rank):
    if stat in stats:
      s = stats[stat]
    else:
      s = [0, 0]      # [all, species only]
      stats[stat] = s
    s[0] += 1
    if rank == "species":
      s[1] += 1

  for z in all_records(AB):
    (x, y) = personae(z)
    if x or y:
      diff = get_difference(z, None)
      rank = ((y and get_rank(y, None)) or (x and get_rank(x, None)))
      tick(diff or 'same', z, rank)
      if mode == 'full':
        doit = True
      elif mode == 'mdd':
        doit = (diff and rank == 'species')
      else:
        assert False   # invalid mode
      if doit:
        yield (blurb(x) if x else '',
               blurb(y) if y else '',
               rank, diff,
               get_match_note(y, '') if y else '')

  for key in sorted(list(stats.keys())):
    s = stats[key]
    log("%7d %7d %s" % (s[0], s[1], key))

def personae(z):
  if is_top(z):
    return (None, None)
  if get_equated(z, None):
    return (get_equated(z).record, z)  # B with equated A
  else:
    rs = get_superior(z)
    if rs.relationship == EQ:        # A equated to B - ignore
      return (None, None)
    else:
      diff = get_difference(z, '')
      # this is a kludge. there ought to be a more direct test
      if 'not in A' in diff:
        return (None, z)
      elif 'not in B' in diff:
        return (z, None)
      else:
        assert False

# -----------------------------------------------------------------------------

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    standard input = the merged checklist we're reporting on
    """)
  parser.add_argument('--mode',
                      default='mdd',
                      help="report type: 'mdd' or 'full'")
  args=parser.parse_args()

  rep = report(csv.reader(sys.stdin), args.mode)
  util.write_rows(rep, sys.stdout)
