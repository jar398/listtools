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

def report(merged_iter, full_report):
  AB = merge.rows_to_merged(merged_iter, {'name': 'AB'})
  ensure_inferiors_indexed(AB)
  for r in all_records(AB): assert get_superior(r)
  theory.ensure_levels(AB)
  merge.resolve_merge_links(AB)
  return generate_report(AB, full_report)

readable_template = "belongs to %s extension"

# Full mode shows every concept in the sum.
# Diff mode only shows changed/new/removed concepts.

readable = {'@@': 'kept, renamed',
            '==': 'kept, not renamed',
            'o@': 'only in B',
            '@o': 'only in A',
            'o': 'not present',
            '@': 'no match',
            '=': readable_template % 'the same',
            '<': readable_template % 'a wider',
            '>': readable_template % 'a narrower',
            '><': readable_template % 'a conflicting',
            '!': readable_template % 'an unrelated',
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
    diff = get_difference(z, None)
    x = y = z
    if diff == "not in A":
      x = None
    elif diff == "matched":
      continue
    else:
      if diff == "not in B":
        y = None
      else:
        # There's an equation, so on this row (the B row) there should be
        #   B -> A equated link
        # and on the corresponding A record there should be
        #   A -> B synonym link
        rx = get_equated(z, None)  # how can this be missing??
        assert rx
        x = rx.record
        y = z
    n1 = status_of_name(x)
    n2 = status_of_name(y)
    if False:
      log("# '%s' in A names a %s extension in B" % (blurb(x), n1))
      log("# '%s' in B names a %s extension in A" % (blurb(y), n2))
    r1 = "A name %s" % readable.get(n1, '?') if x else "not in A"
    r2 = "B name %s" % readable.get(n2, '?') if y else "not in B"
    blob = (readable.get(n1+n2) or "%s, %s" % (r1, r2))
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
      return '@'   # no match (deprecated / new) ⏚
    else:
      r = theory.simple_relationship(y, m.record)
      if r == EQ:
        can1 = get_canonical(y, None)
        can2 = get_canonical(m.record, None)
        if can1 and can2 and can1 != can2:
          return '~'
      return rcc5_symbol(r)

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
