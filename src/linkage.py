#!/usr/bin/env python3

# Record linkage

import util
import parse, rows
import simple
from checklist import *
from workspace import *
from parse import parse_name
from typify import compare_parts, homotypic_comparison, \
  compute_distance, compare_records, pick_better_record, \
  really_get_typification_uf, get_link

# For each record, get list of matching records.... hmm
# Future: read exceptions or whole thing from a file

# Definition of 'subproblem':
# If a node pair (u, v) doesn't have u in us and v in vs for some
# 'subproblem' (us, vs) then it is not going to be considered a
# potential match.

# 'pre_estimate' = LUB estimate from first pass, if this is second pass.


# Informational.  Requires a rewrite for exemplars.

def report_on_links(AB, subproblems):
  full_count = 0
  half_count = 0
  amb_count = 0
  for (key, (us, vs)) in subproblems.items():
    for u in us:
      link = get_link(u, None)
      if link == None:
        pass
      elif link == False:
        amb_count += 1
      else:
        link2 = get_link(link, None)
        if link2 is u:
          full_count += 1
        elif link2 == None:
          half_count += 1
    for v in vs:
      link = get_link(u, None)
      if link == None:
        pass
      elif link == False:
        amb_count += 1
      else:
        link2 = get_link(link, None)
        if link2 is v:
          pass                  # already counted
        elif link2 == None:
          half_count += 1
  log("#   Links: %s mutual, %s one-sided, %s ambiguous" %
      (full_count, half_count, amb_count))

# -----------------------------------------------------------------------------
# Plumbing

def generate_linkage_report(AB):
  yield ('from', 'to', 'comparison')
  for x in all_records(AB.A):
    u = AB.in_left(x)
    v = get_link(u, None)
    if v:
      yield (blurb(u), blurb(v), compare_records(u, v))    # distance
    elif v == False:
      yield (blurb(u), 'ambiguous', MISSING)

import argparse
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Generate list of exemplars proposed for two checklists
    """)
  parser.add_argument('--A', help="the A checklist, as path name or -",
                      default='-')
  parser.add_argument('--B', help="the B checklist, as path name or -",
                      default='-')
  parser.add_argument('--test', action='store_true', help="run smoke tests")
  args=parser.parse_args()

  if args.test:
    print(compare_parts(parse_name("Sturnira angeli"),
                        parse_name("Sturnira magna"),
                        3))
  else:
    a_name = 'A'; b_name = 'B'
    a_path = args.A
    b_path = args.B
    with rows.open(a_path) as a_rows: # rows object
      with rows.open(b_path) as b_rows:
        # compute name matches afresh
        AB = ingest_workspace(a_rows.rows(), b_rows.rows(),
                              A_name=a_name, B_name=b_name)
        find_links(AB)
        report_gen = generate_linkage_report(AB)
        util.write_rows(report_gen, sys.stdout)
