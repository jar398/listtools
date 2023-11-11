#!/usr/bin/env python3

import sys, csv, argparse
import rcc5, rows, checklist, workspace
import theory, exemplar, estimate

from workspace import ingest_workspace
from checklist import *
from rcc5 import rcc5_symbol

def generate_plugin_report(AB):
  yield ("A taxon id",
         "A taxon name",
         "B species that intersect",
         "LUB in B",
         "exemplar ids",
         )
  i = 0
  frequency = 1000
  for x in preorder_records(AB.A):
    if not is_top(x):
      u = AB.in_left(x)

      i += 1
      if i % frequency == 0:
        log("%s %s" % (i, blurb(u)))

      xids = theory.exemplar_ids(AB, u)
      if theory.is_species(u):
        o = theory.get_intersecting_species(u)
        inter = ';'.join(map(lambda s:show_articulation(u, s),
                             o))
        exemplars = ";".join(map(str, sorted(xids)))
      else:
        o = []
        inter = '-'
        exemplars = ";".join(map(str, sorted(xids))) if len(xids) <= 1 else '-'

      if is_accepted(x) or o:
        # filter out uninteresting synonyms
        est = estimate.get_estimate(u, None).record
        # est = local_accepted(est)
        yield (get_primary_key(x),
               blurb(x),
               inter,
               show_articulation(u, est),
               exemplars,
               )

def show_articulation(u, v):
  if v:
    rel = theory.compare(AB, u, v)
    y = get_outject(v)
    rcc5 = rcc5_symbol(rel.relationship)
    # leading space prevents interpretation as excel formula
    return " %s %s" % (rcc5_symbol(rel.relationship),
                      get_primary_key_for_report(y))
  else:
    return MISSING

def get_primary_key_for_report(x):
  assert get_primary_key(x)
  return "%s %s" % (get_primary_key(x), blurb(x))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    """)
  parser.add_argument('--A', help="the A checklist, as path name or -",
                      default='-')
  parser.add_argument('--B', help="the B checklist, as path name or -",
                      default='-')
  parser.add_argument('--exemplars', help="the exemplars table, as path name or -",
                      default=None)
  parser.add_argument('--Aname', help="short name of the A checklist",
                      default='A')
  parser.add_argument('--Bname', help="short name of the B checklist",
                      default='B')
  args=parser.parse_args()
  a_name = args.Aname
  b_name = args.Bname
  d_path = '-'
  with rows.open(args.A) as a_rows:
    with rows.open(args.B) as b_rows:
      AB = ingest_workspace(a_rows.rows(), b_rows.rows(),
                            A_name=a_name, B_name=b_name)
      if args.exemplars:
        exemplar.read_exemplars(rows.open(args.exemplars), AB)
      else:
        exemplar.find_exemplars(AB)
      log("# %s AB.exemplar_ufs" % len(AB.exemplar_ufs))
      theory.theorize(AB, False)
      with rows.open(d_path, "w") as d_gen:
        d_gen.write_rows(generate_plugin_report(AB))
