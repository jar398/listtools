#!/usr/bin/env python3

import sys, csv, argparse
import rcc5, rows, checklist, workspace
import theory, exemplar, estimate

from workspace import ingest_workspace, is_accepted_locally
from checklist import *
from rcc5 import *
from estimate import get_estimate

def generate_plugin_report(AB):
  yield ("A taxon id",
         "A taxon name",
         "operation",
         "intersecting B concepts (species only)",
         "containing B concept",
         # "change from A to B",
         "exemplar ids",
         )
  i = 0
  frequency = 5000
  for x in preorder_records(AB.A):
    if not is_top(x):
      u = AB.in_left(x)

      i += 1
      if i % frequency == 0:
        log("%s %s" % (i, blurb(u)))

      operation = MISSING

      if theory.is_species(u):  # or is_infraspecific(u):  ??
        rels = list(map(lambda v: theory.compare(AB, u, v),
                        theory.get_intersecting_species(u)))
        if len(rels) == 0:
          operation = "removed"
        elif len(rels) > 1:
          operation = "split"
        else:
          rel = rels[0]
          v = rel.record
          if rel.relationship == LT:
            operation = "lumped"
          elif rel.relationship == EQ:
            if (get_canonical(u) and get_canonical(v) and
                get_canonical(u) != get_canonical(v)):
              operation = "renamed"
            if (is_accepted(x) and
                not is_accepted_locally(AB, v)):
              operation = "synonymized"
            if (not is_accepted(x) and
                is_accepted_locally(AB, v)):
              operation = "desynonymized"
            else:
              operation = "unchanged"
          elif ((rel.relationship & EQ) != 0 and
                get_canonical(u) and get_canonical(v) and
                get_canonical(u) == get_canonical(v)):
            operation = "related"
          # Else: LE GE NOINFO  (never DISJOINT, OVERLAP, INCONSISTENT)
          # I think: never GT DISJOINT
          else:
            operation = "other"                # Not sure what to say

        inter = '. ' + ';'.join(map(show_relation, rels))
        xids = estimate.exemplar_ids(AB, u)
        exemplars = ";".join(map(str, sorted(xids)))
      else:
        rels = []
        inter = '-'
        exemplars = '-'
        operation = "de novo"

      # filter out uninteresting synonyms
      if is_accepted(x) or len(rels) > 0:
        est = estimate.get_estimate(u, None).record
        lub = show_articulation(AB, u, est)
        if lub != MISSING: lub = ". " + lub
        # lump    if inter relationship is <
        # rename  if relationship is = and canonical differs
        # ?  if synonym promoted to accepted (sort of a split)
        else:
          status = MISSING
        yield (get_primary_key(x),
               blurb(x),
               operation,
               inter,
               lub,
               exemplars,
               )

def is_infraspecific(u):
  x = get_outject(u)
  r = get_rank(x, None)
  # Kludge.  Will miss some.  Can improve this.
  return r == 'subspecies' or r == 'infraspecific name'    # COL


def show_articulation(AB, u, v):
  return show_relation(theory.compare(AB, u, v))

def show_relation(rel):
  if not is_toplike(rel.record):
    y = get_outject(rel.record)
    rcc5 = rcc5_symbol(rel.relationship)
    # leading space prevents interpretation as excel formula
    return "%s %s" % (rcc5_symbol(rel.relationship),
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
        exemplar.find_exemplars(get_estimate, AB)
      log("# %s AB.exemplar_ufs" % len(AB.exemplar_ufs))
      theory.theorize(AB, False)
      with rows.open(d_path, "w") as d_gen:
        d_gen.write_rows(generate_plugin_report(AB))
