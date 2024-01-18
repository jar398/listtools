#!/usr/bin/env python3

import sys, csv, argparse
import rcc5, rows, checklist, workspace
import theory, exemplar, estimate

from workspace import ingest_workspace, is_accepted_locally, local_sup
from checklist import *
from rcc5 import *
from estimate import get_estimate, get_block, is_empty_block
from property import mep, mep_get, mep_set

def generate_plugin_report(AB):
  yield ("A taxon id",
         "A taxon name",
         "operation",
         "intersecting B concepts (species only)",
         "containing B concept",
         # "change from A to B",
         "kept",
         "removed",
         "added",
         )
  i = [0]
  frequency = 5000

  adoptees = find_adoptees(AB)

  def process_subtree(x):
    yield from process_record(x)
    infs = list(get_inferiors(x))
    infs.sort(key=plugin_sort_key)
    for c in infs:
      yield from process_subtree(c)

    u = AB.in_left(x)
    for v in mep_get(adoptees, u, ()):
      v_sup = local_sup(AB, v)
      v_ids = get_block(v)
      b_not_a = show_xid_set(v_ids)
      yield (MISSING,           # not in A
             MISSING,           # not in A
             "de novo or split off",
             ". " + show_relation(relation(NOINFO, v)),  # inters in B
             ". " + show_relation(v_sup),      # lub in B
             MISSING,                          # in A and B
             MISSING,                          # in A but not B
             b_not_a,                          # in B but not A
             )

  def plugin_sort_key(x):
    if is_accepted(x):
      return (1, blurb(x))
    else:
      return (0, blurb(x))

  def process_record(x):
    if not is_top(x):
      u = AB.in_left(x)

      i[0] += 1
      if i[0] % frequency == 0:
        log("%s %s" % (i[0], blurb(u)))

      operation = None

      a_and_b = MISSING
      a_not_b = MISSING
      b_not_a = MISSING

      p = local_sup(AB, u).record
      if theory.is_species(u) or theory.is_species(p):
        rels = list(map(lambda v: theory.compare(AB, u, v),
                        theory.get_intersecting_species(u)))
        if len(rels) == 0:
          operation = "removed or lumped"
        elif len(rels) > 1:
          operation = "split"
        else:
          operation = impute_operation(u, rels[0])

        # Among rels, find the one that is to a B concept of the same name...?
        # That would v such that get_exemplar(u) is in v's block.
        u_ids = get_block(u)
        v = theory.get_buddy(AB, u)
        if v:
          v_ids = get_block(v)
          a_and_b = show_xid_set(u_ids & v_ids)
          a_not_b = show_xid_set(u_ids - v_ids)
          b_not_a = show_xid_set(v_ids - u_ids)
        if a_not_b == MISSING:
          a_not_b = show_xid_set(u_ids)
        # Sort ???  there are usually only two
        inter = '. ' + ';'.join(map(show_relation, rels))
      else:
        # Not a species
        rels = []
        inter = '-'

      # filter out uninteresting synonyms
      if is_accepted(x) or len(rels) > 0:
        est_rel = estimate.get_estimate(u, None)
        if not operation:
          operation = impute_operation(u, est_rel)
        lub = show_relation(est_rel)
        if lub != MISSING: lub = ". " + lub
        yield (get_primary_key(x),
               blurb(x),
               operation or '-',
               inter,
               lub,
               a_and_b, a_not_b, b_not_a,
               )

  yield from process_subtree(AB.A.top)

def show_xid_set(s):
  # was return ";".join(map(str, sorted(s))) + ";"
  return "{%s}" % ",".join(map(str, sorted(s)))

def impute_operation(u, rel):
  v = rel.record
  if rel.relationship == LT:
    operation = "lumped"
  elif rel.relationship == EQ:
    if (is_accepted_locally(AB, u) and
        not is_accepted_locally(AB, v)):
      operation = "synonymized"
    elif (not is_accepted_locally(AB, u) and
          is_accepted_locally(AB, v)):
      # won't happen since synonyms aren't reported on, but
      operation = "accepted"
    elif get_canonical(u, None) == get_canonical(v, None):
      operation = "unchanged"
    else:
      operation = "renamed"     # ? could be added or removed
  elif ((rel.relationship & EQ) != 0 and
        get_canonical(u, None) == get_canonical(v, None)):
    operation = "related"
  # Else: LE GE NOINFO  (never DISJOINT, OVERLAP, INCONSISTENT)
  # I think: never GT DISJOINT
  else:
    operation = "other"                # Not sure what to say
  return operation

# Find species in B that should be adoped by taxa in A

def find_adoptees(AB):
  ad = mep()
  for y in all_records(AB.B):
    v = AB.in_right(y)
    if theory.is_species(v) and is_empty_block(get_block(v)):
      est = get_estimate(v)
      u = est.record
      ads = mep_get(ad, u, None)
      if ads == None: ads = []
      ads.append(v)
      mep_set(ad, u, ads)
  return ad

def is_infraspecific(u):
  x = get_outject(u)
  r = get_rank(x, None)
  # Kludge.  Will miss some.  Can improve this.
  return r == 'subspecies' or r == 'infraspecific name'    # COL


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
