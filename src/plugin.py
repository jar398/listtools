#!/usr/bin/env python3

import sys, csv, argparse
import rcc5, rows, checklist, workspace
import theory, exemplar, estimate

from workspace import ingest_workspace, is_accepted_locally, local_sup, \
  isinA, isinB
from checklist import *
from rcc5 import *
from estimate import get_estimate, get_block, is_empty_block
from property import mep, mep_get, mep_set
from typify import xid_epithet, explain, compare_records
def generate_plugin_report(AB):
  yield ("A taxon id",
         "A name",
         "operation",
         "B concept via type",
         "intersecting B concepts (species only)",
         "containing B concept",
         # "change from A to B",
         "in both A and B concepts",
         "in A concept but not in B concept",
         "in B concept but not in A concept",
         )
  i = [0]
  frequency = 5000

  # species in B that should be adopted by taxa in A
  adoptees = find_adoptees(AB)

  def process_subtree(AB, x):
    yield from process_record(x)
    infs = list(get_inferiors(x))
    infs.sort(key=plugin_sort_key)
    for c in infs:
      yield from process_subtree(AB, c)
    yield from process_new_species(AB.in_left(x))

  # TBD: Deal with situation where one of these comes from a
  # subspecies or synonym in A
  def process_new_species(u):
    for v in mep_get(adoptees, u, ()):
      # v a non-synonym species
      y = get_outject(v)
      yield (MISSING,           # concept not in A
             MISSING,           # concept not in A
             "added name",
             ". " + show_relation(relation(NOINFO, v)),
             MISSING,           # intersecting
             ". " + show_relation(theory.get_central(AB, v)),     # lub in B?
             MISSING,                          # in A and B
             MISSING,                          # in A but not B
             MISSING,                          # in B but not A
             )

  def plugin_sort_key(x):
    if is_accepted(x):
      return (1, blurb(x))
    else:
      return (0, blurb(x))

  def process_record(x):
    if not is_top(x):
      u = AB.in_left(x)
      if theory.is_species(u):

        i[0] += 1
        if i[0] % frequency == 0:
          log("%s %s" % (i[0], blurb(u)))

        operation = None

        lub = show_relation(estimate.get_estimate(u, None))
        if lub != MISSING: lub = ". " + lub

        inters = theory.get_intersecting_species(u)
        rels = list(map(lambda w: theory.compare(AB, u, w), inters))

        # or, if get_species(u) and bud == v
        bud = theory.get_buddy(AB, u)   # B taxon containing u's type
        v = theory.get_species(bud) if bud else None
        if v:
          if not v in inters:
            log("** Intersection expected but missing: %s ?! %s" %
                (blurb(u), blurb(v)))
          rels2 = (rel for rel in rels if not v is rel.record)
        else:
          rels2 = rels
          if len(rels) == 1:  # and rels[0].relationship == EQ:
            # Curation issue probably (e.g. MSW3)
            v1 = rels[0].record
            #log("** Concepts inferred to be homotypic: %s, %s" %
            #    (blurb(u), blurb(v1)))
            v = v1

        inter = '. ' + ';'.join(map(show_relation, rels2))
        if v:
          assert isinB(AB, v), (blurb(bud), blurb(v))
          operation = impute_operation(AB, u, bud, v)
          u_ids = get_block(u)
          v_ids = get_block(v)
          u_and_v = show_xid_set(AB, u_ids & v_ids)
          u_not_v = show_xid_set(AB, u_ids - v_ids)
          v_not_u = show_xid_set(AB, v_ids - u_ids)
          if bud:
            via_type = ". " + show_relation(theory.compare(AB, u, v))
          else:
            via_type = MISSING

          yield (get_primary_key(x),
                 blurb(x),
                 operation or '-',
                 via_type,
                 inter,
                 lub,
                 u_and_v, u_not_v, v_not_u,
                 )

        else:                   # No bud, no v
          # How can we have intersecting species but no buddy??
          # !!! They match only via synonyms, not via the type.
          # assert len(inters) == 0, blurb(u)

          # or (v and theory.is_species(v) and v is bud)):
          yield (get_primary_key(x),
                 blurb(x),
                 "unclear" if rels else "removed name",
                 MISSING,
                 inter,
                 lub,
                 MISSING, MISSING, MISSING,
                 )

  counts = {}
  for row in process_subtree(AB, AB.A.top):
    op = row[2]
    if op in counts:
      counts[op] += 1
    else:
      counts[op] = 1
    yield row
  for (op, count) in counts.items():
    log("* %6d %s" % (count, op))

# s is a set of exemplar ids

def show_xid_set(AB, s):
  return "{%s}" % ",".join(map(lambda s:"%s %s" % (s, xid_epithet(AB, s)),
                               sorted(s)))

# Could be:
#   Change of rank (promotion/demotion).
#   Change of acceptedness (accepted/synonymized).
#   Change of epithet / gender (regendered).
#   Change of genus (moved).
#   Change of concept (lumped / split / changed).
#   Possible change of concept (perhaps lumped / perhaps split)

def impute_operation(AB, u, bud, v):
  ops = []
  if bud is v:
    if get_rank(u, None) != get_rank(v, None):
      if get_rank(u, None) == "subspecies":
        ops.append("promoted")
      elif get_rank(v, None) == "subspecies":
        ops.append("demoted")
      else:
        ops.append("change of rank")
    if (is_accepted_locally(AB, u) !=
        is_accepted_locally(AB, v)):
      if (is_accepted_locally(AB, u)):
        ops.append("synonymized")
      else:
        ops.append("accepted")

  rel = theory.compare(AB, u, v)
  if rel.relationship == LT:
    ops.append("lumped")
  #elif rel.relationship == LE:
  #  ops.append("perhaps lumped")
  elif rel.relationship == GT:
    ops.append("split")
  #elif rel.relationship == GE:
  #  ops.append("perhaps split")
  elif rel.relationship == EQ:
    p1 = get_parts(get_outject(u))
    p2 = get_parts(get_outject(v))
    if p1.genus != p2.genus:
      ops.append("moved")
    if p1.epithet != p2.epithet:
      ops.append("epithet")
  else:                         # OVERLAP etc.
    ops.append("concept")
  if len(ops) == 0:
    return "unchanged"
  else:
    return ";".join(ops)

# Find species in B that could have been split off from taxa in A

def find_adoptees(AB):
  ad = mep()
  for y in all_records(AB.B):
    v = AB.in_right(y)
    if theory.is_species(v):  # implies is_accepted_locally(v)
      if is_empty_block(get_block(v)):
        est = get_estimate(v)
        u = est.record
        ads = mep_get(ad, u, None)
        if ads == None:
          ads = []
          mep_set(ad, u, ads)
        ads.append(v)
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
