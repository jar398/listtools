#!/usr/bin/env python3

import sys, csv, argparse
import rcc5, rows, checklist, workspace
import theory, exemplar, estimate

from workspace import ingest_workspace, is_accepted_locally, local_sup, \
  local_accepted, isinA, isinB
from checklist import *
from rcc5 import *
from estimate import get_estimate, get_block, is_empty_block
from property import mep, mep_get, mep_set
from typify import xid_epithet, explain, compare_records

counts = {}

def generate_row(AB, u, v, hom, est):
  if u == True: u_field = "A concept"
  elif u: u_field = description_for_report(u)
  else: u_field = MISSING

  if v == True:
    v_field = "B concept"
  elif v:
    v_field = description_for_report(v)
  else:
    v_field = MISSING

  rcc5_field = MISSING
  if v == True:
    v_rel = None
    v_rel_field = "B concept"
    rcc5_field = "RCC5"
  elif u and v:
    v_rel = theory.compare(AB, u, v)
    v_rel_field = ". " + show_relation(v_rel)
    rcc5_field = rcc5_symbol(v_rel.relationship)
  elif v:
    v_rel = relation(NOINFO, v)
    v_rel_field = ". " + show_relation(v_rel)
  else:
    v_rel = None
    v_rel_field = MISSING

  if u == True:
    op_field = "operation"
  else:
    ops = impute_operation(AB, u, v_rel, hom)
    op_field = "; ".join(ops)
    for op in ["side"] if "side" in ops else ops:
      if op in counts:
        counts[op] += 1
      else:
        counts[op] = 1

  est_field = MISSING
  if est == True:
    est_field = "containing B concept"
  elif est:
    if est.relationship != EQ:
      est_field = description_for_report(est.record)

  if u == True:
    u_and_v = "in both A and B concepts"
    u_not_v = "in A concept but not in B concept"
    v_not_u = "in B concept but not in A concept"
  elif u and v:
    u_ids = get_block(u)
    v_ids = get_block(v)
    u_and_v = show_xid_set(AB, u_ids & v_ids)
    u_not_v = show_xid_set(AB, u_ids - v_ids)
    v_not_u = show_xid_set(AB, v_ids - u_ids)
  else:
    u_and_v = MISSING
    u_not_v = MISSING
    v_not_u = MISSING

  return (op_field,
          u_field,
          rcc5_field,
          v_field,
          est_field,
          u_and_v,
          u_not_v,
          v_not_u,
          )

def generate_plugin_report(AB):
  yield generate_row(AB, True, True, True, True) # header
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
    vs = mep_get(adoptees, u, ())
    if len(vs) > 0:
      vs.sort(key=plugin_sort_key)
      for v in vs:
        # v an accepted species
        assert v
        yield generate_row(AB, None, v, False,
                           theory.get_central(AB, v))

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
        if i[0] % frequency == 0: log("# %s %s" % (i[0], blurb(u)))

        est = estimate.get_estimate(u, None)
        lub = ". " + show_relation(est)

        inters = theory.get_intersecting_species(u)
        rels = list(map(lambda w: theory.compare(AB, u, w), inters))

        bud = theory.get_buddy(AB, u)
        if bud:
          v = theory.get_species(bud)
        else:
          v = hom = None  

        if len(inters) > 0:
          if v:
            assert isinB(AB, v), blurb(v)
            hom = "hom"
          else:
            v = min(inters, key=lambda v2:hamming(u, v2))
            hom = "hom?"        # op = "type?"
          yield generate_row(AB, u, v, hom, est)
          for v2 in sorted(inters, key=plugin_sort_key):
            if not v2 is v:
              yield generate_row(AB, u, v2, False, est) # op = "side"
        else:                   # No bud, no v
          yield generate_row(AB, u, None, False, est)

  yield from process_subtree(AB, AB.A.top)

  for (op, count) in counts.items():
    log(" %6d %s" % (count, op))

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

def impute_operation(AB, u, v_rel, hom):
  v = v_rel.record if v_rel else None
  ops = []

  if hom == "hom?":
    hom2 = "type?"
  elif hom:                     # "hom"
    hom2 = None
  else:
    hom2 = "side"

  if u == None:
    ops.append("added name")
  elif v_rel == None:
    ops.append("removed name")
  else:
    if hom2:
      ops.append(hom2)          # ?? what to put here ??

    if v_rel.relationship == LT:
      ops.append("lumped")
    elif v_rel.relationship == GT:
      ops.append("split")
    elif v_rel.relationship == OVERLAP:
      ops.append("overlaps")
    elif v_rel.relationship == DISJOINT:
      ops.append("unrelated")      # shouldn't happen
    elif EQ & v_rel.relationship != 0: # EQ LE GE NOINFO INTERSECT
      ops.append("congruent")
    else:
      assert False

    # Not used in species-only reports
    if hom:
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

    if hom:
      p1 = get_parts(get_outject(u))
      p2 = get_parts(get_outject(v))
      if p1.genus != p2.genus:
        ops.append("moved")
      if p1.epithet != p2.epithet:
        ops.append("epithet")
      if p1.token != p2.token:
        ops.append("author")
      if p1.year != p2.year:
        ops.append("year")

  return ops

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
    rcc5 = rcc5_symbol(rel.relationship)
    # leading space prevents interpretation as excel formula
    return "%s %s" % (rcc5_symbol(rel.relationship),
                      description_for_report(rel.record))
  else:
    return MISSING

def description_for_report(u):
  x = get_outject(u)
  return "%s.%s.%s%s" % (get_source_tag(x),
                         get_primary_key(x),
                         get_canonical(u),
                         "" if get_accepted(x) else "*")

def hamming(u, v):
  p1 = get_parts(get_outject(u))
  p2 = get_parts(get_outject(v))
  h = 0
  if p1.genus != p2.genus: h += 1
  if p1.epithet != p2.epithet: h += 1
  if p1.token != p2.token: h += 1
  if p1.year != p2.year: h += 1
  return h

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
