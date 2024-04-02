#!/usr/bin/env python3

import sys, csv, argparse
import rcc5, rows, checklist, workspace
import theory, exemplar, estimate

from workspace import ingest_workspace, is_accepted_locally, local_sup, \
  local_accepted, isinA, isinB
from checklist import *
from rcc5 import *
from specimen import sid_to_epithet
from estimate import get_block, is_empty_block, get_estimate
from property import mep, mep_get, mep_set
from typify import explain

counts = {}

def generate_row(AB, u, v, hom, est):
  if u == True:
    u_fields = ("A concept taxonID", "A concept name")
  else:
    u_fields = description_for_report(u)

  if v == True:
    v_fields = ("B concept taxonID", "B concept name")
  else:
    v_fields = description_for_report(v)

  rcc5_field = MISSING
  v_rel = None
  if v == True:
    rcc5_field = "RCC5"
  elif u and v:
    v_rel = theory.compare(AB, u, v)
    rcc5_field = rcc5_symbol(v_rel.relationship)
    if rcc5_field.startswith('='):
      rcc5_field = "'" + rcc5_field

  sidep = False
  if u == True:
    op_field = "operation"
  else:
    ops = impute_operation(AB, u, v_rel, hom)
    op_field = "; ".join(ops)
    sidep = "side" in ops
    for op in ["side"] if sidep else ops:
      if op in counts:
        counts[op] += 1
      else:
        counts[op] = 1

  if est == True:
    est_fields = ("containing B concept id", "containing B concept name")
  elif est.relationship != EQ and not sidep:
    est_fields = description_for_report(est.record)
  else:
    est_fields = description_for_report(None)

  if u == True:
    u_and_v = "in both A and B concepts"
    u_not_v = "in A concept but not in B concept"
    v_not_u = "in B concept but not in A concept"
  elif u and v:
    u_ids = get_block(u)
    v_ids = get_block(v)
    u_and_v = show_sid_set(AB, u_ids & v_ids)
    u_not_v = show_sid_set(AB, u_ids - v_ids)
    v_not_u = show_sid_set(AB, v_ids - u_ids)
  else:
    u_and_v = MISSING
    u_not_v = MISSING
    v_not_u = MISSING

  return (op_field,
          u_fields[0], u_fields[1],
          rcc5_field,
          v_fields[0], v_fields[1],
          est_fields[0], est_fields[1],
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

        est = get_estimate(u)

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

# s is a set of exemplar ids...

def show_sid_set(AB, block):
  # xep = (sid, epithet)
  # Sort by epithet
  xeps = sorted(map(lambda sid:(sid, sid_to_epithet(AB, sid)), block),
                key=lambda xep: xep[1])
  return "; ".join(map(lambda xep:"%s %s" % xep, xeps))

# Could be:
#   Change of rank (promotion/demotion).
#   Change of acceptedness (accepted/synonymized).
#   Change of epithet / gender (regendered).
#   Change of genus (moved).
#   Change of concept (lumped / split / changed).
#   Possible change of concept (perhaps lumped / perhaps split)

def impute_operation(AB, u, v_rel, hom):
  ops = []

  if hom == "hom?":
    hom2 = "type?"
  elif hom:                     # "hom"
    hom2 = None
  else:
    hom2 = "side"

  v = v_rel.record if v_rel else None

  if monitor(u):
    x = get_outject(u)
    log("# %s %s has dup_id %s" % (get_primary_key(x), blurb(u),
        get_duplicate_from(x, "(none)")))
  if monitor(v):
    y = get_outject(v)
    log("# %s %s has dup_id %s" % (get_primary_key(y), blurb(v),
        get_duplicate_from(y, "(none)")))

  if u and get_duplicate_from(get_outject(u), None):
    ops.append("dup in A")
  if v and get_duplicate_from(get_outject(v), None):
    ops.append("dup in B")

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


def description_for_report(u):
  if u:
    x = get_outject(u)
    return (get_primary_key(x),
            "%s%s" % (get_canonical(u, None) or get_scientific(u, None),
                      "" if get_accepted(x) else "*"))
  else:
    return (MISSING, MISSING)

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
      theory.theorize(AB, False)
      with rows.open(d_path, "w") as d_gen:
        d_gen.write_rows(generate_plugin_report(AB))
