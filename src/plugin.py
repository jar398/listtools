#!/usr/bin/env python3

import sys, csv, argparse
import rcc5, rows, checklist, workspace
import theory, exemplar, estimate
import jumble
import util

from workspace import ingest_workspace, is_accepted_locally, local_sup, \
  local_accepted, isinA, isinB, is_species
from checklist import *
from rcc5 import *
from specimen import same_type_ufs, sid_to_epithet
from estimate import get_equivalent
from block import get_block, is_empty_block
from property import mep, mep_get, mep_set
from ranks import ranks_dict
from util import reset_log_allowance

counts = {}
def count(tag):
  if tag in counts: counts[tag] += 1
  else: counts[tag] = 1

# Preorder traversal of AB.A

def generate_plugin_report(AB):
  yield generate_row(AB, True, True, True) # header
  i = [0]
  frequency = 5000

  jumble.jumble_workspace(AB)
  spa = spb_spe = spb_nspe = spb_ne = 0

  seen = set()

  for z in preorder_records(AB):
    count("records")
    i[0] += 1
    if i[0] % frequency == 0:
      log("# %s %s" % (i[0], blurb(z)))
    ops = []

    # A -
    # - B
    # A A
    # A B
    
    assert not jumble.is_redundant(AB, z)

    if is_species(z):
      w = choose_partner(AB, z) # may be None
      if isinA(AB, z):
        u = z; v = w
        spa += 1
      else:
        u = w; v = z
        spb_ne += 1

      combo = (get_primary_key(u) if u else None,
               get_primary_key(v) if v else None)
      if not combo in seen:
        seen.add(combo)
        yield generate_row(AB, u, v, ops)

def op_counts_report(counts):
  reset_log_allowance()    # prints to stderr
  for (op, op_count) in counts.items():
    log("  %6d %s" % (op_count, op))

# For each species, we pick one articulation.  The articulation is
# between the species and a "partner".
# z is a species, but partner could be non-species or None.

def choose_partner(AB, z):
  e_rel = get_equivalent(AB, z)
  if e_rel:
    sp = theory.get_species(e_rel.record)
    return sp or e_rel.record

  # ... no congruent record, pick the 'best' one from intersectors ...
  inters = theory.get_intersecting_species(z)
  if len(inters) == 0:
    return None                 # No intersecting species
  else:
    z_ids = get_block(z)
    def intensity(v):
      v_ids = get_block(v)
      # Ad hoc rule:
      # Prefer: 1. same type specimen, 2. maximum overlap, 3. minimum nonoverlap
      same = 0 if same_type_ufs(z, v) else 1
      meet = len(z_ids & v_ids)
      nono = len(z_ids | v_ids) - meet
      return (same, -meet, nono)
    return min(inters, key=intensity)


# Return either the header row or a data row; the code for both is
# together to make sure they stay in sync.
# u or v might be None, or non-species.

def generate_row(AB, u, v, ops):
  # A - rcc5 - B columns
  if u == True:
    u_fields = ("A concept taxonID", "A concept name")
    rcc5_field = "RCC5"         # header
    rcc5_tag = "tag"
    v_fields = ("B concept taxonID", "B concept name")
  else:
    u_fields = description_for_report(AB, u)
    v_fields = description_for_report(AB, v)
    if u and v:
      v_rel = theory.compare(AB, u, v)
      (rcc5_field, rcc5_tag) = impute_relationship(AB, u, v_rel.relationship)
      if rcc5_field.startswith('='):
        # Kludge to appease spreadsheet programs
        rcc5_field = "'" + rcc5_field
    else:
      rcc5_field = MISSING
      if not u:
        tag = "not in A"
      elif not v_rel:
        tag = "not in B"
      else:
        tag = "should not happen"
      v_rel = None

  # change column
  if u == True:
    op_field = "change"      # header.  obsolete
  else:
    # `homotypic` is true iff u is homotypic with v.
    homotypic = u and v and same_type_ufs(u, v)

    if False:
      ops += impute_concept_change(AB, u, v_rel, homotypic)
      ops += impute_name_change(AB, u, v_rel, homotypic)
      # ops += impute_rank_change(AB, u, v_rel, homotypic)
      op_field = "; ".join(ops)
      for op in ops:
        if op in counts:
          counts[op] += 1
        else:
          counts[op] = 1
    else:
      if tag in counts:
        counts[tag] += 1
      else:
        counts[tag] = 1

  # Venn diagram columns
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

  return (;; op_field,
          u_fields[0], u_fields[1],
          rcc5_field,
          v_fields[0], v_fields[1],
          u_and_v,
          u_not_v,
          v_not_u,
          )

# s is a set of exemplar ids...

def show_sid_set(AB, block):
  if True:
    sids = sorted(block, key=lambda sid: sid_to_epithet(AB, sid))
    return "; ".join(map(show_sid, sids))
  else:
    # xep = (sid, epithet)
    xeps = sorted(map(lambda sid:(sid, sid_to_epithet(AB, sid)), block),
                  key=lambda xep: xep[1])
    return "; ".join(map(lambda xep:"%s %s" % xep, xeps))

def show_sid(sid):
  ep = sid_to_epithet(AB, sid)
  return "%s %s" % (sid, ep)

def impute_relationship(AB, u, v_rel):
  ship = v_rel.relationship
  tag = ""
  if ship == EQ: tag = "is congruent with"
  elif ship == LT: tag = "is contained in"
  elif ship == GT: tag = "contains"
  elif ship == OVERLAP: tag = "overlaps"
  return (rcc5_symbol(ship), tag)

def impute_concept_change(AB, u, v_rel, homotypic):
  if not u:
    op = "new"
  elif not v_rel:
    op = "gone"
  else:
    ship = v_rel.relationship
    if ship == EQ:
      # This can happen if e.g. the type is ambiguous or otherwise unmatched
      op = "congruent"
      if not homotypic:
        from specimen import get_specimen_id, get_type_uf
        v = v_rel.record
        log("# u %s   %s" %
            (show_sid(get_specimen_id(get_type_uf(u))),
             blurb(u)))
        log("# v %s   %s" %
            (show_sid(get_specimen_id(get_type_uf(v))),
             blurb(v)))
    elif ship == LT:
      op = "expand" if homotypic else "lump"
    elif ship == GT:
      op = "contract" if homotypic else "split"
    elif ship == OVERLAP:
      op = "contract+expand" if homotypic else "split+lump"
    else:                       # DISJOINT  ?
      assert not homotypic
      op = "should not happen"
  return [op]

def impute_name_change(AB, u, v_rel, homotypic):
  ops = []
  if homotypic:
    v = v_rel.record
    p1 = get_parts(get_outject(u))
    p2 = get_parts(get_outject(v))
    if p1.genus != p2.genus:
      ops.append("change genus")
    if p1.epithet != p2.epithet:
      ops.append("change epithet")
    if p1.token != p2.token:
      ops.append("change author")
    if p1.year != p2.year:
      ops.append("change year")
    r1 = get_rank(get_outject(u), None)
    r2 = get_rank(get_outject(v), None)
    if r1 != r2:
      ops.append("change rank")
  # maybe: promotion, demotion, rank change
  return ops

# u or v might be None.
# Not ready for prime time.

def impute_rank_change(AB, u, v_rel, homotypic):
  ops = []
  if not u or not v_rel: return ops
  v = v_rel.record
  r1 = ranks_dict[get_rank(u, None)]
  r2 = ranks_dict[get_rank(v, None)]
  if r1 and r2:
    if r1 < r2:
      ops.append("promoted")
    elif r1 > r2:
      ops.append("demoted")
  a1 = is_accepted_locally(AB, u)
  a2 = is_accepted_locally(AB, v)
  if a1 < a2:
    ops.append("accepted")
  elif a1 > a2:
    ops.append("synonymized")
  return ops


# Could be:
#   Change of rank (promotion/demotion).
#   Change of acceptedness (accepted/synonymized).
#   Change of epithet / gender (regendered).
#   Change of genus (moved).
#   Change of concept (lumped / split / changed).
#   Possible change of concept (perhaps lumped / perhaps split)

def impute_operation_1(AB, u, v_rel, hom):
  ops = []

  v = v_rel.record if v_rel else None

  if monitor(u):
    x = get_outject(u)
    red = get_redundant(x, None)
    if red:
      log("# %s %s has red_id %s" %
          (get_primary_key(x), blurb(u),
           get_primary_key(red)))
  if monitor(v):
    y = get_outject(v)
    red = get_redundant(y, None)
    if red:
      log("# %s %s has red_id %s" %
          (get_primary_key(y), blurb(v),
           get_primary_key(red)))

  # Does not happen?  Why not?
  if u and get_redundant(get_outject(u), None):
    ops.append("redundant in A")
  if v and get_redundant(get_outject(v), None):
    ops.append("redundant in B")

  return ops

def description_for_report(AB, u):
  if u:
    x = get_outject(u)
    return (get_primary_key(x),
            "%s%s" % (get_canonical(u, None) or get_scientific(u, None),
                      "" if is_accepted_locally(AB, u) else "*"))
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

def generate_report(AB, d_path):
  with rows.open(d_path, "w") as d_gen:
    log("# a")              # doesn't get written
    # writes to io.open(...), a text stream  ??
    d_gen.write_rows(generate_plugin_report(AB))
    print("# resetting logging", file=sys.stderr)    # SUCCESS
    reset_log_allowance()
    print("# wrote rows 2", file=sys.stderr) # SUCCESS
    log("# wrote rows 2a")                   # FAILURE

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
      log("# theorize")         # gets written.
      theory.theorize(AB, False)
      generate_report(AB, d_path)
      log("# c")              # does not gets written
      op_counts_report(counts)  # counts is global.  prints to stderr.
