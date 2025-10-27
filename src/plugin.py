#!/usr/bin/env python3

import sys, csv, argparse
import rcc5, rows, checklist, workspace
import theory, exemplar, estimate
import jumble, util

from workspace import ingest_workspace, is_accepted_locally, local_sup, \
  local_accepted, isinA, isinB
from theory import is_species
from checklist import *
from rcc5 import *
from specimen import same_type_ufs, sid_to_epithet
from estimate import get_block, is_empty_block, get_estimate, get_equivalent
from estimate import get_congruent
from property import mep, mep_get, mep_set
from typify import explain
from ranks import ranks_dict

counts = {}
def count(tag):
  if tag in counts: counts[tag] += 1
  else: counts[tag] = 1

# Preorder traversal of AB.A

def generate_relators(AB):
  # For a node z in AB
  # if is_species(z) then
  #   list of species s that overlap z, sorted by epithet stem
  #   if none then make fake relationship (z, None) or (None, z)
  #   form relationships R = (z, s) or (s, z)
  #   suppress already-seen R
  #   emit relationships R sorted by ... highest overlap then epithet?
  # otherwise:
  #   sort children of Z by (A vs. B) then by name
  # recur
  seen = {}
  def half_key(node):
    if node:
      return get_specimen_id(get_type_uf(node))
    else:
      return None
  def sort_key(node):
    p = get_parts(node)
    return (p.epithet, p.genus)
  def maybe_flip(pair):
    (z, s) = pair
    if isinB(AB, z): return (s, z)
    else: return pair
  def recur(z):
    if is_species(z):
      inters = theory.get_intersecting_species(z)
      inters = sort_list(inters, key=sort_key)
      if len(inters) == 0:
        if isinA(z):
          yield (z, None)
        else:
          yield (None, z)
      else:
        # Modify this sort so that types come up at top
        relators = map(lambda s: (z, s), inters)
        for relator in relators:
          relator = maybe_flip(relator)
          (z, s) = relator
          key = (half_key(z), half_key(s))
          if not key in seen:
            seen.add(key)
            yield relator
    else:
      #   sort children of Z by (A vs. B) then by name
      childs = get_children(z)
      childs = sort_list(childs, key=sort_key)
      for c in childs:
        yield from recur(c)
      
  yield from recur(AB.top)
  return None

def generate_plugin_report(AB):
  yield generate_row(AB, True, True, True) # header
  i = [0]
  frequency = 5000

  jumble.jumble_workspace(AB)
  spa = spb_spe = spb_nspe = spb_ne = 0

  for z in preorder_records(AB):
    count("records in A+B")
    i[0] += 1
    if i[0] % frequency == 0:
      log("# %s %s" % (i[0], blurb(z)))
    ops = []

    # A -
    # - B
    # A A
    # A B
    
    if is_species(z):

      if isinB(AB, z) and get_congruent(AB, z):
        # Skip records that have equivalents in A - they
        # are already handled
        continue

      count("species")          # ? what to call this ?

      w = choose_partner(AB, z) # may be None
      if isinA(AB, z):
        u = z; v = w
        spa += 1
      else:
        u = w; v = z
        spb_ne += 1
        ops.append("no corresponding species")

      yield generate_row(AB, u, v, ops)

  # Kludge!!!  To make sure summary gets displayed
  util.log_allowance = 1000

  for (op, op_count) in counts.items():
    log("  %6d %s" % (op_count, op))

# For each concept, we pick one articulation.  The articulation is
# between the species and a "partner".

def choose_partner(AB, u):
  assert is_species(u)
  e_rel = get_equivalent(AB, u) # not just congruent.  ??
  if e_rel:
    sp = theory.get_species(e_rel.record)
    return sp or e_rel.record

  # ... pick the 'best' one from intersectors ...
  inters = theory.get_intersecting_species(u)
  if len(inters) == 0:
    return None                 # No intersecting species
  else:
    u_ids = get_block(u)
    def intensity(v):
      v_ids = get_block(v)
      hom = 0 if same_type_ufs(u, v) else 1
      meet = len(u_ids & v_ids)
      # Ad hoc rule:
      # Prefer: 1. homotypic, 2. maximum overlap, 3. minimum nonoverlap
      return (hom, -meet, len(u_ids | v_ids) - meet)
    return min(inters, key=intensity)    # Could be congruent


# Return either the header row or a data row; the code for both is
# together to make sure they stay in sync.

def generate_row(AB, u, v, ops):
  # A - rcc5 - B columns
  if u == True:
    u_fields = ("A concept taxonID", "A concept name")
    rcc5_field = "RCC5"         # header
    v_fields = ("B concept taxonID", "B concept name")
  else:
    u_fields = description_for_report(AB, u)
    v_fields = description_for_report(AB, v)
    if u and v:
      v_rel = theory.compare(AB, u, v)
      rcc5_field = rcc5_symbol(v_rel.relationship)
      if rcc5_field.startswith('='):
        # Kludge to appease spreadsheet programs
        rcc5_field = "'" + rcc5_field
    else:
      rcc5_field = MISSING
      v_rel = None

  # change column
  if u == True:
    op_field = "change"      # header
  else:
    # `homotypic` is true iff u and v are homotypic.
    homotypic = u and v and same_type_ufs(u, v)

    ops += impute_concept_change(AB, u, v_rel, homotypic)
    ops += impute_name_change(AB, u, v_rel, homotypic)
    op_field = "; ".join(ops)
    for op in ops:
      if op in counts:
        counts[op] += 1
      else:
        counts[op] = 1

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

  return (op_field,
          u_fields[0], u_fields[1],
          rcc5_field,
          v_fields[0], v_fields[1],
          u_and_v,
          u_not_v,
          v_not_u,
          )

# s is a set of exemplar ids...

def show_sid_set(AB, block):
  # xep = (sid, epithet)
  # Sort by epithet
  xeps = sorted(map(lambda sid:(sid, sid_to_epithet(AB, sid)), block),
                key=lambda xep: xep[1])
  return "; ".join(map(lambda xep:"%s %s" % xep, xeps))

def impute_concept_change(AB, u, v_rel, homotypic):
  if not u:
    op = "new"
  elif not v_rel:
    op = "gone"
  else:
    ship = v_rel.relationship
    if ship == EQ:
      # This can happen if e.g. the type is ambiguous or otherwise unmatched
      op = "congruent" if homotypic else "congruent but heterotypic"
    elif ship == LT:
      op = "expand" if homotypic else "lump"
    elif ship == GT:
      op = "contract" if homotypic else "split"
    elif ship == OVERLAP:
      op = "change" if homotypic else "reorganize"
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
      ops.append("moved")
    if p1.epithet != p2.epithet:
      ops.append("epithet")
    if p1.token != p2.token:
      ops.append("author")
    if p1.year != p2.year:
      ops.append("year")
  # maybe: promotion, demotion, rank change
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
        gen = (row for row in generate_plugin_report(AB))
        kludge = list(gen)
        log("# %s rows" % len(kludge))
        d_gen.write_rows(kludge)
