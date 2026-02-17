#!/usr/bin/env python3

import argparse
import rows
import exemplar, theory, jumble

from util import log
from checklist import preorder_records, get_primary_key, \
  get_canonical, get_scientific, get_parts, \
  reverse_relation, blurb
from workspace import ingest_workspace, is_accepted_locally, isinA, isinB, \
  separated, get_outject
from rcc5 import *

from specimen import same_type_ufs, sid_to_epithet, get_homotypic
from estimate import get_venn, get_estimate, get_block
from theory import is_species, get_species, compare

counts = {}
def count(tag):
  if tag in counts: counts[tag] += 1
  else: counts[tag] = 1


# Generator for csv file rows

def analyze(AB):
  yield emit_articulation(True)
  seen = {}
  # info = (z, rel, tribe).  Generrate one row of csv file
  for info in generate_articulations(AB):
    z = info[0]
    rel = info[1]
    if isinA(AB, z):
      key = (get_primary_key(z),
             get_primary_key(rel.record))
    else:
      key = (get_primary_key(rel.record),
             get_primary_key(z))
    if key in seen:
      pass
    else:
      seen[key] = True
      yield emit_articulation(info)
  for (op, op_count) in counts.items():
    log("  %6d %s" % (op_count, op))

# Generator for articulations ('infos')

def generate_articulations(AB):
  # Set parent pointers / ensure no duplication
  for z in preorder_records(jumble.jumble_workspace(AB)):

    if is_species(z):
      # get_venn(u, v_rel)[1]  = intersecting specimens ... ?
      inters = theory.get_intersecting_species(z)  # list

      # 'Species' carried over from one checklist to the next.
      # I.e. if z in AB is from A, then hom will be in AB but from B.
      hom = get_homotypic(AB, z)
      if hom:
        near = get_species(hom)
      else:
        near = best_topo_match(inters, z)
        # near is often missing

      # Sort intersectors: first element should be 'near' z, i.e. same type or
      # greatest overlap with z;
      # other intersectors should sort by epithet.
      def inter_sort_key(inter):
        return (least_epithet_key
                if inter is near
                else epithet_key(inter))
      inters = sorted(inters,
                      key=inter_sort_key)

      # Classify z according to relations to intersectors
      relations = list(map(lambda i: compare(AB, z, i), inters))
      # get_venn(u, rel)  ???
      tribe = impute_tribe(z, relations)
      assert tribe
      for rel in relations:
        # Generate one articulation for each intersector!
        yield (z, rel, tribe)
      count(tribe)

# info is either True, for header row, or (z, rel, tribe)

def emit_articulation(info):
  if info == True:
    u_fields = ("A concept taxonID", "A concept name")
    rcc5_field = "RCC5"         # header
    v_fields = ("B concept taxonID", "B concept name")
    venn_fields = ("in both A and B concepts",
                   "in A concept but not in B concept",
                   "in B concept but not in A concept")
    op_field = "mode"
  else:
    (z, rel, tribe) = info
    if isinA(AB, z):
      u = z
      if rel: v = rel.record
    else:
      v = z
      if rel:
        u = rel.record
        rel = reverse_relation(u, rel)
    u_fields = description_for_report(AB, u)
    v_fields = description_for_report(AB, v)
    if rel:
      rcc5_field = rcc5_symbol(rel.relationship)
      if rcc5_field.startswith('='):
        # Kludge to appease spreadsheet programs
        rcc5_field = "'" + rcc5_field
    else:
      rcc5_field = MISSING

    ops = [tribe]
    # TBD: Count instances of each tribe!

    # Venn diagram columns
    if u and v:
      (u_and_v, u_not_v, v_not_u) = get_venn(u, rel)
      venn_fields = (show_specimen_id_set(AB, u_and_v),
                     show_specimen_id_set(AB, u_not_v),
                     show_specimen_id_set(AB, v_not_u))
    else:
      venn_fields = (MISSING, MISSING, MISSING)

    # ops += impute_concept_change(AB, u, v_rel, homotypic)
    # ops += impute_name_change(AB, u, v_rel, homotypic)

    op_field = "; ".join(ops)

  return (op_field,
          u_fields[0], u_fields[1],
          rcc5_field,
          v_fields[0], v_fields[1],
          venn_fields[0], venn_fields[1], venn_fields[2],
          )


# split an A species, coalesce a B species, carry over a congruent
# species, A or B species does something 'other' than those

def impute_tribe(z, relations):
  if isinA(AB, z):
    # z in AB comes from A
    if len(relations) == 0:
      tribe = "gone"
    elif len(relations) == 1:
      rel = relations[0]
      if rel.relationship == EQ:
        tribe = "congruent in A"
      else:
        tribe = "pawn in A"
    else:
      # An A species can split
      # e.g. Bradypus torquata
      gt = True
      for rel in relations:
        if rel.relationship != GT and rel.relationship != GE:
          gt = False
          break
      if gt:
        tribe = "split"
      else:
        tribe = "other in A"
  else:
    # z in AB comes from B
    if len(relations) == 0:
      tribe = "new"
    elif len(relations) == 1:
      rel = relations[0]
      if rel.relationship == EQ:
        tribe = "congruent in B"
      else:
        tribe = "pawn in B"
    else:
      # A B species can lump
      lt = True
      for rel in relations:
        if rel.relationship != GT and rel.relationship != GE:
          lt = False
          break
      if lt:
        tribe = "lump"
      else:
        tribe = "other in B"
  return tribe

least_epithet_key = ("aaa", "aaa")

def epithet_key(x):
  p = get_parts(x)
  return (p.epithet, p.genus)

# Used to figure out which intersector to put first.
# Should be either contypic or congruent.

def best_topo_match(inters, u):
  if len(inters) == 0:
    return None                 # No intersecting species
  else:
    u_ids = get_block(u)        # set of exemplar ids
    def intensity(v):
      v_ids = get_block(v)    # Intersectors.
      meet = len(u_ids & v_ids)
      # Ad hoc rule:
      # Prefer: 1. maximum overlap, 2. minimum nonoverlap
      return (-meet, len(u_ids | v_ids) - meet)
    return min(inters, key=intensity)    # Could be congruent

# Given a specimen(exemplar) id set, produce a human-readable string
# that can be put in a column of the report

def show_specimen_id_set(AB, id_set):
  # xep = (specimen_id, epithet)
  # Sort by epithet
  xeps = sorted(map(lambda sid:(sid, sid_to_epithet(AB, sid)), id_set),
                key=lambda xep: xep[1])
  return "; ".join(map(lambda xep:"%s %s" % xep, xeps))

def description_for_report(AB, u):
  if u:
    x = get_outject(u)
    return (get_primary_key(x),
            "%s%s" % (get_canonical(u, None) or get_scientific(u, None),
                      "" if is_accepted_locally(AB, u) else "*"))
  else:
    return (MISSING, MISSING)

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
        gen = (row for row in analyze(AB))
        kludge = list(gen)
        log("# %s rows" % len(kludge))
        d_gen.write_rows(kludge)
