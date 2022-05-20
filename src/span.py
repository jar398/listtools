#!/usr/bin/env python3

# Create a spanning tree for A+B based on B priority

# For each equivalent pair, only the B record will have a superior and inferiors.

import theory
from checklist import rows_to_checklist, checklist_to_rows, get_superior, \
  set_superior, relation, get_children, set_children, get_synonyms, set_synonyms
import newick
import checklist, match_records, workspace
from rcc5 import *
from util import log

def span(AB):
  def traverse(AB, B_priority):
    for y in checklist.preorder_records(AB.B):
      w = AB.in_right(y)
      if B_priority or not theory.get_equivalent(y, None):
        p_rel = possible_insertion(AB, y, B_priority)
        if p_rel:
          p = p_rel.record      # in A
          z = AB.in_left(p)
          sup = relation(p_rel.relationship, z, p_rel.status, p_rel.note)
          plug_superior(w, sup, B_priority)
        else:
          q_rel = get_superior(y, None)
          if q_rel:
            z = AB.in_right(q_rel.record)
            sup = relation(q_rel.relationship, z, q_rel.status, q_rel.note)
            plug_superior(w, sup, B_priority)
  traverse(AB, True)
  traverse(theory.swap(AB), False)
  # checklist.ensure_inferiors_indexed(AB)  no...

def plug_superior(w, sup, B_priority):
  assert len(checklist.get_source_name(w)) != 1
  assert len(checklist.get_source_name(sup.record)) != 1

  if B_priority or not theory.get_equivalent(w, None):
    set_superior(w, sup)
    if sup.relationship == ACCEPTED:
      ch = get_children(sup.record, None)
      if ch != None:
        ch.append(w)
      else:
        set_children(sup.record, [w])
    else:
      ch = get_synonyms(sup.record, None)
      if ch != None:
        ch.append(w)
      else:
        set_synonyms(sup.record, [w])

# y is in B (possibly A/B swapped)
# Returns ... rel to x in A

def possible_insertion(AB, y, B_priority):
  assert len(checklist.get_source_name(y)) == 1
  # Candidate superiors in A+B are {p, q} where p=get_superior(y)
  # Insertion candidate is p
  q_rel = get_superior(y, None)         # default, in B
  p = theory.cross_superior(AB, y) # in A
  if p:
    assert len(checklist.get_source_name(p)) == 1
    p_rel = relation(LT, p, "accepted", "insertion")
  else:
    p_rel = None
  if not q_rel: return p_rel # Insertion (no q)
  if not p: return None
  q = q_rel.record                # in B
  assert len(checklist.get_source_name(y)) == 1
  ship = theory.cross_relationship(p, q).relationship
  if ship == LT:
    # Insertion, y < p < q
    return p_rel
  if ship == GT or ship == EQ:
    return None            #p >= q > y
  # CONFLICT or OVERLAP case.  Priority wins.
  if B_priority:
    log("# making arbitrary parent choice (B, priority)")
    return None    #p is the default
  else:
    log("# making arbitrary parent choice (A, non-priority)")
    return p_rel

# -----------------------------------------------------------------------------

def test():
  def testit(m, n):
    log("\n-- Test: %s + %s --" % (m, n))
    B = rows_to_checklist(newick.parse_newick(n),
                          {"name": "B"})  # meta
    A = rows_to_checklist(newick.parse_newick(m),
                          {"name": "A"})  # meta
    matches_iter = match_records.match_records(checklist_to_rows(A), checklist_to_rows(B))
    AB = workspace.make_workspace(A, B, {"name": "AB"})
    theory.load_matches(matches_iter, AB)
    theory.analyze_tipwards(AB)                # also find tipes
    theory.compute_blocks(AB)
    theory.compute_cross_mrcas(AB)
    theory.find_equivalents(AB)
    span(AB)
    print(newick.compose_newick(checklist_to_rows(AB)))
  testit("a", "a")              # A + B
  testit("(c,d)a", "(c,e)b")
  testit("((a,b)e,c)d", "(a,(b,c)f)D")

if __name__ == '__main__':
  test()
