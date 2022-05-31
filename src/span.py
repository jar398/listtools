#!/usr/bin/env python3

# Create a spanning tree for A+B based on B priority

# For each equivalent pair, only the B record will have a superior and inferiors.

# Would it be worthwhile to have a way to write a spanning tree as a
# checklist, so that it can be re-read without loss? ...
# Actually, would it be worthwhile to be able to do that with a workspace ...

import theory
from checklist import rows_to_checklist, checklist_to_rows, get_superior, \
  set_superior, relation, get_children, set_children, get_synonyms, set_synonyms, \
  get_source, get_outject
import newick
import checklist, match_records, workspace
from rcc5 import *
from util import log

def span(AB):
  assert theory.get_equivalent(AB.in_left(AB.A.top)).record == AB.in_right(AB.B.top)
  AB.top = AB.in_right(AB.B.top)
  def traverse(AB, B_priority):
    for y in checklist.preorder_records(AB.B):    # starts with top
      if B_priority or not theory.get_equivalent(y, None):
        w = AB.in_right(y)
        p_rel = possible_insertion(AB, w, B_priority)
        if p_rel:
          z = p_rel.record      # in A
          sup = relation(p_rel.relationship, z, p_rel.status, p_rel.note)
          plug_superior(w, sup, B_priority)
        else:
          q_rel = get_superior(y, None)
          if q_rel:
            z = AB.in_right(q_rel.record)
            sup = relation(q_rel.relationship, z, q_rel.status, q_rel.note)
            plug_superior(w, sup, B_priority)

  AB.top = AB.in_right(AB.B.top)

  traverse(AB, True)
  traverse(theory.swap(AB), False)
  AB.indexed = True
  return AB

def plug_superior(w, sup, B_priority):
  assert len(checklist.get_source_name(w)) != 1
  assert len(checklist.get_source_name(sup.record)) != 1

  if B_priority or not theory.get_equivalent(w, None):
    checklist.link_superior(w, sup)

# y is in B (possibly A/B swapped)
# Returns Relation to p in AB if insertion, or None

def possible_insertion(AB, w, B_priority):
  y = get_outject(w)
  # Candidate superiors in A+B are {p, q} where p=get_superior(y)
  # Insertion candidate is p
  q_rel = get_superior(y, None)         # default, in B
  p = theory.cross_superior(AB, w)      # AB -> AB
  if p:
    p_rel = relation(LT, p, "accepted", "insertion")
  else:
    p_rel = None
  if not q_rel: return p_rel # Insertion (no q)
  if not p: return None
  q = q_rel.record                # in B
  ship = theory.cross_relationship(AB, p, AB.in_right(q))
  if ship == LT:
    # Insertion, y < p < q
    return p_rel
  if ship == GT or ship == EQ:
    return None            #p >= q > y
  # CONFLICT or OVERLAP case.  Priority wins.
  if B_priority:
    #log("# making arbitrary parent choice (B, priority)")
    return None    #p is the default
  else:
    #log("# confirming arbitrary parent choice?")
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
    theory.find_equivalents(AB)
    theory.compute_cross_mrcas(AB)
    span(AB)
    print(newick.compose_newick(checklist.preorder_rows(AB)))
  testit("a", "a")              # A + B
  testit("(c,d)a", "(c,e)b")
  testit("((a,b)e,c)d", "(a,(b,c)f)D")

if __name__ == '__main__':
  test()
