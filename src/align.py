#!/usr/bin/env python3

import sys, csv, argparse
import util, property as prop
import checklist, workspace, match_records
import theory, span, rows, linkage

from util import windex, MISSING
from property import mep_get, mep_set
from rcc5 import *
from checklist import *
from workspace import *
from theory import local_sup, \
  get_block, is_empty_block
from lub import get_equivalent, get_estimate

def align(A_iter, B_iter, A_name='A', B_name='B', matches_iter=None):
  AB = ingest_workspace(a_rows.rows(), b_rows.rows(),
                        A_name=A_name, B_name=B_name)
  find_links(AB, matches_iter)
  al = generate_alignment(AB)
  return generate_alignment_report(al, AB)

# -----------------------------------------------------------------------------
# An alignment is a set of species/species articulations.
# An articulation is a pair (record, relation).

def make_alignment(AB):
  unsorted = list(generate_alignment(AB))
  log("-- %s articulations in alignment" % len(unsorted))
  return sorted(unsorted, key=articulation_order)

# B has priority.  u in A part of AB, v in B part
def articulation_order(art):
  (u, rel) = art
  v = rel.record
  ship = rel.relationship
  if is_species(v) or not is_species(u):
    return (blurb(v), ship)
  else:
    return (blurb(u), ship)

def articulation(z, rel):
  assert get_workspace(z)
  return (z, rel)

# Return generator for alignment; each row is an articulation

def generate_alignment(AB, matches=None):
  # Assumes that name matches are already stored in AB.
  theory.theorize(AB)
  #span.span(AB)

  seen = set()

  def unseen(u, v):
    ku = u.id
    kv = v.id
    key = (kv, ku) if kv < ku else (ku, kv)
    if key in seen:
      return False
    else:
      seen.add(key)
      return True

  def doit(AB, swapped):

    def process_articulator(u, v_art, rd):
      v = get_acceptable(AB, v_art)        # Compare v to art
      if unseen(u, v):
        rel = theory.cross_compare(AB, u, v)
        if (is_empty_block(get_block(u)) and
            not is_empty_block(get_block(v))): # Peripheral attachment
          rel = relation(rel.relationship | DISJOINT,    # Nico's request
                         rel.record,
                         rel.note,
                         rel.span)
        if rd or rel.relationship == OVERLAP: # can't happen
          assert get_workspace(u)
          assert get_workspace(rel.record)
          yield (u, rel)

    def traverse(x):            # u in AB
      u = AB.in_left(x)
      assert get_workspace(u)
      if not get_acceptable(AB, u): # ?????
        return

      # 1. Report the estimate (= or <)
      v = theory.get_estimate(u, None).record
      yield from process_articulator(u, v, 1)

      # 2. Report estimate vs. parent - may overlap
      if False:
       sup = local_sup(AB, u)
       if sup:
        p = sup.record
        assert get_workspace(p)
        yield from process_articulator(p, v, 2)

      # 3. Report on children of estimate looking for ><
      for c in get_inferiors(get_outject(v)):
        v2 = AB.in_right(c)
        if not get_equivalent(AB, v2): # premature optimization
          # compare(AB, u, v2).relationship == OVERLAP ...
          yield from process_articulator(u, v2, False)

      # 4. Report on name match
      vl = get_link(u, False)         # Match by name
      if vl: yield from process_articulator(u, vl, 4)

      if not is_species(u):
        for c in get_children(x, ()): # subspecies and synonyms
          yield from traverse(c)
    yield from traverse(AB.A.top)

  yield from doit(swap(AB), True) # B has priority
  yield from doit(AB, False)

def obvious(AB, u, ship, v):
  if is_toplike(u) or is_toplike(v):
    return True
  if ((get_equivalent(AB, u) or get_equivalent(AB, v)) and
      (ship == LT or ship == GT or ship == DISJOINT)):
    return True

def is_species(u):              # z local
  x = get_outject(u)
  return get_rank(x, None) == 'species' and is_accepted(x)

# Returns ancestor

def get_acceptable(AB, u):
  assert get_workspace(u)
  x = get_outject(u)
  # Ascend until an acceptable is found
  while not is_acceptable(x):
    sup = get_superior(x, None)
    if not sup:
      break
    x = sup.record
  return AB.in_left(x) if isinA(AB, u) else AB.in_right(x)

# Acceptable = accepted and not infraspecific.
# Approximation!

def is_acceptable_locally(AB, u):  # in AB
  return is_acceptable(get_outject(u))

def is_acceptable(x):
  assert not get_workspace(x)
  if is_accepted(x):
    if get_rank(x, None) == 'species':
      return True
    sup = get_superior(x, None)
    if not sup:
      return True
    if get_rank(sup.record, None) == 'species':
      # Taxa that are below species (variety, etc.) are not acceptable
      return False
    return True
  else:
    return False

def sort_key(c):
  return (get_canonical(c, ''),
          get_scientific(c, ''),
          get_primary_key(c, ''))

# -----------------------------------------------------------------------------

def test():
  def testit(m, n):
    log("\n-- Test: %s + %s --" % (m, n))
    report = align(newick.parse_newick(n),
                   newick.parse_newick(m))
    util.write_rows(report, sys.stdout)
  testit("a", "a")              # A + B
  testit("(c,d)a", "(c,e)b")
  testit("((a,b)e,c)d", "(a,(b,c)f)D")

def simple_report(al):
  def artrow(art):
    (subj, pred) = art
    return (blurb(subj),
            rcc5_symbol(pred.relationship),
            blurb(pred.record),
            pred.note)
  yield ("A record", "relationship", "B record", "note")
  for art in al:
    yield artrow(art)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    """)
  parser.add_argument('--A', help="the A checklist, as path name or -",
                      default='-')
  parser.add_argument('--B', help="the B checklist, as path name or -",
                      default='-')
  parser.add_argument('--Aname', help="short name of the A checklist",
                      default='A')
  parser.add_argument('--Bname', help="short name of the B checklist",
                      default='B')
  parser.add_argument('--test', action='store_true', help="run smoke tests")
  args=parser.parse_args()
  if args.test:
    test()
  else:
    a_path = args.A
    b_path = args.B
    assert a_path != b_path
    with rows.open(a_path) as a_rows:
      with rows.open(b_path) as b_rows:
        # compute name matches afresh
        AB = ingest_workspace(a_rows.rows(), b_rows.rows(),
                              A_name=args.Aname, B_name=args.Bname)
        linkage.find_links(AB)
        al = generate_alignment(AB)
        report = simple_report(al)
        util.write_rows(report, sys.stdout)
