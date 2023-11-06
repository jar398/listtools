#!/usr/bin/env python3

import sys, csv, argparse
import util, property as prop
import checklist, workspace, exemplar #match_records
import estimate, theory, span, rows, linkage

from util import windex, MISSING
from property import mep_get, mep_set
from rcc5 import *
from checklist import *
from workspace import *
from theory import cross_compare, is_species
from estimate import get_estimate, get_equivalent

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
  def normalize(art):
    (u, rel, kind) = art
    (u, rel) = normalize_articulation((u, rel))
    return (u, rel, kind)

  return sort_alignment(map(normalize, generate_alignment(AB)))

def sort_alignment(al):
  unsorted = list(al)
  log("-- %s articulations in alignment" % len(unsorted))
  unsorted.sort(key=articulation_sort_order)
  return unsorted

# B (rel) has sort priority.  u in A part of AB, v in B part
def articulation_sort_order(art):
  (u, rel, kind) = art
  if rel:
    v = rel.record
    return (blurb(v), blurb(u), kind)
  else:
    return (blurb(u), '~~~~~', kind)

def articulation(z, rel):
  assert get_workspace(z)
  return (z, rel)

def get_acceptables(AB, u, v):
  return (get_acceptable(AB, u), get_acceptable(AB, v))

# Return generator for alignment; each row is an articulation

def generate_alignment(AB, matches=None):
  # Assumes that name matches are already stored in AB.
  theory.theorize(AB)

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

    def traverse(x):            # u in AB
      u = AB.in_left(x)
      assert get_workspace(u)
      #if not get_acceptable(AB, u): return   # ?????

      # Show species overlaps
      e = exemplar.get_bare_exemplar(u)
      if e:
        (_, ue, ve) = e
        if ue is u:
          v1 = ve if in_same_tree(AB, u, ue) else ue
          (u2, v2) = get_acceptables(AB, u, v1)
          if unseen(u2, v2):
            rel2 = cross_compare(AB, u2, v2)
            yield (u2, rel2, 'related')

      if is_acceptable_locally(AB, u):

        est_rel = get_estimate(u)
        v = est_rel.record

        if is_acceptable_locally(AB, est_rel.record):
          # Equivalence - for long report
          if unseen(u, v):
            if est_rel.relationship == EQ:
                yield (u, est_rel, 'equivalent')
            else:
              # De novo, unassigned split, or retracted
              if get_rank(u, None) == 'species':
                if not estimate.get_cross_mrca(u, None):
                  yield (u,
                         relation(est_rel.relationship | DISJOINT,
                                  v, est_rel.note, est_rel.span),
                         'peripheral')

          # 'Insertions' / change of parent   vc <? u < v
          # u < v    v is smallest such
          for c in get_children(get_outject(v), ()):
            # c < v
            if is_acceptable(c):
              vc = AB.in_right(c)
              # vc < v
              if unseen(vc, u):
                rel = cross_compare(AB, vc, u)
                ship = rel.relationship
                assert not ship == EQ
                if ship != DISJOINT and ship != GT:  # LT, OVERLAP
                  yield (vc, rel, 'topological change')

        # Report on change of circumscription
        vl = get_link(u, None)         # Match by name
        if vl == False:
          yield (u, relation(NOINFO, vl), 'ambiguous')
        elif vl and is_acceptable_locally(AB, vl):
          if unseen(u, vl):
            rel = cross_compare(AB, u, vl)
            if rel.relationship != EQ:
              yield (u, rel, 'change in circumscription')
            elif get_canonical(u) != get_canonical(vl):
              yield (u, rel, 'rename')

      for c in get_inferiors(x): # subspecies and synonyms
        yield from traverse(c)

    yield from traverse(AB.A.top)

  yield from doit(swap(AB), True) # B has priority
  yield from doit(AB, False)

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
  assert get_workspace(u)
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
  return get_best_name(c)

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
        # linkage.find_links(AB)  -- not needed!
        al = generate_alignment(AB)
        report = simple_report(al)
        util.write_rows(report, sys.stdout)
