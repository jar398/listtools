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
from theory import local_sup, get_equivalent

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

# Return generator for alignment; each row is an articulation

def generate_alignment(AB, matches=None):
  # Assumes that name matches are already stored in AB.
  theory.theorize(AB)
  #span.span(AB)

  seen = {}

  def doit(AB, swapped):
    def traverse(x):            # u in AB
      assert is_accepted(x)
      count = [0, -1]
      u = AB.in_left(x)
      assert is_accepted_locally(AB, u)
      arts = list(taxon_articulators(AB, u))
      assert len(arts) > 0, blurb(u)
      for art in arts: # Normalized
        v = get_acceptable(AB, art)        # Compare v to art
        assert separated(u, v)
        assert is_accepted_locally(AB, v)
        if swapped: (u, v) = (v, u)
        count[0] += 1
        count[1] = max(count[0], count[1])
        if monitor(u) or monitor(v):
          if count[0] % 100 == 0:
            log("# articulator %s: %s, %s" % (blurb(art), blurb(v), count[1]))
        key = (get_primary_key(u), get_primary_key(v))
        if key in seen:
          pass
        else:
          seen[key] = True
          assert separated(u, v)
          rel = theory.cross_compare(AB, u, v)

          # Suppress if entailed by either checklist
          ship = rel.relationship
          if obvious(AB, u, rel.relationship, v):
            pass
          else:
            yield (u, rel)
      if not is_species(u):
        for c in get_children(x, ()): # subspecies and synonyms
          yield from traverse(c)
    yield from traverse(AB.A.top) # B has priority

  yield from doit(swap(AB), True) # B has priority
  yield from doit(AB, False)

def obvious(AB, u, ship, v):
  if is_toplike(u) or is_toplike(v):
    return True
  if ((get_equivalent(AB, u) or get_equivalent(AB, v)) and
      (ship == LT or ship == GT or ship == DISJOINT)):
    return True

# Returns generator of records v to compare to u.
# u and v are in different checklists.
# There may be duplicates.
# Need to filter out articulations involving synonyms.

def taxon_articulators(AB, u):
  assert isinA(AB, u)
  assert get_workspace(u)
  rel = theory.get_estimate(u, None)
  assert get_workspace(rel.record)
  assert is_accepted(get_outject(u))
  if rel and rel.record:
    v = rel.record
    assert separated(u, v)
    assert is_accepted_locally(AB, v)
    yield v                   # u <= v
    if True or rel.relationship != EQ:
      for c in get_children(get_outject(v), ()):
        assert is_accepted(c)
        v2 = AB.in_right(c)
        if not get_equivalent(AB, v2):
          yield v2
  v = get_link(u, None)         # Match by name
  if v:
    assert is_accepted_locally(AB, v)
    yield v

def is_species(v):              # z local
  z = get_outject(v)
  return get_rank(z, None) == 'species' and is_accepted(z)

# Acceptable = accepted and not infraspecific.
# Approximation!

def is_acceptable(u):           # in AB
  return is_acceptable_locally(get_outject(u))

# Returns ancestor

def get_acceptable(AB, u):
  x = get_outject(u)
  # Ascend until an acceptable is found
  while not is_acceptable_locally(x):
    sup = get_superior(x, None)
    if not sup:
      break
    x = sup.record
  return AB.in_left(x) if isinA(AB, u) else AB.in_right(x)

# This is inaccurate - doesn't allow nested infraspecific ranks

def is_acceptable_locally(x):
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
