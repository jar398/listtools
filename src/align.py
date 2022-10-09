#!/usr/bin/env python3

import sys, csv, argparse
import util, property as prop
import checklist, workspace, match_records
import theory, span

from util import windex, MISSING
from property import mep_get, mep_set
from rcc5 import *
from checklist import *
from workspace import *
from theory import get_tipward, local_sup

def align(A_iter, B_iter):
  A = rows_to_checklist(A_iter, {'tag': "A"})  # meta
  B = rows_to_checklist(B_iter, {'tag': "B"})  # meta
  AB = workspace.make_workspace(A, B, {'tag': "AB"})
  al = make_alignment(AB)
  return generate_alignment_report(al)

# Returns a row generator

def generate_alignment_report(al):
  yield ("A id", "A name", "rcc5", "B name", "B id",
         "action", "note", "comment")
  for art in al:
    (v, ship, w, note, comment) = art
    assert note
    x = get_outject(v)
    y = get_outject(w)
    yield  (get_primary_key(x),
            blurb(x),
            rcc5_symbol(ship),
            blurb(y),
            get_primary_key(y),
            describe(v, ship, w),
            note,
            comment)

def describe(v, ship, w):
  mat = (get_matched(v) == w)
  if ship == GT:
    action = 'split' if mat else 'add'
  elif ship == EQ:
    action = '' if mat else 'equivalent'
  elif ship == LT:
    action = 'lump' if mat else 'remove'
  elif ship == DISJOINT:
    action = 'move' if mat else 'disjoint'
  elif ship == CONFLICT:        # Euler/X calls it 'overlaps'
    action = 'reorganize' if mat else 'overlaps'
  # Non-RCC5
  elif ship == OVERLAP:
    action = 'intersect in some way'
  elif ship == NOINFO:
    action = 'no information'
  elif ship == SYNONYM:
    action = 'synonym of'
  else:
    action = rcc5_symbol(ship)
  return action

# -----------------------------------------------------------------------------
# An alignment is a generator (an Iterable) of articulations.
# An articulation is a tuple (v, ship, w, note, comment).

def make_alignment(AB):
  matches_iter = match_records.match_records(checklist_to_rows(AB.A),
                                             checklist_to_rows(AB.B))
  theory.load_matches(matches_iter, AB)
  theory.theorize(AB)
  span.span(AB)

  seen = {}

  # Preorder
  def traverse(z, infra):
    for art in taxon_articulations(AB, z, infra):
      (v, ship, w, note, comment) = art
      if theory.isinA(AB, v):
        key = (get_primary_key(v), get_primary_key(w))
      else:
        key = (get_primary_key(w), get_primary_key(v))
      if key in seen or (v == AB.top or w == AB.top):
        pass
      else:
        yield art
        seen[key] = True
    if get_rank(z, None) == 'species': infra = z
    for c in get_inferiors(z):
      yield from traverse(c, infra)
  yield from traverse(AB.top, False)

# *** METHOD 1 ***
# Articulations of interest:
#    w (the cross_relation)
#    m (the name match)

def taxon_articulations(AB, z, infra):
  if not is_species(z): return

  if theory.isinA(AB, z) and theory.get_equivalent(z, None):
    # Redundant, generated previously
    return

  # Topology-based partner
  rel = partner(AB, z)
  w = rel.record
  (w, comment) = get_acceptable(AB, w) # Accepted and not infraspecific
  yield articulate(AB, z, w, comment)

  # Match-based partner
  m = get_match(z, None)      # name match

  if m and m.record:
    yield articulate(AB, z, m.record, comment)

# Is this the place to filter out redundancies and unwanteds?

def articulate(AB, v, w, comment):
  (v, comment1) = get_acceptable(AB, v)
  (w, comment2) = get_acceptable(AB, w)
  if NORMALIZE and theory.isinB(AB, v):   # Normalize order?
    (v, w) = (w, v)
    comment = rev_comment(comment)
  rel = theory.cross_relation(AB, v, w)
  return (v, rel.relationship, w, rel.note, comment)

NORMALIZE = True

def rev_comment(comment):
  if comment:
    return "rev(%s)" % comment
  else:
    return comment

# Partner species.  Returns (record, comment)

def partner(AB, v):
  assert isinstance(v, prop.Record)
  equ = theory.get_equivalent(v, None)
  if equ:
    return equ
  else:
    xsup = theory.cross_superior(AB, v)
    if xsup:
      return xsup
    return None

def is_species(x):
  return get_rank(x, None) == 'species' and is_accepted(x)

# Not used !?  Ascend the hierarchy until species is found

def get_species(AB, x):
  if is_species(x):
    return x
  sup = local_sup(AB, x, None)
  if sup == None:
    return None
  else:
    return get_species(AB, sup.record)

# Acceptable = accepted and not infraspecific.
# Approximation!

def is_acceptable(y):
  a = get_outject(y)
  if (is_accepted(a) and
      get_rank(a, None) != 'subspecies'):
    return True
  else:
    return False

def get_acceptable(AB, w):
  comment = ""
  a = get_outject(w)
  b = get_accepted(a)
  if b and b != a:
    a = b
    comment = "via synonym"
  if get_rank(a, None) == 'subspecies':
    sup = get_superior(a, None)
    a = sup.record
    comment = "via subspecies"
  if theory.isinA(AB, w):
    z = AB.in_left(a)
  else:
    z = AB.in_right(a)
  return (z, comment)

def sort_key(c):
  return (get_canonical(c, ''),
          get_scientific(c, ''),
          get_primary_key(c, ''))

# -----------------------------------------------------------------------------
# Requires spanning tree

def specimens_table(AB):
  yield("checklist", "species_id", "species_canonical", "specimen_id", "specimen_canonical")

  def traverse(z):
    if is_species(z):
      # The species that the specimen belongs to
      tag = get_tag(AB.A) if theory.isinA(AB, z) else get_tag(AB.B)
      x = get_outject(z)

      for r in theory.specimen_B_records(AB, z):
        y = get_outject(r)
        yield(tag,
              get_primary_key(x),
              blurb(x),
              get_primary_key(y),
              "^%s" % blurb(y))
    for c in get_inferiors(z):
      yield from traverse(c)
  yield from traverse(AB.top)

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

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    """)
  parser.add_argument('--A', help="the A checklist, as path name or -",
                      default='-')
  parser.add_argument('--B', help="the B checklist, as path name or -",
                      default='-')
  parser.add_argument('--test', action='store_true', help="run smoke tests")
  args=parser.parse_args()
  if args.test:
    test()
  else:
    a_path = args.A
    b_path = args.B
    assert a_path != b_path
    with util.stdopen(a_path) as a_file:
      with util.stdopen(b_path) as b_file:
        report = align(csv.reader(a_file),
                       csv.reader(b_file))
        util.write_rows(report, sys.stdout)
