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
  return generate_alignment_report(AB, al)

# Returns a row generator

def generate_alignment_report(AB, al):
  yield ("A id", "A name", "rcc5", "B name", "B id",
         "action", "note", "comment")
  for art in al:
    (v, ship, w, note, comment) = art
    assert note
    if theory.isinB(AB, v):   # Normalize order?
      (v, w) = (w, v)
      comment = rev_comment(comment)
      ship = reverse_relationship(ship)
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
  rel = get_matched(v)
  mat = rel and (rel.record == w)
  stay = (get_rank(get_accepted(get_outject(v)), None) ==
          get_rank(get_accepted(get_outject(w)), None))
  if ship == GT:
    if stay: action = ('shrink' if mat else 'split')
    else: action = 'add'
  elif ship == EQ:
    action = ('' if mat else 'equivalent')
  elif ship == LT:
    if stay: action = ('expand' if mat else 'lump')
    else: action = 'remove'
  elif ship == DISJOINT:
    action = ('false match' if mat else 'disjoint')
  elif ship == CONFLICT:        # Euler/X calls it 'overlaps'
    action = ('reorganize' if mat else 'overlaps')
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
  find_acceptables(AB)

  seen = {}

  # Preorder
  def traverse(z, infra):
    if is_species(z):       # Policy
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
      infra = z
    for c in get_inferiors(z):
      # sort ??
      yield from traverse(c, infra)
  yield from traverse(AB.top, False)

# *** METHOD 1 ***
# Articulations of interest:
#    w (the cross_relation)
#    m (the name match)

OLD_METHOD = False

def taxon_articulations(AB, z, infra):
  rr = list(theory.specimen_records(AB, z))
  if len(rr) > 0:
    for (a, b) in rr:
      assert theory.separated(a, b)
      (s, comment) = get_acceptable(AB, a if theory.isinB(AB, z) else b)
      assert theory.separated(z, s)
      if b != z and b != s:
        comment = '%s; via %s' % (comment, blurb(b),)
      yield articulate(AB, z, s, comment)
  else:
    rel = theory.get_reflection(AB, z)
    (w, comment) = get_acceptable(AB, rel.record) # Accepted and not infraspecific
    yield articulate(AB, z, w, comment)

# Is this the place to filter out redundancies and unwanteds?

def articulate(AB, v, w, comment):
  (v, comment1) = get_acceptable(AB, v)
  (w, comment2) = get_acceptable(AB, w)
  if NORMALIZE and theory.isinB(AB, v):   # Normalize order?
    (v, w) = (w, v)
    comment = rev_comment(comment)
  rel = theory.cross_relation(AB, v, w)
  return (v, rel.relationship, w, rel.note, comment)

NORMALIZE = False

def rev_comment(comment):
  if comment:
    return "rev(%s)" % comment
  else:
    return comment

def is_species(z):              # x in AB
  return get_rank(z, None) == 'species' and is_accepted(get_outject(z))

# Acceptable = accepted and not infraspecific.
# Approximation!

(see_acceptable, set_acceptable) = prop.get_set(prop.declare_property('acceptable'))

def is_acceptable(y):           # in AB
  return not see_acceptable(y, None)

def get_acceptable(AB, y):
  return see_acceptable(y, None) or (y, 'acceptable')

# Acceptable = accepted AND not infraspecific
# A synonym could be a synonym of a genus, species, subspecies, etc.

def find_acceptables(AB):
  def find(C, in_lr):
    def traverse(x, a):
      if is_accepted(x):
        # Infraspecific if a is a species
        if a and get_rank(a, None) == 'species':
          # assert get_rank(x, None) != 'variety', 1
          set_acceptable(in_lr(x), (in_lr(a), 'infraspecific'))
        else:
          a = x
      else:
        set_acceptable(in_lr(x), (in_lr(a), 'synonym'))
      for c in get_inferiors(x):
        traverse(c, a)
    traverse(C.top, None)
  find(AB.A, lambda x: AB.in_left(x))
  find(AB.B, lambda y: AB.in_right(y))

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

      for (r, s) in theory.specimen_records(AB, z):
        y = get_outject(s)
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
