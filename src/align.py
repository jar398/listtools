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
from theory import local_sup

def ingest(A_iter, A_name, B_iter, B_name):
  A = rows_to_checklist(A_iter, {'tag': A_name or "A"})  # meta
  B = rows_to_checklist(B_iter, {'tag': B_name or "B"})  # meta
  AB = workspace.make_workspace(A, B, {'tag': "AB"})
  al = make_alignment(AB)
  return (al, AB)

# Returns a row generator

def generate_alignment_report(al, AB):
  yield ("A id", "A name", "rcc5", "B name", "B id",
         "category", "note", "witnesses", "match")
  for (v, rel) in al:
    ship = rel.relationship
    w = rel.record
    x = theory.get_outject(v)
    y = theory.get_outject(w)
    d = category(v, ship, w)

    if d.startswith('remove'):
      # The following seems to have zero effect.
      mrel = get_match(x, None)          # Mutual match - we expect this to be included
      if mrel:
        if mrel.record:
          rel = theory.cross_relation(AB, v, mrel.record)
          log("# removing [%s %s %s], %s" % \
              (blurb(v), rcc5_symbol(rel.relationship),
               blurb(rel.record), rel.note))
        else:
          log("# removing %s, match note %s" %
              (blurb(v), rel.note))

    yield  (get_primary_key(x),
            blurb(x),
            rcc5_symbol(ship),
            blurb(y),
            get_primary_key(y),
            d, rel.note,
            witness_comment(AB, v, rel),
            match_comment(AB, v, rel))

def generate_short_report(al, AB):
  yield ("A name", "B name", "category", "witnesses", "match")
  for (v, rel) in al:
    w = rel.record
    ship = rel.relationship
    assert is_acceptable(v), blurb(v)
    assert is_acceptable(w), blurb(w)
    d = category(v, ship, w)
    if not d.startswith('retain'):
      # Show mismatch / ambiguity ... ?
      # m = get_match(z); m.note if m else ...
      yield  (blurb(get_outject(v)) if ship != IREP else "NA",
              blurb(get_outject(w)) if ship != PERI else "NA",
              d,
              witness_comment(AB, v, rel),
              match_comment(AB, v, rel))

# This column is called "category" for consistency with MDD.

def category(v, ship, w):
  rel = get_match(v, None)      # Does the name persist in B?
  mat = rel and (rel.record == w)
  x = get_outject(v)
  y = get_outject(w)
  stay = (get_rank(get_accepted(x), None) ==
          get_rank(get_accepted(y), None))
  if ship == GT:
    if stay: action = ('retain (split)' if mat else 'split')
    else: action = 'add (but it was sort of there?)'
  elif ship == LT:
    if stay: action = ('retain (lump)' if mat else 'lump')
    else: action = 'remove (but sort of staying?)'
  elif ship == EQ:
    # No match - ambiguous match - unique match
    # Different names - same name
    if get_canonical(v) != get_canonical(w):
      action = 'rename'
    else:
      action = 'retain'
  elif ship == DISJOINT:
    action = ('false match' if mat else 'disjoint')
  elif ship == CONFLICT:        # Euler/X calls it 'overlaps'
    action = 'overlap'          # ('reorganize' if mat else ...)
  # Non-RCC5
  elif ship == OVERLAP:
    action = 'intersect in some way'
  elif ship == NOINFO:
    action = 'no information'
  elif ship == SYNONYM:              # a synonym of v <= w
    # x = y1* <= y   i.e. x is being demoted from accepted to synonym (y1)
    t = get_tipe(x, None)
    if t and t == get_tipe(y, None): # ?
      action = 'rename (demotion)'
    else:
      action = 'lump (demotion)'
  elif ship == MYNONYS:              # a synonym of w <= v
    # x >= x1* = y   i.e. x1 is being promoted from synonym to accepted (y)
    t = get_tipe(x, None)
    if t and t == get_tipe(y, None):
      action = 'rename (promotion)'
    else:
      action = 'split (promotion)'
  elif ship == IREP:
    action = 'add'
  elif ship == PERI:
    action = 'remove'
  else:
    action = rcc5_symbol(ship)
  return action

# rel.note will often be the match type, so don't duplicate it

def witness_comment(AB, v, rel):
  w = rel.record
  ship = rel.relationship
  (flush, keep, add) = witnesses(AB, v, w)
  if flush or keep or add and not ship == PERI and not ship == IREP:
    return ("%s|%s|%s" %
            (blurb(flush) if flush else '',
             blurb(keep) if keep else '',
             blurb(add) if add else ''))


def match_comment(AB, v, rel):
  w = rel.record
  comments = []

  # info = match_info(v, None)
  # if info: (us, dir, kind, basis) = info ...
  m1 = get_match(v, None) 
  if m1 and m1.record != w and m1.note and ' ' in m1.note:
    # v and w are matches, so don't need info in both directions
    comments.append(m1.note)

  m2 = get_match(w, None)
  if m2 and m2.record != v and m2.note and ' ' in m2.note:
    # Shouldn't happen, since matches are symmetric
    comments.append(m2.note)

  return '|'.join(comments)

# Report on block(v) - block(w) (synonyms removed)
# and block(w) - block(v) (synonyms added)

def witnesses(AB, v, w):
  b1 = theory.get_block(v, theory.BOTTOM_BLOCK)
  b2 = theory.get_block(w, theory.BOTTOM_BLOCK)
  # Names are all according to B
  flush = choose_witness(AB, b1.difference(b2), w)
  keep  = choose_witness(AB, b1.intersection(b2), w)
  add   = choose_witness(AB, b2.difference(b1), w)
  return (flush, keep, add)

def choose_witness(AB, ids, v):
  have = None
  for id in ids:
    e = theory.exemplar_record(AB, id, v) 
    if not have:
      have = e
    if is_acceptable(e):
      # maybe pick earliest year or something??
      have = e
  return have

# -----------------------------------------------------------------------------
# An alignment is a generator (an Iterable) of species/species articulations.
# An articulation is a pair (record, relation)

def make_alignment(AB):
  return sorted(list(alignment_iter(AB)), key=articulation_order)

# B has priority.  v in A part of AB, w in B part
def articulation_order(art):
  (v, rel) = art
  w = rel.record
  ship = rel.relationship
  if is_species(w) or not is_species(v):
    return (blurb(w), ship)
  else:
    return (blurb(v), ship)

def alignment_iter(AB):
  matches_iter = match_records.match_records(checklist_to_rows(AB.A),
                                             checklist_to_rows(AB.B))
  theory.load_matches(matches_iter, AB)
  theory.theorize(AB)
  span.span(AB)
  #find_acceptables(AB)

  seen = {}

  # Preorder... 
  def traverse(z, infra):
    for (v, w) in taxon_articulators(AB, z, infra): # Normalized
      if is_species(z):
        key = (get_primary_key(v), get_primary_key(w))
        if key in seen or (v == AB.top or w == AB.top):
          pass
        else:
          rel = theory.cross_relation(AB, v, w)
          yield (v, rel)
          seen[key] = True
        infra = z
    for c in get_inferiors(z):
      yield from traverse(c, infra)
  yield from traverse(AB.in_right(AB.B.top), False) # B priority
  yield from traverse(AB.in_left(AB.A.top), False)

# There may be duplicates

def taxon_articulators(AB, z, infra):
  rr = list(theory.opposite_exemplar_records(AB, z))
  if len(rr) > 0:
    for r in rr:        # in same checklist as z
      yield articulator(AB, z, r)
  else:
    # Peripheral
    yield articulator(AB, z, theory.get_reflection(z).record)

def articulator(AB, u, r):
  v = get_acceptable(AB, r) # Accepted and not infraspecific
  if theory.isinA(AB, u):
    return (u, v)
  else:
    return (v, u)

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

def is_acceptable(u):           # in AB
  if False:
    return not see_acceptable(u, None)
  else:
    x = get_outject(u)
    aa = not(unacceptable_locally(x))
    if aa: assert is_accepted(x)
    return aa

# This is inaccurate - doesn't allow nested infraspecific ranks
# Returns reason

def unacceptable_locally(x):
  if not is_accepted(x):
    return "not accepted"
  if get_rank(x, None) == 'species':
    return False
  sup = get_superior(x, None)
  if not sup:
    return False
  if get_rank(sup.record, None) == 'species':
    # Taxa that are subspecies, variety, etc. are not acceptable
    return 'infraspecific'
  return False

# Returns (ancestor, comment)

def get_acceptable(AB, u):
  if False:
    return see_acceptable(y, None) or (u, 'acceptable')
  else:
    x = get_outject(u)
    # Ascend until a species is found
    reason = ''
    while True:
      reason2 = unacceptable_locally(x)
      if not reason2: break
      reason = reason2
      sup = get_superior(x, None)
      if not sup:
        break
      x = sup.record
    assert is_accepted(x), blurb(x)
    # discard reason
    if isinA(AB, u):
      return AB.in_left(x)
    else:
      return AB.in_right(x)

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

def align(A_iter, A_name, B_iter, B_name):
  (al, AB) = ingest(A_iter, A_name, B_iter, B_name)
  return generate_alignment_report(al, AB)

def test():
  def testit(m, n):
    log("\n-- Test: %s + %s --" % (m, n))
    report = align(newick.parse_newick(n), "A",
                   newick.parse_newick(m), "B")
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
        (al, AB) = ingest(csv.reader(a_file), "A",
                          csv.reader(b_file), "B")
        generate_report
        util.write_rows(report, sys.stdout)
