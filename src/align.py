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

    seen = {}
    if (d.startswith('remove') and get_match(v, None) != None and
        not get_primary_key(v) in seen):
      log("# Removing %s, why?" % blurb(v))
      diagnose_match(v)         # logs
      seen[get_primary_key(v)] = True

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
      yield  (blurb(get_outject(v)) if ship != IREP else MISSING,
              blurb(get_outject(w)) if ship != PERI else MISSING,
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
  probe = witnesses(AB, v, w)
  if not probe: return ''
  (flush, keep, add) = probe
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
  if len(b1) == 0 or len(b2) == 0:
    return None
  # Names are all according to B
  flush = choose_witness(AB, b1.difference(b2), w)
  keep  = choose_witness(AB, b1.intersection(b2), w)
  add   = choose_witness(AB, b2.difference(b1), w)
  if len(b1) + len(b1) > 1000:
    debug_withnesses(AB, v, w)
  return (flush, keep, add)

def debug_withnesses(AB, v, w):
  b1 = theory.get_block(v, theory.BOTTOM_BLOCK)
  b2 = theory.get_block(w, theory.BOTTOM_BLOCK)
  def foo(v, which, w, s):
    print("# choose_witness: %s %s %s has %s exemplars" %
          (blurb(v), which, blurb(w), len(s)),
          file=sys.stderr)
  foo(v, "-", w, b1.difference(b2))
  foo(v, "∩", w, b1.intersection(b2))
  foo(w, "-", v, b2.difference(b1))

def choose_witness(AB, ids, v):
  have = None
  for id in ids:
    e = theory.exemplar_record(AB, id, v) 
    if not have:
      have = e
    if is_acceptable(e):     # ?
      # maybe pick earliest year or something??
      have = e
      break
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
    vs = get_matches(z)
    for v in vs:
      yield articulator(AB, z, v)
    # Peripheral - has to attach somewhere
    re = theory.get_reflection(z).record
    yield articulator(AB, z, re)

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
        report = generate_alignment_report(al, AB)
        util.write_rows(report, sys.stdout)
