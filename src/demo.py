#!/usr/bin/env python3

# Report generation for alignment procedure

import sys, csv, argparse
import util, workspace, align, linkage, theory, exemplar
import rows

from checklist import *
from workspace import *
from eulerx import generate_eulerx
from theory import get_equivalent
from align import is_species, is_acceptable

# Returns a row generator

def generate_alignment_report(al, AB):
  yield ("A id", "A name", "rcc5", "B name", "B id",
         "category", "note", "witnesses", "match")
  seen = set()
  for (u, rel) in al:
    ship = rel.relationship
    v = rel.record
    x = theory.get_outject(u)
    y = theory.get_outject(v)
    d = category(u, ship, v)

    # What is this about ???
    if (d.startswith('drop') and get_match(u, None) != None and
        not get_primary_key(u) in seen):
      if monitor(u):
        log("# Removing %s, why?" % blurb(u))
        diagnose_match(u)         # logs
      seen.add(get_primary_key(u))

    yield  (get_primary_key(x),
            blurb(x),
            rcc5_symbol(ship),
            blurb(y),
            get_primary_key(y),
            d, rel.note,
            witness_comment(AB, u, rel),
            match_comment(AB, u, rel))

def generate_short_report(al, AB):
  yield ("A name", "rcc5", "B name", "category", "witnesses", "match")
  for (u, rel) in al:
    w = rel.record

    # Elide equalities lacking name change ...
    ship = rel.relationship
    if ship == EQ and get_canonical(u) == get_canonical(w):
      continue

    if ((is_species(u) or is_species(w)) and
        is_acceptable(u) and is_acceptable(w)):  # filter out subspecies too
      if plausible(AB, u, ship, w):
        continue
      d = category(u, ship, w)
      if (not d.startswith('retain')) or get_canonical(u) != get_canonical(w):
        # Show mismatch / ambiguity ... ?
        # m = get_match(z); m.note if m else ...
        yield  (blurb(get_outject(u)) ,# if ship != IREP else MISSING,
                rcc5_symbol(ship),
                blurb(get_outject(w)) ,# if ship != PERI else MISSING,
                d,
                witness_comment(AB, u, rel),
                match_comment(AB, u, rel))

# Is relationship mmerely copied over from a source checklist?

def plausible(AB, u, ship, v):
  if (ship == DISJOINT):
    sup1 = get_superior(u, None)
    sup2 = get_superior(v, None)
    if sup1 and sup2:
      eq1 = get_equivalent(AB, sup1.record)
      eq2 = get_equivalent(AB, sup2.record)
      r1 = eq1.record if eq1 else None
      r2 = eq2.record if eq2 else None
      if (r1 == r2):
        #log("# unsurprising! %s %s" % (blurb(u), blurb(v)))
        return True
  return False

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
    else: action = 'add'
  elif ship == LT:
    if stay: action = ('retain (lump)' if mat else 'lump')
    else: action = 'drop'
  elif ship == EQ:
    # No match - ambiguous match - unique match
    # Different names - same name
    if get_canonical(v) != get_canonical(w):
      action = 'rename'
    else:
      action = 'retain'
  elif ship == DISJOINT:
    action = ('false match' if mat else 'disjoint')
  elif ship == OVERLAP:        # Euler/X calls it 'overlaps'
    action = 'overlap'          # ('reorganize' if mat else ...)
  # Non-RCC5
  elif ship == INTERSECT:
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
    action = 'drop'
  else:
    action = rcc5_symbol(ship) + ' ??'
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


# Ideally, display only surprises: where name match and hierarchy 
# match give different results.  That would only be the = case.
# And = is suggested by any positive matches where either match set
# contains the opposite node.

def match_comment(AB, v, rel):
  w = rel.record
  ws = get_matches(v)
  vs = get_matches(w)
  x1 = v in vs
  x2 = w in ws
  if (x1 or x2) == (rel.relationship == EQ):
    return MISSING
  else:            # Surprising
    return ("a -> %s %s b -> %s %s a" %
            (len(ws), ';' if x1 else '/', len(vs), ';' if x2 else '/'))

def match_comment1(AB, v, rel):
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
  b1 = theory.get_block(v)
  b2 = theory.get_block(w)
  if len(b1) == 0 or len(b2) == 0:
    return None
  # Names are all according to B
  flush = choose_witness(AB, b1.difference(b2), w)
  keep  = choose_witness(AB, b1.intersection(b2), w)
  add   = choose_witness(AB, b2.difference(b1), w)
  if len(b1) + len(b1) > 10000:
    debug_withnesses(AB, v, w)
  return (flush, keep, add)

def debug_withnesses(AB, v, w):
  b1 = theory.get_block(v)
  b2 = theory.get_block(w)
  def foo(v, which, w, s):
    if False:
      log("# choose_witness: %s %s %s has %s exemplars" %
          (blurb(v), which, blurb(w), len(s)))
    pass
  foo(v, "-", w, b1.difference(b2))
  foo(v, "∩", w, b1.intersection(b2))
  foo(w, "-", v, b2.difference(b1))

def choose_witness(AB, ids, v):
  have = None
  count = 0
  for id in ids:
    e = exemplar.xid_to_record(AB, id, v) 
    if not have:
      have = e
    if is_acceptable(e):     # ?
      # maybe pick earliest year or something??
      have = e
      break
    if count == 100:
      have = e
      break
    count += 1
  return have


def test():
  def testit(m, n):
    log("\n-- Test: %s + %s --" % (m, n))
    (report, eulerx, tipwards) = demo(newick.parse_newick(n), "A",
                                      newick.parse_newick(m), "B")
    util.write_rows(report, sys.stdout)
    for line in eulerx:
      print(line, file=sys.stdout)
  testit("a", "a")              # A + B
  #testit("(c,d)a", "(c,e)b")
  #testit("((a,b)e,c)d", "(a,(b,c)f)D")
  #testit("(a,b*)c", "(a)c")
  testit("((b*)a)c", "(a,b)c")     # WRONG WRONG

def demo(A_iter, A_name, B_iter, B_name):
  AB = ingest_workspace(A_iter, A_name, B_iter, B_name)
  al = align.generate_alignment(AB)
  return (generate_alignment_report(al, AB),
          generate_short_report(al, AB),
          generate_eulerx(AB, al))

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
  parser.add_argument('--eulerx', help="where to put the Euler/X version of the alignment",
                      default=None)
  parser.add_argument('--short', help="where to put the differential alignment report",
                      default=None)
  parser.add_argument('--long', help="where to put the full alignment report",
                      default=None)
  parser.add_argument('--test', action='store_true', help="run smoke tests")
  args=parser.parse_args()
  if args.test:
    test()
  else:
    a_name = args.Aname
    b_name = args.Bname
    a_path = args.A
    b_path = args.B
    e_path = args.eulerx
    d_path = args.short         # diff
    l_path = args.long          # long alignment report
    assert a_path != b_path
    with rows.open(a_path) as a_rows:
      with rows.open(b_path) as b_rows:
        AB = ingest_workspace(a_rows.rows(), b_rows.rows(),
                              A_name=a_name, B_name=b_name)
        linkage.find_links(AB)
        al = list(align.generate_alignment(AB))
    if l_path or (not d_path and not e_path):
      with rows.open(l_path, "w") as l_gen:
        l_gen.write_rows(generate_alignment_report(al, AB))
    if d_path:
      with rows.open(d_path, "w") as d_gen:
        d_gen.write_rows(generate_short_report(al, AB))
    if e_path:
      with util.stdopen(e_path, "w") as e_gen:
        for line in generate_eulerx(AB, al):
          print(line, file=e_gen)
