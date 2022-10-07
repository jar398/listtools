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
from theory import is_accepted, get_accepted, get_tipward, local_sup

def align(A_iter, B_iter):
  A = rows_to_checklist(A_iter, {"name": "A"})  # meta
  B = rows_to_checklist(B_iter, {"name": "B"})  # meta
  return generate_alignment(A, B)

# Returns an Iterable

def generate_alignment(A, B):
  AB = workspace.make_workspace(A, B, {"name": "AB"})
  matches_iter = match_records.match_records(checklist_to_rows(A),
                                             checklist_to_rows(B))
  theory.load_matches(matches_iter, AB)
  theory.theorize(AB)
  span.span(AB)

  seen = {}
  report = []
  report.append(("kind", "A id", "rcc5", "B id",
                 "action", "note", "comment"))

  def do_species(v, w, comment):
    # Three articulations of interest:
    #    w (the cross_relation)
    #    m (the name match)
    # Show the *distinct* articulations among these.

    m = get_match(v, None)      # name match
    if m and m.record:
      w2 = m.record
      if w2 == w:
        # Topology- and name-inspired articulation
        # If block is empty then name only, not topo+name
        if get_tipward(v, None):
          do_row('tipward', v, w, comment)
        else:
          do_row('match', v, w, comment)
      else:
        # Topology-inspired articulation
        do_row('topo', v, w, comment)

        # Name-inspired articulation
        # Compose a nice comment explaining what's going on.
        # Three cases: v -> w2 -> v, v -> w2 -> v3, v -> w2 -> ---

        m_back = get_match(w2, None) # reciprocal name match
        tag = None
        if m_back:
          if m_back.relationship != NOINFO:
            if m_back.record == v:
              # v -> w2 -> v
              tag = "reciprocal name match"
              pass
            else:
              # v -> w2 -> v'
              # Should not happen; EQ implies it's reciprocal.
              tag = "nonreciprocal name match (%s)" % (blurb(m_back.record),)
          else:
            if m_back.record:
              # v -> w2 -> v'
              tag = "matches record that matches ambiguously (1)"
            else:
              # Doesn't happen?
              tag = "matches record that matches ambiguously (2)"
        else:
          # v -> w2 -> nothing
          # Doesn't happen?
          tag = "matches record that matches ambiguously (3)"
        if comment:
          if tag: comment = "%s; %s" % (comment, tag)
        else:
          comment = tag
        do_row('name', v, m.record, comment)
    else:
      # Topology-inspired articulation only
      do_row('topo', v, w, comment)

  def do_row(kind, v, w, comment):
    key = (get_primary_key(v), get_primary_key(w))
    if key in seen: return
    seen[key] = True

    x_info = ('', '')
    y_info = ('', '')
    sym = ''
    if v:
      action = 'deleted'
      x = get_outject(v)
      note = '-'
    if w:
      action = 'new'
      y = get_outject(w)
      note = '-'
    if v and w:
      rel = theory.cross_relation(AB, v, w)
      if rel.relationship == NOINFO:
        return
      action = verbalize(rel, get_rank(v, None) == get_rank(w, None))
      sym = rcc5_symbol(rel.relationship)
      note = rel.note
    assert note
    report.append((kind,
                   get_primary_key(x),
                   sym,
                   get_primary_key(y),
                   action,
                   note,
                   comment))

  for z in checklist.preorder_records(AB):
    assert isinstance(z, prop.Record), blurb(z)
    # Only report on species
    if not is_acceptable(get_outject(z)):
      continue
    rel = partner(AB, z)
    w = rel.record
    (w, comment) = get_acceptable(AB, w) # Acceptable and not a subspecies
    if z == AB.top or w == AB.top:
      pass
    elif theory.isinB(AB, z):
      do_species(w, z, rev_comment(comment)) # w in A, z in B
    elif theory.get_equivalent(z, None):
      # Redundant, generated on previous pass
      pass
    else:
      # In A but not in B
      do_species(z, w, comment)        # z in A, w in B
  return (report, specimens_table(AB))

# Requires spanning tree

def specimens_table(AB):
  yield("checklist", "species_id", "species_canonical", "specimen_id", "specimen_canonical")

  def traverse(z, species):
    if not theory.get_equivalent(z, None) or theory.isinA(AB, z):
      if get_rank(z, None) == 'species': species = z
      specimen_id = theory.get_specimen_id(AB, z)
      if specimen_id:
        # Specimen in some species.  Prefer info from B
        (t0, t) = AB.specimen_taxa[specimen_id]

        # The species that the specimen belongs to
        tag = "A" if theory.isinA(AB, z) else "B"
        x = get_outject(species) if species else None

        yield(tag,
              get_primary_key(x) if x else None,
              blurb(x) if x else None,
              specimen_id,
              blurb(t))
    for c in get_inferiors(z):
      yield from traverse(c, species)
  yield from traverse(AB.top, None)

def get_species(x):
  sup = get_superior(x, None)
  if sup == None:
    return None
  elif sup.relationship == ACCEPTED and get_rank(x, None) == 'species':
    return x
  else:
    return get_species(sup.record)

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

# Result is accepted and not a species ... could be a genus

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
    b = get_accepted(sup.record)
    if b and b != a:
      # 'Canis indica' has rank subspecies and is in 'Canis'
      log("# %s subspecies of %s synonym of %s" % (blurb(w), blurb(a), blurb(b)))
      comment = "via subspecies of synonym"
      a = b
  if theory.isinA(AB, w):
    w = AB.in_left(a)
  else:
    w = AB.in_right(a)
  return (w, comment)

def is_acceptable(y):
  if get_accepted(y) != y:
    return False
  if (get_rank(y, None) == 'species' or
      get_rank(y, None) == None):
    return True
  else:
    return False

def sort_key(c):
  return (get_canonical(c, ''),
          get_scientific(c, ''),
          get_primary_key(c, ''))

def verbalize(rel, same_rank):
  if rel.relationship == GT:
    action = 'split' if same_rank else 'added'
  elif rel.relationship == EQ:
    action = ''
  elif rel.relationship == LT:
    action = 'lump' if same_rank else 'removed'
  elif rel.relationship == DISJOINT:
    action = 'move'
  elif rel.relationship == CONFLICT:
    action = 'reorganize'
  elif rel.relationship == OVERLAP:
    action = 'overlap'
  else:
    action = rcc5_symbol(rel.relationship)
  return action

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
