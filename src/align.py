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
from theory import is_accepted, get_accepted, get_tipward

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
  theory.analyze_tipwards(AB)                # also find the 'tipes'
  theory.compute_blocks(AB)                  # sets of 'tipes'
  theory.find_equivalents(AB)
  theory.compute_cross_mrcas(AB)
  span.span(AB)

  seen = {}
  report = []
  report.append(("kind", "A id", "rcc5", "B id",
                 "action", "comment"))

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
          do_row('name', v, w, comment)
        else:
          do_row('topo+name', v, w, comment)
      else:
        # Topology-inspired articulation
        do_row('topo', v, w, comment)

        # Name-inspired articulation
        # Compose a nice comment explaining what's going on.
        # Three cases: v -> w2 -> v, v -> w2 -> v3, v -> w2 -> ---

        m_back = get_match(w2, None)
        if m_back:
          if m_back.record == v:
            if theory.get_equivalent(v, None) == w2:
              comment = "reciprocal name match and equivalent"
            else:
              comment = "reciprocal name match but inequivalent"
          else:
            # this case can't happen, actually
            comment = "one of several matches to target??"
        else:
          comment = "one of several matches to target"
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
    if w:
      action = 'new'
      y = get_outject(w)
    if v and w:
      rel = theory.cross_relation(AB, v, w)
      if rel.relationship == NOINFO:
        return
      action = verbalize(rel, get_rank(v, None) == get_rank(w, None))
      sym = rcc5_symbol(rel.relationship)
    report.append((kind,
                   get_primary_key(x),
                   sym,
                   get_primary_key(y),
                   action,
                   comment))

  for z in checklist.preorder_records(AB):
    assert isinstance(z, prop.Record), blurb(z)
    # Only report on species
    if not is_acceptable(get_outject(z)):
      continue
    w = partner(AB, z)
    (w, comment) = get_acceptable(AB, w) # Acceptable and not a subspecies
    if z == AB.top or w == AB.top:
      pass
    elif theory.isinA(AB, z):
      if theory.get_equivalent(z, None):
        pass
      else:
        do_species(z, w, comment)        # z in A, w in B
    else:
      do_species(w, z, rev_comment(comment))          # w in A, z in B

  return report

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
    w = equ.record
  else:
    xsup = theory.cross_superior(AB, v)
    if not xsup:
      return (None, "top")
    w = xsup.record
  return w

# Result is accepted and not a species ... could be a genus

def get_acceptable(AB, w):
  comment = ""
  y = get_outject(w)
  a = get_accepted(y)
  if a != y:
    y = a
    comment = "via synonym"
  if get_rank(a, None) == 'subspecies':
    sup = get_superior(a, None)
    y = sup.record
    assert is_acceptable(y)
    comment = "via subspecies"
  if theory.isinA(AB, w):
    w = AB.in_left(y) if y else None
  else:
    w = AB.in_right(y) if y else None
  return (w, comment)

def is_acceptable(y):
  if get_accepted(y) != y:
    log("# not accepted: %s" % blurb(y))
    return False
  if (get_rank(y, None) == 'species' or
      get_rank(y, None) == None):
    return True
  else:
    log("# not acceptable: %s" % blurb(y))
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
