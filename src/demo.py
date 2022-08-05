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
from theory import is_accepted, get_accepted

def demo(A_iter, B_iter):
  A = rows_to_checklist(A_iter, {"name": "A"})  # meta
  B = rows_to_checklist(B_iter, {"name": "B"})  # meta
  AB = workspace.make_workspace(A, B, {"name": "AB"})
  matches_iter = match_records.match_records(checklist_to_rows(A),
                                             checklist_to_rows(B))
  theory.load_matches(matches_iter, AB)
  theory.analyze_tipwards(AB)                # also find the 'tipes'
  theory.compute_blocks(AB)                  # sets of 'tipes'
  theory.find_equivalents(AB)
  theory.compute_cross_mrcas(AB)
  span.span(AB)
  return generate_report(AB)

def generate_distinct(AB):
  def generate(z):
    children = tuple(get_inferiors(z))
    if is_acceptable(get_outject(z)):
      rel = get_match(z, None)
      w = rel.record if rel else None
      eq = w and is_acceptable(get_outject(w))
      if eq and theory.isinB(AB, z):
        return     # *** Has been processed differently
      if z != AB.top:
        yield z
      if eq and theory.isinA(AB, z):
        children = children + tuple(get_inferiors(w))
    for c in sorted(children, key=sort_key):
      yield from generate(c)
  yield from generate(AB.top)

# Returns an Iterable

def generate_report(AB):
  # Too difficult to do a generator
  report = []

  report.append(("A id", "A name", "rcc5", "B name", "B id",
                 "action", "comment"))

  def do_row(v, w, comment):
    x_info = ('', '')
    y_info = ('', '')
    sym = ''
    if v:
      action = 'deleted'
      x = get_outject(v)
      x_info = (get_primary_key(x),
                get_canonical(x, None) or get_scientific(x))
    if w:
      action = 'new'
      y = get_outject(w)
      y_info = (get_primary_key(y),
                get_canonical(y, None) or get_scientific(y))
    if v and w:
      rel = theory.cross_relation(AB, v, w)
      if rel.relationship == NOINFO:
        return
      action = verbalize(rel)
      sym = rcc5_symbol(rel.relationship)
    report.append((x_info[0],
                   x_info[1],
                   sym,
                   y_info[1],
                   y_info[0],
                   action,
                   comment))

  def promote(z):
    a = get_accepted(z)
    if a:
      c = "via synonym")
      if get_rank(a.record, None) == 'subspecies':
        a = get_superior(a.record).record
        c = "via subspecies"
    else:
      a = z
      c = ""
    if get_rank(a.record, None) == 'species':
      return (a, c)
    else:
      return None

  for z in preorder_rows(AB):
    # Skip if z is not a species
    if get_rank(z, None) == 'species':
      if isinB(AB, z):
        equ = get_equivalent(z, None)
        if equ:
          do_row(equ.record, z, "foo")
        else:
          rel = theory.cross_superior(AB, z)
          (w, comment) = promote(rel.record)
          do_row(w, z, comment)
      elif not get_equivalent(z, None):
        rel = theory.cross_superior(AB, z)
        (v, comment) = promote(rel.record)
        do_row(z, v, comment)


  # Eight cases, of which 3 are to be ignored.

  #   yy    y     x     xx
  # 1       acc = acc
  # 2       acc = syn < acc
  # 4 acc > syn = acc
  # 5 acc > syn = syn < acc   IGNORE

  # 3       acc.      < acc
  # 6 acc > syn.      < acc   IGNORE

  # 7 acc >       acc
  # 8 acc >       syn.        IGNORE
  # 9     >                   CAN'T HAPPEN

  for z in generate_distinct(AB):
    # N.b. z is acceptable...

    def when_A(x):
      # cases 7-8.  x is acceptable
      # case 7, acc >       acc.
      if monitor(z): log("when_A")
      rel = get_match(z, None) # in B
      if rel:
        if rel.record:
          if rel.relationship == EQ:
            pass            # handled under when_B case
          else:             # NOINFO
            do_row(z, rel.record, 'multiple A records match B name')
        else:
          (p, comment) = partner(AB.swap(), z)
          # 'A name has multiple matches in B'
          do_row(z, p, comment)

    def when_B(y):
      # cases 1-6.  y is acceptable
      rel = theory.get_equivalent(z, None)
      if rel:
        #       acc = acc
        #       acc = syn < acc
        # acc > syn = acc
        v = rel.record
        x = get_outject(rel.record)
        if is_acceptable(x):
          # case 1, acc = acc
          if monitor(z): log("when_B case 1")
          do_row(v, z, '')
        else:
          # case 2
          if monitor(z): log("when_B case 2")
          do_row(AB.in_left(get_acceptable(x)), z,
                 'via A synonym')

      else:
        # case 3, acc. ... < acc
        if monitor(z): log("when_B case 3")
        rel = get_match(z, None)
        if rel:
          if rel.record:
            if rel.relationship == EQ:
              do_row(rel.record, z, 'record match')
            else:
              do_row(rel.record, z, 'B name has multiple matches in A')
          else:
            rel = theory.cross_superior(AB, z)
            do_row(rel.record, z,
                   'B name has multiple matches in A')
        else:
          (p, comment) = partner(AB, z)
          do_row(p, z, comment)

    AB.case(z, when_A, when_B)

  return report

def partner(AB, v):
  x = get_outject(v)
  a = get_accepted(x)
  if a != x:
    comment = "via synonym"
    v = AB.in_A(a)
  rel = theory.get_equivalent(v, None)
  w = rel.record if rel else None
  if w:
    comment = "equivalent"
  else:
    if False:
      rel = get_match(v, None)
      w = rel.record if rel else None
      if w:
        comment = "name match"
    rel = theory.cross_superior(AB, v)
    w = rel.record
    if w:
      comment = "nearest enclosing"

  if x:
    return get_acceptable(x)       # .... hmm ....
  else:
    if x:
      return get_acceptable(x)
    else:
      rel = theory.cross_superior(AB, v)
      return (rel.record, '(something) only')

def get_acceptable(x):
  a = AB.in_B(x)

  return (a, 'comment')

def is_acceptable(x):
  return is_accepted(x) and get_rank(x, 'species') == 'species'

def sort_key(c):
  return (get_canonical(c, ''),
          get_scientific(c, ''),
          get_primary_key(c, ''))

def verbalize(rel):
  if rel.relationship == GT:
    action = 'split'
  elif rel.relationship == EQ:
    action = ''
  elif rel.relationship == LT:
    action = 'lumped'
  elif rel.relationship == DISJOINT:
    action = 'moved'
  elif rel.relationship == CONFLICT:
    action = 'reorganized'
  elif rel.relationship == OVERLAP:
    action = 'overlap'
  else:
    action = rcc5_symbol(rel.relationship)
  return action

# -----------------------------------------------------------------------------

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    """)
  parser.add_argument('--A', help="the A checklist, as path name or -",
                      default='-')
  parser.add_argument('--B', help="the B checklist, as path name or -",
                      default='-')
  args=parser.parse_args()
  a_path = args.A
  b_path = args.B
  assert a_path != b_path
  with util.stdopen(a_path) as a_file:
    with util.stdopen(b_path) as b_file:
      report = demo(csv.reader(a_file),
                    csv.reader(b_file))
      util.write_rows(report, sys.stdout)
