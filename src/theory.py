#!/usr/bin/env python3

import types
import property as prop, checklist, workspace

from util import log
from checklist import *
from workspace import *
from match_records import match_records

# So how does this work.
# We determine an alignment heuristically.
# Concurrently, given a developing alignment, we determine a spanning
# tree.
# The alignment is represented as...
#   - Relateds values of the 'equivalent' property
#   - Relateds values of the 'superior' property
#   - A set of dropped nodes (conflicting or garbage)
# umm, that's pretty much it??

def analyze(AB):
  m = match_records(checklist_to_rows(AB.A), checklist_to_rows(AB.B))
  load_matches(m, AB)
  mtrm(AB)

# -----------------------------------------------------------------------------
# Find tipward matches

(get_equivalent, set_equivalent) = prop.get_set(prop.get_property("TRM"))
(get_trm, set_trm) = prop.get_set(prop.get_property("TRM"))

def mtrm(AB):
  trm(AB.A, AB.in_left, set_trm)
  def set_if_mutual(z, m):
    set_trm(z, m)           # Nonmutual, might be of interest

    z2 = m.other
    m2 = get_trm(z2, None)      # tipward match to/in A
    if m2 and m2.other == z:
      if monitor(z): log("# MTRM: %s :=: %s" % (blurb(z), blurb(z2)))
      set_equivalent(z, m)
      set_equivalent(z2, m2)

  trm(AB.B, AB.in_right, set_if_mutual)    # of AB.flip()

def trm(A, in_left, set_trm):
  ensure_inferiors_indexed(A)
  empty = [] # want ()
  def traverse(x):
    seen = False
    # Work around bug in property.py
    for c in get_children(x, None) or empty:
      seen = traverse(c) or seen
    for c in get_synonyms(x, None) or empty:
      seen = traverse(c) or seen
    if not seen and not is_top(x):
      z = in_left(x)
      m = get_match(z)
      if m and m.relation == EQ:
        if monitor(z): log("# TRM: %s = %s" % (blurb(z), blurb(m.other)))
        set_trm(z, m)
        seen = z
    return seen
  traverse(A.top)

def loser():
  if False: yield "lose"

# -----------------------------------------------------------------------------
# Load/dump a set of provisional matches (could be either record match
# or taxonomic matches... but basically, record matches)

# Fields of match records <A (matched), relation, B (taxon), remark>
get_match_relation = prop.getter(prop.get_property("relation"))
get_matched_key = prop.getter(prop.get_property("matched_id"))
get_match_note = prop.getter(prop.get_property("match_note"))

def load_matches(row_iterator, AB):

  header = next(row_iterator)
  plan = prop.make_plan_from_header(header)
  match_count = 0
  miss_count = 0
  for row in row_iterator:
    # row = [matchID, rel, taxonID, remark]
    match = prop.construct(plan, row)
    x = y = None
    xkey = get_matched_key(match, None)
    if xkey:
      x_in_A = look_up_record(AB.A, xkey)
      if x_in_A:
        x = AB.in_left(x_in_A)
    ykey = get_primary_key(match, None)
    if ykey:
      y_in_B = look_up_record(AB.B, ykey)
      if y_in_B:
        y = AB.in_right(y_in_B) 

    # The columns of the csv file
    rel = rcc5_relation(get_match_relation(match))
    note = get_match_note(match, MISSING)

    # x or y might be None with rel=NOINFO ... hope this is OK

    if x:
      set_match(x, Related(rel, y, "record match", note))
    if y:
      set_match(y, Related(reverse_relation(rel), x, "record match",
                           reverse_note(note)))
    if x and y: match_count += 1
    else: miss_count += 1

  log("# %s matches, %s misses" % (match_count, miss_count))

def reverse_note(note):
  return "↔ " + note            # tbd: deal with 'coambiguous'

def record_match(x):
  ship = get_match(x)
  if ship.relation == EQ: return ship.other
  return None

"""
TBD: filter out seniors
  seniors = 0
      # Filter out senior synonyms here
      if is_senior(item, accepted_item):
        seniors += 1
      else:
  if seniors > 0:     # Maybe interesting
    print("-- Suppressed %s senior synonyms" % seniors,
          file=sys.stderr)

"""

# -----------------------------------------------------------------------------

if __name__ == '__main__':
  def testit(m, n):
    AB = workspace.workspace_from_newicks(m, n)
    analyze(AB)
    workspace.show_workspace(AB)
  testit("a", "a")
  testit("(c,d)a", "(c,e)b")
