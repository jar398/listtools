#!/usr/bin/env python3

import sys, csv
from typing import NamedTuple, Any

import property as prop
import checklist

from util import log, MISSING
from checklist import *
from coproduct import *

# -----------------------------------------------------------------------------
# Sum / coproduct / merged checklist / theory workspace
# Could be either just the matches (not a checklist), or the matches
# plus overrides (synthetic checklist)

# Returns B side

REUSE_KEYS = True

def make_workspace(A, B, meta=None):
  Q = prop.make_context()       # allows overriding A and/or B
  (get_inject, set_inject) = prop.get_set(inject_prop, context=Q)

  register = prop.get_registrar(primary_key_prop, Q)
  pk_counter = [0]
  missing = False

  def ensure_injected(x):
    z = get_inject(x, None)
    if not z:
      z = prop.clone(x)
      set_inject(x, z)          # contextual
      set_outject(z, x)
      set_source(z, AB)
      have_key = get_primary_key(z)
      if not have_key or not REUSE_KEYS:
        pk_counter[0] += 1
        set_primary_key(z, str(pk_counter[0]))
      elif not ('!' in have_key):
        key = "%s!%s" % (get_source_name(x), have_key)
        set_primary_key(z, key)
      register(z)
    return z

  def _in_left(x):
    assert get_source(x) == A
    return ensure_injected(x)
  def _in_right(y):
    assert get_source(y) == B
    return ensure_injected(y)
  def _case(z, when_left, when_right):
    w = get_outject(z, None)
    assert w
    if get_source(w) == A:
      return when_left(w)
    else:
      assert get_source(w) == B
      return when_right(w)
  AB = Coproduct(_in_left, _in_right, _case)
  AB.context = Q

  atop = AB.in_left(A.top)
  btop = AB.in_right(B.top)

  set_equated(atop, Related(EQ, btop, "top"))
  set_equated(btop, Related(EQ, atop, "top"))
  AB.top = btop

  AB.topship = Related(ACCEPTED, AB.top, "uninitialized")
  AB.indexed = False
  AB.meta = meta

  AB.A = A           # need ??
  AB.B = B

  # Force local copies of all source records
  for y in all_records(B): AB.in_right(y) # not including top
  for x in all_records(A): AB.in_left(x)  # not including top
  #log("# taxonID counter: %s" % pk_counter[0])

  return AB

# Is given synonym usage a senior synonym of its accepted usage?
# In the case of splitting, we expect the synonym to be a senior
# synonym of the item.

# We could also look to see if taxonomicStatus is 'senior synonym'.

# Hmm, if there's exactly one senior synonym, we should remember it as
# being the 'split from' taxon.

def is_senior(synonym_item, accepted_item):
  syn_year = get_year(synonym_item, None)  # the synonym
  acc_year = get_year(accepted_item, None)
  if syn_year and acc_year:
    # Older name (senior) has priority over newer (junior)
    # but if junior is accepted we don't care about the senior.
    if syn_year <= acc_year:
      # synonym is older than accepted, so syn > acc.  Shouldn't
      # happen.  (these are generated by MDD)
      print("# Flushing senior synonym '%s' of '%s'" %
            (get_scientific(synonym_item),
             get_scientific(accepted_item)))
      return True
    else:
      # synonym is newer than accepted, so syn < acc.  Junior.
      #  (split)
      return False
  else:
    return False

# -----------------------------------------------------------------------------

# Output with additional columns needed by report.py

def workspace_to_rows(AB, props=usual_props):

  # Do we need these at all?

  def get_id_a(z, default=None):
    x = left_persona(AB, z)
    return get_primary_key(x) if x else default

  def get_id_b(z, default=None):
    y = right_persona(AB, z)
    return get_primary_key(y) if y else default

  workspace_props = props + \
      (prop.get_property("taxonID_A", getter=get_id_a),
       prop.get_property("taxonID_B", getter=get_id_b),
       prop.get_property("match_key", getter=recover_match_key),
       prop.get_property("match_note", getter=recover_match_note))

  return checklist_to_rows(AB, props)

def recover_match_key(z, default=None):
  m = get_match(z, None)
  return get_primary_key(m.record) if m else default

def recover_match_note(z, default=None):
  m = get_match(z, None)
  return m.note if m else default

# -----------------------------------------------------------------------------

import newick

def workspace_from_newicks(m, n):
  B = rows_to_checklist(newick.parse_newick(n),
                        {"name": "B"})  # meta
  A = rows_to_checklist(newick.parse_newick(m),
                        {"name": "A"})  # meta
  return make_workspace(A, B, {"name": "AB"})

def show_workspace(AB):
  writer = csv.writer(sys.stdout)
  rows = list(preorder(AB))
  for row in rows:
    writer.writerow(row)
  print(' = ' + newick.compose_newick(rows))

if __name__ == '__main__':
  def testit(m, n):
    AB = workspace_from_newicks(m, n)
    show_workspace(AB)
  testit("(c,d)a", "(c,e)b")
