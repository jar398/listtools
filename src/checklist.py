#!/usr/bin/env python3

# Source and sum checklists

import sys, csv
import property as prop

from coproduct import *
from util import log

checklist_prop = prop.get_property("checklist")    # which checklist does this belong to?
(get_source_checklist, set_source_checklist) = prop.get_set(checklist_prop)

# ---- Common code for both kinds of checklist

# This approach makes checklists rather heavyweight... that's OK for now

primary_key_prop = prop.get_property("taxonID")
canonical_prop = prop.get_property("canonicalName")
scientific_prop = prop.get_property("scientificName")
rank_prop = prop.get_property("taxonRank")
year_prop = prop.get_property("year")  # http://rs.tdwg.org/dwc/terms/year
parent_key_prop = prop.get_property("parentNameUsageID")
accepted_key_prop = prop.get_property("acceptedNameUsageID")
  
(get_source_checklist, set_source_checklist) = prop.get_set(checklist_prop)
(get_parent_key, set_parent_key) = prop.get_set(parent_key_prop)
(get_accepted_key, set_accepted_key) = prop.get_set(accepted_key_prop)
(get_year, set_year) = prop.get_set(year_prop)

(get_canonical, set_canonical) = prop.get_set(canonical_prop)

def init_checklist(C):
  C.columns = {}                # Needed by property.py
  C.index = Index()

  # Properties that are usually CSV columns.  Ambient but can be
  # overridden contextually.
  (C.get_primary_key, C.set_primary_key) = \
    prop.get_set(primary_key_prop, context=C)
  (C.get_canonical, C.set_canonical) = \
    prop.get_set(canonical_prop, context=C)
  (C.get_scientific, C.set_scientific) = \
    prop.get_set(scientific_prop, context=C)
  (C.get_rank, C.set_rank) = \
    prop.get_set(rank_prop, context=C)
  # Not CSV columns.  - tbd, use this ambiently as well?
  (C.get_checklist, C.set_checklist) = \
    prop.get_set(checklist_prop, context=C)       # ?????
  (C.get_same, C.set_same) = \
    prop.get_set(prop.get_property("same"), context=C)
  (C.get_parent, C.set_parent) = \
    prop.get_set(prop.get_property("parent"), context=C)
  (C.get_accepted, C.set_accepted) = \
    prop.get_set(prop.get_property("accepted"), context=C)
  (C.get_children, C.set_children) = \
    prop.get_set(prop.get_property("children"), context=C)
  (C.get_synonyms, C.set_synonyms) = \
    prop.get_set(prop.get_property("synonyms"), context=C)

  # Convenience methods
  def get_superior(x):
    return C.get_parent(x, None) or C.get_accepted(x, None)
  def get_inferiors(p):
    return C.get_children(p, []) + C.get_synonyms(p, [])
  C.get_superior = get_superior
  C.get_inferiors = get_inferiors

class Index:
  def __init__(self):
    self.by_primary_key_dict = {}     # from DwC taxonID column.  maybe others
    self.by_canonical_dict = {}     # from DwC taxonID column.  maybe others

# -----------------------------------------------------------------------------
# Sum / coproduct / merged checklist

def contains(AB, x):
  return AB.coproduct.contains(x)

# Returns B side

def make_sum(A, B, kind, rm_sum=None):
  def _join(x, y):
    if x and y:
      co.set_same(x, y)         # x is junior synonym of y
      co.set_same(y, x)         # y is senior synonym of x
    return y or x               # prefer B to A
  def _split(z):
    w = co.get_same(z, None)
    if B.contains(z):
      return (w, z)       # z is senior synonym if w isn't None
    else:
      assert A.contains(z)
      return (z, w)       # z is junior synonym if w isn't None

  def _contains(x):
    return A.contains(x) or B.contains(x)

  """
  # Another way to do it, avoiding contextual properties
  A_injections = {}
  B_injections = {}
  def _join(x, y):
    if x: assert A.contains(x)
    if y: assert in_checklist(y, B)
    z = (y and B_injections.get(y, None)) or (x and A_injections.get(x, None))
    if not z:
      z = make_injected(x, y)
      if x: prop.mep_set(A_injections, x.id, z)
      if y: prop.mep_set(B_injections, y.id, z)
    return z
  def _split(z):
    return (out_A(z, None), out_B(z, None))
  """
  co = Coproduct(_join, _split)
  co.A = A
  co.B = B
  co.contains = _contains    # accessed via Side
  co.rm_sum = rm_sum
  co.kind = kind
  init_checklist(co)

  (co.AB.A, co.AB.B) = (A, B)
  (co.BA.A, co.BA.B) = (B, A)

  return co.AB

"""
  def record_match(self, x):
    # This doesn't feel right
    y = self.out_left(u)             # y in A
    r = self.out_right(self.rm_sum.in_left(y))
    return self.rm_sum.in_right(r) if r else None
"""

# -----------------------------------------------------------------------------
# Source checklists

class Source:
  def __init__(self, name):
    init_checklist(self)
    self.name = name
  def contains(self, x):
    return get_source_checklist(x) == self

# Read and write Darwin Core files
#   Stream of row <-> Source structure

# Darwin Core CSV format hierarchy file ingest

def load_source(iterator, name):
  S = Source(name)
  for usage in load_usages(S, iterator):
    pass
  link_to_superiors(S)   # sets parents, accepteds, and S.roots
  collect_inferiors(S)
  return S

# Two ways to make these things.  One is using plans (with `construct`).
# The other is by using the custom constructor (`make_usage`, two arguments).

make_usage = prop.constructor(primary_key_prop, checklist_prop)

def load_usages(S, iterator):
  header = next(iterator)
  plan = prop.make_plan_from_header(header)
  fixup = []
  for row in iterator:
    usage = prop.construct(plan, row)
    set_source_checklist(usage, S)
    assert get_source_checklist(usage)
    assert S.get_checklist(usage)
    # S.set_checklist(usage, S)
    key = S.get_primary_key(usage)
    S.index.by_primary_key_dict[key] = usage # Register
    yield usage

def lookup_usage(S, key, comes_from):
  probe = S.index.by_primary_key_dict.get(key, None)
  if not probe:
    print("-- Dangling taxonID reference: %s (from %s)" %
          (get_primary_key(key, '?'), comes_from),
          file=sys.stderr)
  return probe

def link_to_superiors(S):
  seniors = 0
  for usage in S.index.by_primary_key_dict.values():
    parent_key = get_parent_key(usage, None)
    accepted_key = get_accepted_key(usage, None)
    if accepted_key:
      accepted_usage = lookup_usage(S, accepted_key, usage)
      if accepted_usage:
        S.set_accepted(usage, accepted_usage)
        # Filter out senior synonyms here
        if seniority(C, usage, accepted_usage) == "senior synonym":
          seniors += 1
          del S.index.by_primary_key_dict[S.get_primary_key(usage)] # ????
        else:
          if monitor(usage) or monitor(accepted_usage):
            print("> accepted %s := %s" %
                  (get_blurb(usage), get_blurb(accepted_usage)),
                  file=sys.stderr)
    elif parent_key:
      parent_usage = lookup_usage(S, parent_key, usage)
      if parent_usage:
        S.set_parent(usage, parent_usage)
        if monitor(usage) or monitor(parent_usage):
          print("> parent %s := %s" % (get_blurb(usage), get_blurb(parent_usage)),
                file=sys.stderr)
  if seniors > 0:     # Maybe interesting
    print("-- Suppressed %s senior synonyms" % seniors,
          file=sys.stderr)

# Is given synonym usage a senior synonym of its accepted usage?
# In the case of splitting, we expect the synonym to be a senior
# synonym of the item.

# We could also look to see if taxonomicStatus is 'senior synonym'.

# Hmm, if there's exactly one senior synonym, we should remember it as
# being the 'split from' taxon.

def seniority(C, item, accepted_item):
  year1 = S.get_year(item, None)  # the synonym
  year2 = S.get_year(accepted_item, None)
  if year1 and year2:
    if year1 <= year2:
      return "senior synonym"
    else:
      return "junior synonym"
  else:
    return "synonym"

# -----------------------------------------------------------------------------
# Collect parent/child and accepted/synonym relationships
# Set children and synonyms properties

# E.g. Glossophaga bakeri (in 1.7) split from Glossophaga commissarisi (in 1.6)

def collect_inferiors(C):
  items = C.index.by_primary_key_dict.values()
  roots = []
  for item in items:

    # Add item to list of accepted's synonyms
    accepted_item = C.get_accepted(item, None)
    if accepted_item:
      ch = C.get_synonyms(accepted_item, None)
      if ch:
        ch.append(item)
      else:
        C.set_synonyms(accepted_item, [item])

    else:
      # Add item to list of parent's children
      parent_item = C.get_parent(item, None)
      if parent_item:
        ch = C.get_children(parent_item, None)
        if ch:
          ch.append(item)
        else:
          C.set_children(parent_item, [item])
      else:
        roots.append(item)

  # We really only want one root (this is so that mrca can work)
  if True or len(roots) > 1:
    top = make_usage(TOP, C)
    C.set_canonical(top, TOP)
    C.index.by_primary_key_dict[TOP] = top
    for root in roots: C.set_parent(root, top)
    C.set_children(top, roots)
    roots = [top]
  else:
    top = roots[0]
  C.top = top
  C.roots = roots
  return C

# -----------------------------------------------------------------------------
# Report on differences between record matches and hierarchy matches
# r could be either a base usage record or a union

def get_blurb(r):
  if isinstance(r, prop.Instance):
    u = get_unique(r)      # string
    return u
  elif r:
    return "[not a dict]"
  else:
    return "[no match]"

global_unique_dict = {}  # human readable unique name

def pick_unique(x):
  for name in generate_names(x):
    if name in global_unique_dict:
      log("# %s is in global_unique_dict !?? %s" % name)
    else:
      global_unique_dict[name] = name
      return name

# Generate names similar to x's but globally unique

def generate_names(x):
  stem = get_canonical(x, None) or get_primary_key(x)
  yield stem
  check = get_name(get_source_checklist(x))
  yield "%s sec. %s" % (stem, check)
  while True:
    yield "%s (%x)" % (stem, b)
    yield "%s (%x) sec. %s" % (stem, b, check)

get_unique = prop.getter(prop.get_property("unique", filler=pick_unique))


def monitor(x):
  return (x and len(get_canonical(x, None)) == 1) #.startswith("Metachirus"))

# -----------------------------------------------------------------------------
# Generate csv rows for a checklist

def emit_rows(C, props):
  props = list(props)
  yield [prp.label for prp in props]
  getters = tuple(map(prop.getter, props))
  for usage in C.index.by_primary_key_dict.values():
    if not is_top(usage):
      yield [get(usage, prop.MISSING) for get in getters]

TOP = "[top]"

def is_top(x):
  return get_primary_key(x) == TOP

get_primary_key = prop.getter(primary_key_prop)

# -----------------------------------------------------------------------------

# Test with, say, src/newick.py --newick "(a,b)c" | src/checklist.py 

if __name__ == '__main__':
  src = load_source(csv.reader(sys.stdin), "A")

  props = (prop.get_property(label)
           for label in ("taxonID", "canonicalName", "parentNameUsageID"))
  writer = csv.writer(sys.stdout)
  for row in emit_rows(src, props): writer.writerow(row)

