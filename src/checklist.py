#!/usr/bin/env python3

# Source and sum checklists

import sys, csv
import property as prop

from coproduct import *
from util import log

primary_key_prop = prop.get_property("taxonID")
get_primary_key = prop.getter(primary_key_prop)

source_prop = prop.get_property("source")    # which checklist does this belong to?
(get_source, set_source) = prop.get_set(source_prop)

# -----------------------------------------------------------------------------
# Common code... for both alignments and checklists...

# every object using this must do  foo.context = Context()

def look_up_usage(C, key, comes_from=None):
  if not key: return None
  col = prop.get_column(primary_key_prop, C.context) # we could cache this
  probe = prop.get_instance(col, key, default=None)
  if not probe:
    print("-- Dangling taxonID reference: %s (from %s)" %
          (key, comes_from),
          file=sys.stderr)
  return probe

# -----------------------------------------------------------------------------
# Source checklists

class Source:
  def __init__(self, meta):
    self.context = None
    self.meta = meta
    init_checklist(self, None)
  def contains(self, x):
    return get_source(x) == self
  def all_usages(self):
    col = prop.get_column(primary_key_prop, self.context)
    return prop.get_instances(col)
  def get_usage(self, pk):
    g = prop.getter(primary_key_prop, context=self.context)
    return g(pk, None)

# -----------------------------------------------------------------------------
# Sum / coproduct / merged checklist
# Could be either just the matches (not a checklist), or the matches
# plus overrides (synthetic checklist)

# Returns B side

def make_sum(A, B, kind, rm_sum=None):
  def _in_left(x):
    assert A.contains(x)
    return x
  def _in_right(y):
    assert B.contains(y)
    return y
  def _case(z, when_left, when_right):
    if A.contains(z):
      return when_left(z)
    else:
      assert B.contains(z)
      return when_right(z)
  AB = Coproduct(_in_left, _in_right, _case)
  def contains(x):
    return A.contains(x) or B.contains(x)
  AB.contains = contains

  AB.kind = kind     # record match vs. extension match
  AB.A = A           # need ??
  AB.B = B

  # Usage coproducts are about matching

  AB.rm_sum = rm_sum
  if rm_sum:
    def record_match(x):
      assert False # TBD
      return rm_sum.get_match(x, None)
    AB.record_match = record_match

  Q = prop.make_context()       # allows overriding A and/or B
  init_checklist(AB, Q)    # we do this for any kind of checklist

  return AB

# Does this checklist contain this particular instance?

def contains(C, x):
  return AB.contains(x)

def record_match(AB, x):
  return AB.record_match(x)

"""
  def record_match(self, x):
    # This doesn't feel right
    y = self.out_left(u)             # y in A
    r = self.out_right(self.rm_sum.in_left(y))
    return self.rm_sum.in_right(r) if r else None
"""

# -----------------------------------------------------------------------------
# Common code for both kinds of checklist

# This approach makes checklists rather heavyweight... that's OK for now

canonical_prop = prop.get_property("canonicalName")
scientific_prop = prop.get_property("scientificName")
rank_prop = prop.get_property("taxonRank")
year_prop = prop.get_property("year")  # http://rs.tdwg.org/dwc/terms/year
parent_key_prop = prop.get_property("parentNameUsageID")
accepted_key_prop = prop.get_property("acceptedNameUsageID")
parent_prop = prop.get_property("parent")
accepted_prop = prop.get_property("accepted")
equivalent_prop = prop.get_property("match")
  
# Ambient
(get_source, set_source) = prop.get_set(source_prop)
(get_parent_key, set_parent_key) = prop.get_set(parent_key_prop)
(get_accepted_key, set_accepted_key) = prop.get_set(accepted_key_prop)
(get_year, set_year) = prop.get_set(year_prop)
(get_canonical, set_canonical) = prop.get_set(canonical_prop)

def init_checklist(C, Q):
  C.context = Q

  (C.get_match, C.set_match) = prop.get_set(equivalent_prop, Q)    # ???

  # Properties that are usually CSV context.  Ambient but can be
  # overridden contextually.
  (C.get_primary_key, C.set_primary_key) = \
    prop.get_set(primary_key_prop, context=Q)

  (C.get_canonical, C.set_canonical) = \
    prop.get_set(canonical_prop, context=Q)
  (C.get_scientific, C.set_scientific) = \
    prop.get_set(scientific_prop, context=Q)
  (C.get_rank, C.set_rank) = \
    prop.get_set(rank_prop, context=Q)
  # Not CSV context.  - tbd, use this ambiently as well?
  (C.get_source, C.set_checklist) = \
    prop.get_set(source_prop, context=Q)       # ?????
  (C.get_same, C.set_same) = \
    prop.get_set(prop.get_property("same"), context=Q)
  (C.get_parent, C.set_parent) = \
    prop.get_set(parent_prop, context=Q)
  (C.get_accepted, C.set_accepted) = \
    prop.get_set(accepted_prop, context=Q)
  (C.get_children, C.set_children) = \
    prop.get_set(prop.get_property("children"), context=Q)
  (C.get_synonyms, C.set_synonyms) = \
    prop.get_set(prop.get_property("synonyms"), context=Q)

  # Convenience methods
  def get_superior(x):
    return C.get_parent(x, None) or C.get_accepted(x, None)
  def get_inferiors(p):
    return C.get_children(p, []) + C.get_synonyms(p, [])
  C.get_superior = get_superior
  C.get_inferiors = get_inferiors

# -----------------------------------------------------------------------------
# Read and write Darwin Core files
#  (reading always yield source checklist not sum??)

#   Stream of row <-> Source structure

# Darwin Core CSV format hierarchy file ingest
# iterable -> source checklist

def rows_to_checklist(iterabl, meta):
  S = Source(meta)
  # Three passes: read, link, reverse link
  Q = prop.table_to_context(iterabl, primary_key_prop)
  S.context = Q
  column = prop.get_column(primary_key_prop, Q)
  for usage in prop.get_instances(column):
    set_source(usage, S)
  link_to_superiors(S)   # sets parents, accepteds, and S.roots
  collect_inferiors(S)
  return S

# Two ways to make these things.  One is using plans (with `construct`).
# The other is by using the custom constructor (`make_usage`, two arguments).

make_usage = prop.constructor(primary_key_prop, source_prop)

# Convert references to instances rather than ids

def link_to_superiors(S):
  col = prop.get_column(primary_key_prop, S.context)
  for usage in prop.get_instances(col):
    parent_key = get_parent_key(usage, None)
    accepted_key = get_accepted_key(usage, None)
    if accepted_key:
      accepted_usage = look_up_usage(S, accepted_key, usage)
      if accepted_usage:
        S.set_accepted(usage, accepted_usage)
        if monitor(usage) or monitor(accepted_usage):
          print("> accepted %s := %s" %
                (get_blurb(usage), get_blurb(accepted_usage)),
                file=sys.stderr)
    elif parent_key:
      parent_usage = look_up_usage(S, parent_key, usage)
      if parent_usage:
        S.set_parent(usage, parent_usage)
        if monitor(usage) or monitor(parent_usage):
          print("> parent %s := %s" % (get_blurb(usage), get_blurb(parent_usage)),
                file=sys.stderr)

# Is given synonym usage a senior synonym of its accepted usage?
# In the case of splitting, we expect the synonym to be a senior
# synonym of the item.

# We could also look to see if taxonomicStatus is 'senior synonym'.

# Hmm, if there's exactly one senior synonym, we should remember it as
# being the 'split from' taxon.

def is_senior(C, synonym_item, accepted_item):
  syn_year = S.get_year(synonym_item, None)  # the synonym
  acc_year = S.get_year(accepted_item, None)
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
# Load/dump sets of matches (either record or taxonomic)

def load_matches(row_iterator, A, B):

  AB = make_sum(A, B, "matches")

  (AB.get_equivalent, AB.set_equivalent) = \
    prop.get_set(equivalent_prop, context=AB.context)
  (AB.get_equivalence_note, AB.set_equivalence_note) = \
    prop.get_set(equivalence_note_prop, context=AB.context)

  header = next(iterator)
  plan = prop.make_plan_from_header(header)
  for row in row_iterator:
    # [equivalent taxonID remark]
    match = prop.construct(plan, row)

    xkey = get_equivalent_key(match, None)
    ykey = get_primary_key(match, None)
    remark = get_equivalence_note(match)

    x = look_up_usage(A, xkey) if xkey else None
    y = look_up_usage(B, ykey, x) if ykey else None

    if x:
      AB.set_equivalent(x, y)
      AB.set_equivalence_note(x, remark)
    if y:
      AB.set_equivalent(y, x)
      AB.set_equivalence_note(y, remark)

get_equivalent_key = prop.getter(prop.get_property("equivalentID"))
get_equivalence_note = prop.getter(prop.get_property("equivalence_note"))

equivalent_prop = prop.get_property("equivalent")
equivalence_note_prop = prop.get_property("equivalence_note")


# -----------------------------------------------------------------------------
# Collect parent/child and accepted/synonym relationships;
# set children and synonyms properties.
# **** Checklist could be either source or coproduct. ****

# E.g. Glossophaga bakeri (in 1.7) split from Glossophaga commissarisi (in 1.6)

def collect_inferiors(C):
  roots = []
  seniors = 0

  for item in C.all_usages():

    # Add item to list of accepted's synonyms
    accepted_item = C.get_accepted(item, None)
    if accepted_item:

      # Filter out senior synonyms here
      if is_senior(C, item, accepted_item):
        seniors += 1
      else:
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

  if seniors > 0:     # Maybe interesting
    print("-- Suppressed %s senior synonyms" % seniors,
          file=sys.stderr)

  # We really only want one root (this is so that mrca can work)
  if True or len(roots) > 1:
    top = make_usage(TOP, C)
    C.set_canonical(top, TOP)
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
  check = get_name(get_source(x))
  yield "%s sec. %s" % (stem, check)
  while True:
    yield "%s (%x)" % (stem, b)
    yield "%s (%x) sec. %s" % (stem, b, check)

get_unique = prop.getter(prop.get_property("unique", filler=pick_unique))


def monitor(x):
  return (x and len(get_canonical(x, None)) == 1) #.startswith("Metachirus"))

# -----------------------------------------------------------------------------
# Convert a checklist to csv rows (as an iterable); inverse of rows_to_checklist, above

# Version 1

def checklist_to_rows(C, props=None):
  if props == None:
    def get_parent_key(x):
      p = get_parent(x, None)
      return get_primary_key(p) if p else MISSING
    def get_accepted_key(x):
      p = get_accepted(x, None)
      return get_accepted_key(p) if p else MISSING
    props = (primary_key_prop,
             canonical_prop,
             scientific_prop,
             prop.get_property("parentNameUsageID",
                               getter=get_parent_key),
             prop.get_property("acceptedNameUsageID",
                               getter=get_accepted_key))
  yield [prp.label for prp in props]
  getters = tuple(map(prop.getter, props))
  for usage in C.all_usages():
    if not is_top(usage):
      yield [get(usage, prop.MISSING) for get in getters]



TOP = "[top]"

def is_top(x):
  return get_primary_key(x) == TOP



# -----------------------------------------------------------------------------

# Test with, say, src/newick.py --newick "(a,b)c" | src/checklist.py 

if __name__ == '__main__':
  import newick
  src = rows_to_checklist(newick.parse_newick("a"),
                          {"name": "A"})  # meta
  writer = csv.writer(sys.stdout)
  rows = list(checklist_to_rows(src))
  for row in rows:
    writer.writerow(row)
  print(newick.compose_newick(rows))
