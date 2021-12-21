#!/usr/bin/env python3

# Source and sum checklists

import sys, csv
import property as prop

from util import log

primary_key_prop = prop.get_property("taxonID")

canonical_prop = prop.get_property("canonicalName")
scientific_prop = prop.get_property("scientificName")
year_prop = prop.get_property("year")  # http://rs.tdwg.org/dwc/terms/year
rank_prop = prop.get_property("taxonRank")
managed_id_prop = prop.get_property("managed_id")
type_prop = prop.get_property("type")

source_prop = prop.get_property("source")    # which checklist does this belong to?
parent_key_prop = prop.get_property("parentNameUsageID")
accepted_key_prop = prop.get_property("acceptedNameUsageID")
parent_prop = prop.get_property("parent")
accepted_prop = prop.get_property("accepted")

# For workspaces
inject_prop = prop.get_property("inject") # Contextual only!
outject_prop = prop.get_property("outject")
equivalent_prop = prop.get_property("equivalent")

(get_primary_key, set_primary_key) = prop.get_set(primary_key_prop)
(get_parent_key, set_parent_key) = prop.get_set(parent_key_prop)
(get_accepted_key, set_accepted_key) = prop.get_set(accepted_key_prop)
(get_source, set_source) = prop.get_set(source_prop)
(get_canonical, set_canonical) = prop.get_set(canonical_prop)
(get_year, set_year) = prop.get_set(year_prop)
(get_managed_id, set_managed_id) = prop.get_set(managed_id_prop)

(get_parent, set_parent) = prop.get_set(parent_prop)
(get_accepted, set_accepted) = prop.get_set(accepted_prop)
(get_children, set_children) = prop.get_set(prop.get_property("children"))
(get_synonyms, set_synonyms) = prop.get_set(prop.get_property("synonyms"))

(get_outject, set_outject) = prop.get_set(outject_prop)
(get_equivalent, set_equivalent) = prop.get_set(equivalent_prop)


# -----------------------------------------------------------------------------
# Common code... for both alignments and checklists...

# every object using this must do  foo.context = Context()

def look_up_record(C, key, comes_from=None):
  if not key: return None
  col = prop.get_column(primary_key_prop, C.context) # we could cache this
  probe = prop.get_record(col, key, default=None)
  if not probe:
    print("-- Dangling taxonID reference: %s (from %s)" %
          (key, comes_from),
          file=sys.stderr)
  return probe

# -----------------------------------------------------------------------------
# Source checklists

class Source:
  def __init__(self, meta):
    self.context = prop.make_context()  # for lookup by primary key
    self.meta = meta

def all_records(C):
  col = prop.get_column(primary_key_prop, C.context)
  return prop.get_records(col)

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
  for record in prop.get_records(column):
    set_source(record, S)
  link_to_superiors(S)   # sets parents, accepteds, and S.roots
  collect_inferiors(S)
  return S

# Two ways to make these things.  One is using plans (with `construct`).
# The other is by using the custom constructor (`make_record`, two arguments).

make_record = prop.constructor(primary_key_prop, source_prop)

# Convert references to records rather than ids

def link_to_superiors(S):
  col = prop.get_column(primary_key_prop, S.context)
  for record in prop.get_records(col):
    parent_key = get_parent_key(record, None)
    accepted_key = get_accepted_key(record, None)
    if accepted_key:
      accepted_record = look_up_record(S, accepted_key, record)
      if accepted_record:
        S.set_accepted(record, accepted_record)
        if monitor(record) or monitor(accepted_record):
          print("> accepted %s := %s" %
                (blurb(record), blurb(accepted_record)),
                file=sys.stderr)
    elif parent_key:
      parent_record = look_up_record(S, parent_key, record)
      if parent_record:
        set_parent(record, parent_record)
        if monitor(record) or monitor(parent_record):
          print("> parent %s := %s" % (blurb(record), blurb(parent_record)),
                file=sys.stderr)

# -----------------------------------------------------------------------------
# Collect parent/child and accepted/synonym relationships;
# set children and synonyms properties.
# **** Checklist could be either source or coproduct. ****

# E.g. Glossophaga bakeri (in 1.7) split from Glossophaga commissarisi (in 1.6)

def collect_inferiors(C):
  roots = []

  for item in all_records(C):

    # Add item to list of accepted's synonyms
    accepted_item = get_accepted(item, None)
    if accepted_item:

      ch = get_synonyms(accepted_item, None)
      if ch:
        ch.append(item)
      else:
        set_synonyms(accepted_item, [item])

    else:
      # Add item to list of parent's children
      parent_item = get_parent(item, None)
      if parent_item:
        ch = get_children(parent_item, None)
        if ch:
          ch.append(item)
        else:
          set_children(parent_item, [item])
      else:
        roots.append(item)

  # We really only want one root (this is so that mrca can work)
  if True or len(roots) > 1:
    top = make_record(TOP, C)
    set_canonical(top, TOP)
    for root in roots: set_parent(root, top)
    set_children(top, roots)
    roots = [top]
  else:
    top = roots[0]
  C.top = top
  C.roots = roots
  return C

# -----------------------------------------------------------------------------
# For debugging

def blurb(r):
  if isinstance(r, prop.Record):
    return (get_canonical(r, None) or     # string
            get_scientific(r, None) or
            get_managed_id(r, None) or
            get_primary_key(r, None) or
            "[no identifier]")
  elif r:
    return "[not a record]"
  else:
    return "[no match]"

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
             prop.get_property("parentNameUsageID",
                               getter=get_parent_key),
             prop.get_property("acceptedNameUsageID",
                               getter=get_accepted_key),
             canonical_prop,
             scientific_prop,
             # year?  type?  ...
             rank_prop)
  yield [prp.label for prp in props]
  getters = tuple(map(prop.getter, props))
  for record in all_records(C):
    if not is_top(record):
      yield [get(record, prop.MISSING) for get in getters]


TOP = "[top]"

def is_top(x):
  return get_primary_key(x) == TOP


# -----------------------------------------------------------------------------

if __name__ == '__main__':
  import newick

  def testit(n):
    src = rows_to_checklist(newick.parse_newick(n),
                            {"name": "A"})  # meta
    writer = csv.writer(sys.stdout)
    rows = list(checklist_to_rows(src))
    for row in rows:
      writer.writerow(row)
    print(newick.compose_newick(rows))

  testit("a")
  testit("(a,b)c")
