#!/usr/bin/env python3

# Source and sum checklists

import sys, csv
from typing import NamedTuple, Any

import property as prop
from util import log
from rcc5 import *

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

(get_superior, set_superior) = prop.get_set(parent_prop)
(get_children, set_children) = prop.get_set(prop.get_property("children"))
(get_synonyms, set_synonyms) = prop.get_set(prop.get_property("synonyms"))
(get_taxonomic_status, set_taxonomic_status) = \
  prop.get_set(prop.get_property("taxonomicStatus"))

(get_outject, set_outject) = prop.get_set(outject_prop)
(get_equivalent, set_equivalent) = prop.get_set(equivalent_prop)

class Related(NamedTuple):
  relation : Any    # < (LT), <= (LE), =, NOINFO, maybe others
  other : Any       # the record that we're relating this one to
  status : str      # status (taxonomicStatus) relative to other ('synonym' etc)

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
    set_source(record, S)       # not same as EOL "source" column
  S.top = make_top(S)             # Superior of last resort
  resolve_superior_links(S)   # sets superior relateds
  collect_inferiors(S)
  return S

# Two ways to make these things.  One is using plans (with `construct`).
# The other is by using the custom constructor (`make_record`, two arguments).

make_record = prop.constructor(primary_key_prop, source_prop)

# Convert references to records rather than ids

def resolve_superior_links(S):
  topship = Related(LT, S.top, "root")

  for record in all_records(S):
    parent_key = get_parent_key(record, None)
    accepted_key = get_accepted_key(record, None)
    if accepted_key:
      accepted_record = look_up_record(S, accepted_key, record)
      if accepted_record:
        status = get_taxonomic_status(record, "synonym")
        set_superior(record, Related(LE, accepted_record, status))
        if monitor(record) or monitor(accepted_record):
          print("> accepted %s := %s" %
                (blurb(record), blurb(accepted_record)),
                file=sys.stderr)
      else:
        set_superior(record, Related(LE, top, "dangling reference"))
    elif parent_key:
      parent_record = look_up_record(S, parent_key, record)
      if parent_record:
        status = get_taxonomic_status(record, "accepted")
        # If it's not accepted or valid or something darn close, we're confused
        set_superior(record, Related(LT, parent_record, status))
        if monitor(record) or monitor(parent_record):
          print("> parent %s := %s" % (blurb(record), blurb(parent_record)),
                file=sys.stderr)
      else:
        set_superior(record, Related(LT, top, "dangling reference"))
    else:
      set_superior(record, topship)

# -----------------------------------------------------------------------------
# Collect parent/child and accepted/synonym relationships;
# set children and synonyms properties.
# **** Checklist could be either source or coproduct. ****

# E.g. Glossophaga bakeri (in 1.7) split from Glossophaga commissarisi (in 1.6)

def collect_inferiors(C):
  for record in all_records(C):
    ship = get_superior(record)

    if ship.relation == LE:     # synonym
      accepted = ship.other
      # Add record to list of accepted record's synonyms
      ch = get_synonyms(accepted, None)
      if ch:
        ch.append(record)
      else:
        set_synonyms(accepted, [record])

    elif ship.relation == EQ:   # accepted
      parent = ship.other
      # Add record to list of parent's children
      ch = get_children(parent, None)
      if ch:
        ch.append(record)
      else:
        set_children(parent, [record])

    else:
      assert "surprising superior relationship", rcc5_symbol(ship.relation)

  return C

def make_top(C):
  top = make_record(TOP, C)
  set_canonical(top, TOP)

TOP = "⊤"

def is_top(x):
  return get_primary_key(x) == TOP

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
  if props == None: props = usual_props
  yield [prp.label for prp in props]
  getters = tuple(map(prop.getter, props))
  for record in all_records(C):
    if not is_top(record):
      yield [get(record, prop.MISSING) for get in getters]

def _get_parent_key(x):
  sup = get_superior(x)
  if sup.relation == LT:
    return get_primary_key(sup.other)
  else: return MISSING
def _get_accepted_key(x):
  sup = get_superior(x)
  if sup.relation == LE:
    return get_primary_key(sup.other)
  else: return MISSING
def _get_status(x):
  sup = get_superior(x)
  if sup.status:
    return sup.status
  elif sup.relation == LE:
    return "synonym"
  elif sup.relation == LT:
    return "accepted"
  else:
    assert False

usual_props = \
            (primary_key_prop,
             prop.get_property("parentNameUsageID",
                               getter=_get_parent_key),
             prop.get_property("acceptedNameUsageID",
                               getter=_get_accepted_key),
             prop.get_property("taxonomicStatus",
                               getter=_get_status),
             canonical_prop,
             scientific_prop,
             # year?  type?  ...x
             rank_prop)

def preorder(C, props=None):
  yield [prp.label for prp in props]
  def traverse(x):
    for c in get_children(x): traverse(c)
    for c in get_synonyms(x): traverse(c)
    if not is_top(x):
      yield [get(x, prop.MISSING) for get in getters]
  traverse(C.top)

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
