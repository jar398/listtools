#!/usr/bin/env python3

# Source and sum checklists

import sys, csv
from typing import NamedTuple, Any

import property as prop
from util import log, MISSING
from rcc5 import *

# Strings (field values)
primary_key_prop = prop.get_property("taxonID")
canonical_prop = prop.get_property("canonicalName")
scientific_prop = prop.get_property("scientificName")
year_prop = prop.get_property("year")  # http://rs.tdwg.org/dwc/terms/year
rank_prop = prop.get_property("taxonRank")
managed_id_prop = prop.get_property("managed_id")
type_prop = prop.get_property("type")
parent_key_prop = prop.get_property("parentNameUsageID")
accepted_key_prop = prop.get_property("acceptedNameUsageID")

# Other checklist properties
source_prop = prop.get_property("source", inherit=False)    # which checklist does this belong to?
superior_prop = prop.get_property("superior", inherit=False)    # value is a Related
equated_prop = prop.get_property("equated", inherit=False)    # value is a Related

# For workspaces
inject_prop = prop.get_property("inject") # Contextual only!
outject_prop = prop.get_property("outject")

# For making use of name matching output
match_prop = prop.get_property("match", inherit=False)
match_key_prop = prop.get_property("match_id")
match_relation_prop = prop.get_property("relation")
match_note_prop = prop.get_property("match_note")

(get_primary_key, set_primary_key) = prop.get_set(primary_key_prop)
(get_parent_key, set_parent_key) = prop.get_set(parent_key_prop)
(get_accepted_key, set_accepted_key) = prop.get_set(accepted_key_prop)
(get_source, set_source) = prop.get_set(source_prop)
(get_canonical, set_canonical) = prop.get_set(canonical_prop)
(get_scientific, set_scientific) = prop.get_set(scientific_prop)
(get_rank, set_rank) = prop.get_set(rank_prop)
(get_year, set_year) = prop.get_set(year_prop)
(get_managed_id, set_managed_id) = prop.get_set(managed_id_prop)

(get_superior, set_superior) = prop.get_set(superior_prop)
(get_children, set_children) = prop.get_set(prop.get_property("children", inherit=False))
(get_synonyms, set_synonyms) = prop.get_set(prop.get_property("synonyms", inherit=False))
(get_taxonomic_status, set_taxonomic_status) = \
  prop.get_set(prop.get_property("taxonomicStatus"))
(get_equated, set_equated) = prop.get_set(equated_prop)

(get_outject, set_outject) = prop.get_set(outject_prop)

(get_match, set_match) = prop.get_set(match_prop)
get_match_key = prop.getter(match_key_prop)
get_match_relation = prop.getter(match_relation_prop)
get_match_note = prop.getter(match_note_prop)

class Related(NamedTuple):
  relation : Any    # < (LT, ACCEPTED), <= (LE, SYNONYM), =, NOINFO, maybe others
  record : Any       # the record that we're relating this one to
  status : str      # status (taxonomicStatus) relative to record ('synonym' etc)
  note : str = ''   # further comments justifying this relationship

# -----------------------------------------------------------------------------
# Common code... for both alignments and checklists...

# every object using this must do  foo.context = Context()

def look_up_record(C, key, comes_from=None):
  if not key: return None
  col = prop.get_column(primary_key_prop, C.context) # we could cache this
  probe = prop.get_record(col, key, default=None)
  if not probe:
    log("-- Dangling taxonID reference: %s (from %s)" %
        (key, comes_from))
  return probe

# -----------------------------------------------------------------------------
# Source checklists

class Source:
  def __init__(self, meta):
    self.context = prop.make_context()  # for lookup by primary key
    self.meta = meta
    self.indexed = False    # get_children, get_synonyms set?

def all_records(C):             # not including top
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
  Q = rows_to_context(iterabl, primary_key_prop)
  S.context = Q
  column = prop.get_column(primary_key_prop, Q)
  for record in prop.get_records(column):
    set_source(record, S)       # not same as EOL "source" column
  S.top = make_top(S)             # Superior of last resort
  resolve_superior_links(S)
  return S

def rows_to_context(row_iterable, primary_key_prop):
  Q = prop.make_context()
  register = prop.get_registrar(primary_key_prop, Q)
  row_iterator = iter(row_iterable)
  plan = prop.make_plan_from_header(next(row_iterator))
  for row in row_iterator:
    register(prop.construct(plan, row))
  return Q

# Two ways to make these things.  One is using plans (with `construct`).
# The other is by using the custom constructor (`make_record`, two arguments).

make_record = prop.constructor(primary_key_prop, source_prop)

# Create direct references to records, rather than leaving links as
# ids

def resolve_superior_links(S):
  for record in all_records(S): # not including top
    resolve_superior_link(record)

def resolve_superior_link(record):
  sup = get_superior(record, None)
  if sup != None: return sup
  parent_key = get_parent_key(record, None)
  accepted_key = get_accepted_key(record, None)
  S = get_source(record)
  if accepted_key:
    accepted_record = look_up_record(S, accepted_key, record)
    if accepted_record:
      status = get_taxonomic_status(record, "synonym")
      sup = Related(SYNONYM, accepted_record, status)
    else:
      sup = Related(SYNONYM, top, "dangling reference")
  elif parent_key:
    parent_record = look_up_record(S, parent_key, record)
    if parent_record:
      status = get_taxonomic_status(record, "accepted")
      # If it's not accepted or valid or something darn close, we're confused
      sup = Related(ACCEPTED, parent_record, status)
    else:
      sup = Related(ACCEPTED, top, "dangling reference")
  else:
    sup = Related(ACCEPTED, S.top, "root")
  set_superior_carefully(record, sup)
  if False and (monitor(record) or monitor(sup.record)):
    log("> %s %s := %s" % (sup.status, blurb(record), blurb(sup.record)))
  return sup

def set_superior_carefully(x, ship):
  assert_local(x, ship.record)
  set_superior(x, ship)

def assert_local(x, y):
  assert get_source(x) == get_source(y), \
    (blurb(x), get_source_name(x), blurb(y), get_source_name(y))

# -----------------------------------------------------------------------------
# Collect parent/child and accepted/synonym relationships;
# set children and synonyms properties.
# **** Checklist could be either source or coproduct. ****

# E.g. Glossophaga bakeri (in 1.7) split from Glossophaga commissarisi (in 1.6)

def ensure_inferiors_indexed(C):
  if C.indexed: return
  #log("# indexing inferiors")

  for record in all_records(C):
    assert get_source(record) == C
    if record == C.top: continue

    ship = get_superior(record, None)
    if not ship:
      # Default relation, if none given, is child of top ('orphan')
      # log("# %s is child of top" % blurb(record))
      ship = Related(ACCEPTED, C.top, "root")
      set_superior_carefully(record, ship)

    assert_local(record, ship.record)

    if ship.relation == ACCEPTED:   # accepted
      #log("# accepted %s -> %s" % (blurb(record), blurb(ship.record)))
      parent = ship.record
      # Add record to list of parent's children
      ch = get_children(parent, None) # list of records
      if ch != None:
        ch.append(record)
      else:
        set_children(parent, [record])

    else:         # synonym (LE)
      assert ship.relation == SYNONYM, rcc5_symbol(ship.relation)
      accepted = ship.record
      # Add record to list of accepted record's synonyms
      ch = get_synonyms(accepted, None)
      if ch != None:
        ch.append(record)
      else:
        set_synonyms(accepted, [record])
      if get_taxonomic_status(record, None) == "equivalent":
        # let's also note that superior is equivalent to record
        set_equated(accepted, Related(EQ, record,   "equivalent"))
        set_equated(record,   Related(EQ, accepted, "equivalent"))

  C.indexed = True
  #log("# Top has %s child(ren)" % len(get_children(C.top)))

def get_source_name(x):
  return get_source(x).meta['name']

def make_top(C):
  top = make_record(TOP, C)     # key is TOP, source is C
  # registrar(primary_key_prop)(top)  # Hmmph
  set_canonical(top, TOP)
  return top

TOP = "⊤"

def is_top(x):
  return x == get_source(x).top

def is_toplike(x):
  return get_canonical(x) == TOP

def add_child(c, x, status="accepted"):
  ch = get_children(x, None)
  if ch != None:
    ch.append(c)
  else:
    ch = [c]
  set_superior_carefully(c, Related(ACCEPTED, x, status))

def get_inferiors(x):
  e = get_equated(x, None)
  if e and e.relation == SYNONYM: yield e
  yield from get_children(x, ())
  for syn in get_synonyms(x, ()):
    y = get_equated(syn, None)
    def passing(m): clog(m, x, syn, y)
    # TBD: Also keep it if canonicalName differs
    if not y:
      pass #ing("keep because not equated to anything")
    elif y.record != x:
      passing("keep because equated but to wrong record")
    elif get_canonical(y.record) != get_canonical(syn):
      passing("keep because provides a different name")
    else:
      continue
    yield syn

# ----- Functions for filling columns in output table

def recover_parent_key(x, default=MISSING):
  sup = get_superior(x, None)
  if sup and not is_toplike(sup.record) and sup.relation == ACCEPTED:
    return get_primary_key(sup.record)
  else: return default

def recover_accepted_key(x, default=MISSING):
  sup = get_superior(x, None)
  if sup and not is_toplike(sup.record) and sup.relation != ACCEPTED:
    return get_primary_key(sup.record)
  else: return default

def recover_status(x, default=MISSING):
  sup = get_superior(x, default)
  if not sup:
    return default
  elif sup.status:
    return sup.status
  elif sup.relation == ACCEPTED:
    return "accepted"
  else:
    assert sup.relation == SYNONYM
    return "synonym"

# -----------------------------------------------------------------------------
# Convert a checklist to csv rows (as an iterable); inverse of rows_to_checklist, above

usual_props = \
    (primary_key_prop,
     prop.get_property("parentNameUsageID",
                       getter=recover_parent_key),
     prop.get_property("acceptedNameUsageID",
                       getter=recover_accepted_key),
     prop.get_property("taxonomicStatus",
                       getter=recover_status),
     canonical_prop,
     scientific_prop,
     rank_prop)

# Alternative ways to order the rows

def checklist_to_rows(C, props=None):
  (header, record_to_row) = begin_table(C, props)
  yield header
  for x in all_records(C):
    if not is_toplike(x):
      yield record_to_row(x)

def preorder(C, props=None):
  (header, record_to_row) = begin_table(C, props)
  yield header
  ensure_inferiors_indexed(C)
  def traverse(x):
    if not is_toplike(x):
      yield record_to_row(x)
    for c in get_inferiors(x):
      yield from traverse(c)
  yield from traverse(C.top)

def begin_table(C, props):
  if props == None: props = usual_props
  getters = tuple(map(prop.getter, props))
  def record_to_row(x):
    return [get(x, prop.MISSING) for get in getters]
  return ([prp.label for prp in props], record_to_row)

# -----------------------------------------------------------------------------
# For debugging

def clog(x, *records):
  log(x + " " + " ".join(map(blurb, records)))

def blurb(r):
  if isinstance(r, prop.Record):
    name = (get_canonical(r, None) or     # string
            get_scientific(r, None) or
            get_managed_id(r, None) or
            get_primary_key(r, None) or
            "[no identifier]")
    sup = get_superior(r, None)
    if sup and sup.status == SYNONYM:
      return name + "*"
    else:
      return name
  elif isinstance(r, Related):
    return "[%s %s]" % (rcc5_symbol(r.relation), blurb(r.record))
  elif isinstance(r, str):
    return "'%s'" % r           # kludge
  elif r:
    return "[not a record]"
  else:
    return "[no match]"

def monitor(x):
  return (x and len(get_canonical(x, None)) == 1) #.startswith("Metachirus"))

# -----------------------------------------------------------------------------

if __name__ == '__main__':
  import newick

  def testit(n):
    src = rows_to_checklist(newick.parse_newick(n),
                            {"name": "A"})  # meta
    writer = csv.writer(sys.stdout)
    if False:  # Test this if second fails
      rows = list(checklist_to_rows(src))
    else:
      rows = list(preorder(src))
    for row in rows:
      writer.writerow(row)
    log(' = ' + newick.compose_newick(rows))

  #testit("a")
  testit("(a,b)c")
