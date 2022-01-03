#!/usr/bin/env python3

# Source and sum checklists

import sys, csv
from typing import NamedTuple, Any

import property as prop
from util import log, MISSING
from rcc5 import *

# Strings (field values)
primary_key_prop = prop.declare_property("taxonID")
canonical_prop = prop.declare_property("canonicalName")
scientific_prop = prop.declare_property("scientificName")
year_prop = prop.declare_property("year")  # http://rs.tdwg.org/dwc/terms/year
rank_prop = prop.declare_property("taxonRank")
managed_id_prop = prop.declare_property("managed_id")
type_prop = prop.declare_property("type") # a.k.a. 'tipe'

# Other checklist properties
source_prop = prop.declare_property("source", inherit=False)    # which checklist does this belong to?

# Links to other records, sometimes with explanation
parent_key_prop = prop.declare_property("parentNameUsageID", inherit=False)    # LT
accepted_key_prop = prop.declare_property("acceptedNameUsageID", inherit=False) # LE
superior_note_prop = prop.declare_property("superior_note", inherit=False)
superior_prop = prop.declare_property("superior", inherit=False)    # value is a Related

# For A/B identifications
equated_key_prop = prop.declare_property("equated_id", inherit=False)    # value is a Related
equated_note_prop = prop.declare_property("equated_note", inherit=False)    # value is a Related
equated_prop = prop.declare_property("equated", inherit=False)    # value is a Related

# For record matches made by name(s)
match_key_prop = prop.declare_property("match_id", inherit=False)
basis_of_match_prop = prop.declare_property("basis_of_match", inherit=False)
match_prop = prop.declare_property("match", inherit=False)
get_match_relationship = prop.getter(prop.declare_property("relationship", inherit=False))

# For workspaces
inject_prop = prop.declare_property("inject") # Contextual only!
outject_prop = prop.declare_property("outject")

(get_primary_key, set_primary_key) = prop.get_set(primary_key_prop)
(get_source, set_source) = prop.get_set(source_prop)
(get_canonical, set_canonical) = prop.get_set(canonical_prop)
(get_scientific, set_scientific) = prop.get_set(scientific_prop)
(get_rank, set_rank) = prop.get_set(rank_prop)
(get_year, set_year) = prop.get_set(year_prop)
(get_managed_id, set_managed_id) = prop.get_set(managed_id_prop)

# Links
(get_parent_key, set_parent_key) = prop.get_set(parent_key_prop)
(get_accepted_key, set_accepted_key) = prop.get_set(accepted_key_prop)
(get_taxonomic_status, set_taxonomic_status) = \
  prop.get_set(prop.declare_property("taxonomicStatus"))
get_superior_note = prop.getter(superior_note_prop)
(get_superior, set_superior) = prop.get_set(superior_prop)
(get_children, set_children) = prop.get_set(prop.declare_property("children", inherit=False))
(get_synonyms, set_synonyms) = prop.get_set(prop.declare_property("synonyms", inherit=False))

# Merge related links
get_equated_key = prop.getter(equated_key_prop)
get_equated_note = prop.getter(equated_note_prop)
(get_equated, set_equated) = prop.get_set(equated_prop)

get_match_key = prop.getter(match_key_prop)
get_basis_of_match = prop.getter(basis_of_match_prop)
(get_match, set_match) = prop.get_set(match_prop)

(get_outject, set_outject) = prop.get_set(outject_prop)

class Relative(NamedTuple):
  relationship : Any    # < (LT, ACCEPTED), <= (LE, SYNONYM), =, NOINFO, maybe others
  record : Any       # the record that we're relating this one to
  status : str      # status (taxonomicStatus) relative to record ('synonym' etc)
  note : str = ''   # further comments justifying this relationship

def relation(ship, record, status, note=''):
  assert ((ship == NOINFO and not record) or \
          isinstance(record, prop.Record))
  assert isinstance(ship, int)
  return Relative(ship, record, status, note)

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

def preorder_records(C):
  ensure_inferiors_indexed(C)
  def traverse(x):
    yield x
    for c in get_inferiors(x):
      yield from traverse(c)
  yield from traverse(C.top)

# -----------------------------------------------------------------------------
# Read and write Darwin Core files
#  (reading always yield source checklist not sum??)

#   Stream of row <-> Source structure

# Darwin Core CSV format hierarchy file ingest
# iterable -> source checklist

def rows_to_checklist(iterabl, meta):
  S = Source(meta)
  # Three passes: read, link, reverse link
  (Q, register) = rows_to_context(iterabl, primary_key_prop)
  S.context = Q
  S.register = register
  S.top = make_top(S)             # Superior of last resort
  column = prop.get_column(primary_key_prop, Q)
  for record in prop.get_records(column):
    set_source(record, S)       # not same as EOL "source" column
    resolve_superior_link(S, record)
  return S

def rows_to_context(row_iterable, primary_key_prop):
  Q = prop.make_context()
  register = prop.get_registrar(primary_key_prop, Q)
  row_iterator = iter(row_iterable)
  plan = prop.make_plan_from_header(next(row_iterator))
  for row in row_iterator:
    assert (isinstance(row, tuple) or isinstance(row, list)), row
    register(prop.construct(plan, row))
  return (Q, register)

# Two ways to make these things.  One is using plans (with `construct`).
# The other is by using the custom constructor (`make_record`, two arguments).

make_record = prop.constructor(primary_key_prop, source_prop)

# Create direct references to records, rather than leaving links as
# ids.  Sets a superior for every record other than top.

def resolve_superior_link(S, record):
  if record == S.top: return None

  assert not get_superior(record, None)
  sup = None

  accepted_key = get_accepted_key(record, None)
  if accepted_key and accepted_key != get_primary_key(record):
    accepted = look_up_record(S, accepted_key, record)
    if accepted:
      status = get_taxonomic_status(record, "synonym")
      if status == "equivalent":
        rel = relation(EQ, record, "equivalent")
        set_equated(accepted, rel)
        if False:
          log("# Set equated for %s to %s" %
              (blurb(accepted), blurb(rel)))
        sup = relation(EQ, accepted, status)
      else:
        sup = relation(SYNONYM, accepted, status)
    else:
      log("# Dangling accepted reference in %s: %s -> %s" %
          (S.meta['name'], blurb(record), accepted_key))
      sup = relation(SYNONYM, S.top, "root", "unresolved accepted id")
  else:
    parent_key = get_parent_key(record, None)
    if parent_key:
      parent = look_up_record(S, parent_key, record)
      if parent:
        status = get_taxonomic_status(record, "accepted")
        # If it's not accepted or valid or something darn close, we're confused
        sup = relation(ACCEPTED, parent, status)
      else:
        log("# Dangling parent reference in %s: %s -> %s" %
            (S.meta['name'], blurb(record), parent_key))
        sup = relation(ACCEPTED, S.top, "unresolved parent id")
    else:
      # This is a root.  Hang on to it so that preorder can see it.
      assert isinstance(record, prop.Record)
      #log("# No accepted, no parent: %s '%s' '%s'" %
       #   (blurb(record), parent_key, accepted_key))
      sup = relation(ACCEPTED, S.top, "root", "no parent")

  assert sup.record != record
  set_superior_carefully(record, sup)

def set_superior_carefully(x, sup):
  have = get_superior(x, None)
  if have:
    if False:                   # GBIF fails if we insist
      assert have.record == sup.record, \
        (blurb(x), blurb(get_superior(x)), blurb(sup)) # record
    if have.relationship != sup.relationship:
      log("**** Changing sup of %s from %a to %s" %
          (blurb(x), blurb(have), blurb(sup)))
  assert x != sup.record, (blurb(x), blurb(sup))        # no self-loops

  # If sup itself is a synonym, we're in trouble
  if sup.relationship != EQ:
    grand = get_superior(sup.record, None)
    if grand and grand.relationship != ACCEPTED:
      log("** synonym of synonym: %s %s %s" %
          (blurb(x), blurb(sup), blurb(grand)))
      assert False

  # OK go for it
  set_superior(x, sup)
  if False:
    if (monitor(x) or monitor(sup.record)):  #False and 
      log("> superior of %s is %s" % (blurb(x), blurb(sup)))

def look_up_record(C, key, comes_from=None):
  if not key: return None
  col = prop.get_column(primary_key_prop, C.context) # we could cache this
  return prop.get_record(col, key, default=None)

# -----------------------------------------------------------------------------
# Check for children of synonyms, cycles, etc.

def validate(C):

  seen = prop.mep()
  def traverse(x):
    assert not prop.mep_get(seen, x, False)
    prop.mep_set(seen, x, True)
    for c in get_inferiors(x):
      traverse(c)
  traverse(C.top)
  log("# %s records reachable from top" % len(seen))

  count = 0
  for x in all_records(C):
    count += 1
    if not prop.mep_get(seen, x, False):
      log("** %s not in hierarchy; has to be a cycle" % (blurb(x),))
    sup = get_superior(x, None)
    if sup and sup.relationship == SYNONYM:
      if len(get_children(x, ())) > 0:
        log("** synonym %s has child %s" % (blurb(x), blurb(get_children(x)[0])))
      if len(get_synonyms(x, ())) > 0:
        if monitor(x):
          log("> synonym %s has synonym %s" % (blurb(x), blurb(get_synonyms(x)[0])))
  log("# %s total records" % len(seen))

# -----------------------------------------------------------------------------
# Collect parent/child and accepted/synonym relationships;
# set children and synonyms properties.
# **** Checklist could be either source or coproduct. ****
# Run this AFTER all superior links have been set in other ways (e.g. merge).

# E.g. Glossophaga bakeri (in 1.7) split from Glossophaga commissarisi (in 1.6)

def ensure_inferiors_indexed(C):
  if C.indexed: return
  #log("# --indexing inferiors")

  count = 0

  for record in all_records(C):
    if record == C.top: continue
    # top has no superior, all others

    assert get_source(record) == C
    sup = get_superior(record, None)
    assert sup, blurb(record)
    assert_local(record, sup.record)

    if sup.relationship == ACCEPTED:
      #log("# accepted %s -> %s" % (blurb(record), blurb(sup.record)))
      parent = sup.record
      # Add record to list of parent's children
      ch = get_children(parent, None) # list of records
      if ch != None:
        ch.append(record)
      else:
        set_children(parent, [record])
        #log("# child of %s is %s" % (blurb(parent), blurb(record)))
      count += 1

    else:                       # SYNONYM or EQ
      accepted = sup.record
      # Add record to list of accepted record's synonyms
      ch = get_synonyms(accepted, None)
      if ch != None:
        ch.append(record)
      else:
        set_synonyms(accepted, [record])
      #log("# %s of %s is %s" % (sup.status, blurb(accepted), blurb(record)))
      # if EQ then also set equivalent ??
      # but the target may be a ghost.
      count += 1

  C.indexed = True
  #log("# %s inferiors indexed.  Top has %s child(ren)" %
  #   (count, len(get_children(C.top, ()))))
  validate(C)

def assert_local(x, y):
  assert get_source(x) == get_source(y), \
    (blurb(x), get_source_name(x), blurb(y), get_source_name(y))

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
  return get_canonical(x, None) == TOP

# Not used I think
def add_inferior(r, status="accepted"):
  if r.relationship == LT:
    ch = get_children(x, None)
    s = set_children
  else:
    ch = get_synoyms(x, None)
    s = set_synonyms
  if ch != None:
    ch.append(r.record)
  else:
    s(x, [r.record])
  set_superior_carefully(c, r)

def get_inferiors(x):
  yield from get_children(x, ())
  yield from get_synonyms(x, ())


"""
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
"""

# ----- Functions for filling columns in output table

def recover_parent_key(x, default=MISSING):
  sup = get_superior(x, None)
  if sup and sup.relationship == ACCEPTED:
    y = get_accepted(sup.record)
    if not is_top(y):
      return get_primary_key(y)
  return default

def recover_accepted_key(x, default=MISSING):
  sup = get_superior(x, None)
  if sup and sup.relationship != ACCEPTED:
    y = get_accepted(sup.record)
    if not is_top(y):
      return get_primary_key(y)
  return default

def get_accepted(x):
  rp = get_superior(x, None)
  if rp and rp.relationship != ACCEPTED:
    if monitor(x):
      log("> snap link %s %s" % (blurb(x), blurb(rp)))
    return get_accepted(rp.record) # Eeeek!
    # return rp.record
  return x

def recover_status(x, default=MISSING): # taxonomicStatus
  sup = get_superior(x, None)
  if not sup:
    return default
  elif sup.status:
    return sup.status
  elif sup.relationship == ACCEPTED:
    return "accepted"
  else:
    assert sup.relationship == SYNONYM
    return "synonym"

# Non-DwC relationships, use optionally, see workspace

def recover_equated_key(x, default=MISSING):
  m = get_equated(x, None)
  if m and m.record:
    return get_primary_key(m.record)
  else: return default

def recover_equated_note(x, default=MISSING):
  m = get_equated(x, None)
  return m.note if m else default

def recover_match_key(x, default=MISSING):
  m = get_matched(x)
  if m: return get_primary_key(m)
  else: return default

def recover_basis_of_match(x, default=MISSING):
  m = get_match(x, None)
  return m.note if m else default

def get_matched(x):
  if x:
    rel = get_match(x, None)
    if rel and rel.relationship == EQ:
      return rel.record
  return None

# -----------------------------------------------------------------------------
# Convert a checklist to csv rows (as an iterable); inverse of rows_to_checklist, above

# This is a basic set of Darwin Core only properties.

usual_props = \
    (primary_key_prop,
     canonical_prop,
     scientific_prop,
     rank_prop,
     prop.declare_property("parentNameUsageID",
                       getter=recover_parent_key),
     prop.declare_property("acceptedNameUsageID",
                       getter=recover_accepted_key),
     prop.declare_property("taxonomicStatus",
                       getter=recover_status),
     managed_id_prop)

# Alternative ways to order the rows

def checklist_to_rows(C, props=None):
  return records_to_rows(C, all_records(C), props)

def preorder_rows(C, props=None):
  return records_to_rows(C, preorder_records(C), props)

def records_to_rows(C, records, props=None):
  (header, record_to_row) = begin_table(C, props)
  yield header
  for x in records:
    if not is_toplike(x):       # ?
      yield record_to_row(x)

# Filter out unnecessary equivalent A records!

def keep_record_notused(x):
  sup = get_superior(x, None)
  # if not sup: assert False
  if sup.relationship != EQ:
    return True

  # Hmm.  we're an A node that's EQ to some B node.
  # Are they also a record match?
  m = get_matched(x)
  if not m or m.record != sup.record:
    # If record is unmatched, or matches something not equivalent,
    # then keep it
    return True
  # They're a record match, are they also a name match?
  can1 = get_canonical(x)       # in A
  can2 = get_canonical(m.record) # in B
  if can1 and can1 != can2:
    # Similarly, keep if canonical differs
    return True
  # Can't figure out a way for them to be different.  Flush it.
  return False

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
    if sup and sup.relationship != ACCEPTED:
      return name + "*"
    else:
      return name
  elif isinstance(r, Relative):
    if r.note:
      return ("[%s %s '%s']" %
              (rcc5_symbol(r.relationship),
               blurb(r.record), r.note))
    else:
      return "[%s %s]" % (rcc5_symbol(r.relationship), blurb(r.record))
  elif isinstance(r, str):
    return "'%s'" % r           # kludge
  elif r:
    return "[not a record]"
  else:
    return "[no match]"

def monitor(x):
  if not x: return False
  name = get_canonical(x, '')
  return ((x and len(name) == 1)
          or name.startswith("Muri")
          )

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
      rows = list(preorder_rows(src))
    for row in rows:
      writer.writerow(row)
    log(' = ' + newick.compose_newick(rows))

  #testit("a")
  testit("(a,b)c")
