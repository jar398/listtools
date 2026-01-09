#!/usr/bin/env python3

# Source and sum checklists

import sys, csv
from typing import NamedTuple, Any

import property as prop
from util import log, MISSING
from rcc5 import *
from ranks import ranks_dict     # not very abstract
import parse

# Strings (field values)
primary_key_prop = prop.declare_property("taxonID")
canonical_prop = prop.declare_property("canonicalName")
scientific_prop = prop.declare_property("scientificName")
authorship_prop = prop.declare_property("scientificNameAuthorship")
year_prop = prop.declare_property("namePublishedInYear")  # http://rs.tdwg.org/dwc/terms/
rank_prop = prop.declare_property("taxonRank")
managed_id_prop = prop.declare_property("managed_id")
stemmed_prop = prop.declare_property("canonicalStem")
taxonomic_status_prop = prop.declare_property("taxonomicStatus")
nomenclatural_status_prop = prop.declare_property("nomenclaturalStatus")

# Other checklist properties
source_prop = prop.declare_property("source", inherit=False)    # which checklist does this belong to?

# Links to other records, sometimes with explanation
parent_key_prop = prop.declare_property("parentNameUsageID", inherit=False)
accepted_key_prop = prop.declare_property("acceptedNameUsageID", inherit=False)
superior_note_prop = prop.declare_property("superior_note", inherit=False)
superior_prop = prop.declare_property("superior", inherit=False)    # value is a Relation

# For A/B identifications
equated_key_prop = prop.declare_property("equated_id", inherit=False)    # value is a Relation
equated_note_prop = prop.declare_property("equated_note", inherit=False)    # value is a Relation
equated_prop = prop.declare_property("equated", inherit=False)    # value is a Relation

# For record matches made by name(s)
match_prop = prop.declare_property("match", inherit=False)

# For workspaces
inject_prop = prop.declare_property("inject") # Contextual only!
outject_prop = prop.declare_property("outject")

(get_primary_key, set_primary_key) = prop.get_set(primary_key_prop)
(get_source, set_source) = prop.get_set(source_prop)     # checklist, A or B or AB
(get_canonical, set_canonical) = prop.get_set(canonical_prop)
(get_scientific, set_scientific) = prop.get_set(scientific_prop)
(get_authorship, set_authorship) = prop.get_set(authorship_prop)
(get_rank, set_rank) = prop.get_set(rank_prop)
(get_year, set_year) = prop.get_set(year_prop)
(get_managed_id, set_managed_id) = prop.get_set(managed_id_prop)
(get_stemmed, set_stemmed) = prop.get_set(stemmed_prop)
(get_nomenclatural_status, set_nomenclatural_status) = \
  prop.get_set(nomenclatural_status_prop)
(get_taxonomic_status, set_taxonomic_status) = \
  prop.get_set(taxonomic_status_prop)

# One column from the matches table?
(get_match, set_match) = prop.get_set(match_prop)

# Links
(get_parent_key, set_parent_key) = prop.get_set(parent_key_prop)
(get_accepted_key, set_accepted_key) = prop.get_set(accepted_key_prop)
get_superior_note = prop.getter(superior_note_prop)
(get_superior, set_superior) = prop.get_set(superior_prop)
(get_children, set_children) = prop.get_set(prop.declare_property("children", inherit=False))
(get_synonyms, set_synonyms) = prop.get_set(prop.declare_property("synonyms", inherit=False))

get_gn_full = prop.getter(prop.declare_property("gn_canonical_full"))
get_gn_stem = prop.getter(prop.declare_property("gn_canonical_stem"))
get_gn_auth = prop.getter(prop.declare_property("gn_authorship"))

(get_redundant, set_redundant) = \
  prop.get_set(prop.declare_property("redundant"))

# Merge related links
get_equated_key = prop.getter(equated_key_prop)
get_equated_note = prop.getter(equated_note_prop)
(get_equated, set_equated) = prop.get_set(equated_prop)

# Workspaces
(get_outject, set_outject) = prop.get_set(outject_prop)

# Records in the matches table
get_match_key = prop.getter(prop.declare_property("match_id", inherit=False))
get_match_direction = prop.getter(prop.declare_property("direction", inherit=False))
get_match_kind = prop.getter(prop.declare_property("kind", inherit=False))
get_basis_of_match = prop.getter(prop.declare_property("basis_of_match", inherit=False))
get_match_relationship = prop.getter(prop.declare_property("relationship", inherit=False))

# -----------------------------------------------------------------------------
# Relations and relationships

class Relation(NamedTuple):
  relationship : Any    # < (LT, HAS_PARENT), <= (LE), =, NOINFO, maybe others
  record : Any       # the record that we're relating this one to
  span : int
  note : Any          # comments justifying this relationship

def relation(ship, record, note=MISSING, span=None):
  assert isinstance(ship, int), ship
  assert note != None
  assert record == False or isinstance(record, prop.Record), blurb(record)
  if span == None:
    if ship == EQ: span = 0
    #elif ship == SYNONYM or MYNONYS: span = 1
    else: span = 2
  return Relation(ship, record, span, note)

# rel: origin -> target

def reverse_relation(origin, rel):
  assert isinstance(rel, Relation)
  return Relation(reverse_relationship(rel.relationship),
                  origin,
                  rel.span,
                  reverse_note(rel.note))

def reverse_articulation(art):
  (origin, rel) = art
  return (rel.record, Relation(reverse_relationship(rel.relationship),
                               origin,
                               rel.span,
                               reverse_note(rel.note)))

def compose_relations(rel1, rel2):
  assert rel1
  assert rel2
  assert (not get_workspace(rel1.record)) == (not get_workspace(rel2.record))
  if is_identity(rel1): return rel2
  if is_identity(rel2): return rel1
  # Consider special relationships for synonyms 
  # if rel1.span <= 1 and rel2.span <= 1: ...
  return Relation(compose_relationships(rel1.relationship, rel2.relationship),
                  rel2.record,
                  rel1.span + rel2.span,
                  compose_notes(rel1.note, rel2.note))

def is_identity(rel):
  return rel.relationship == EQ and rel.note == MISSING

# Algebra of 'notes'

def reverse_note(note):
  if note:
    if note.startswith("←"):
      return note[1:]
    else:
      return "←" + note
  else:
    return note

def compose_notes(note1, note2):
  if note1 == MISSING: return note2
  if note2 == MISSING: return note1
  if note1 == note2: return note1     # for LT/GT chains
  return "%s; %s" % (note1, note2)

# -----------------------------------------------------------------------------
# Source checklists

class Source:
  def __init__(self, meta):
    self.context = prop.make_context()  # for lookup by primary key
    self.meta = meta
    self.indexed = False    # get_children, get_synonyms set?
    make_top(self)          # Superior of last resort
    self.workspace = None

def all_records(C):             # not including top
  col = prop.get_column(primary_key_prop, C.context)
  # wait, this includes records that should get discarded. careful.
  yield from prop.get_records(col)

def all_records_inclusive(C):   # including top
  yield from all_records(C)
  yield C.top

def preorder_records(C):        # starting from top
  assert C.top
  def traverse(x):
    assert isinstance(x, prop.Record)
    yield x
    for c in get_inferiors(x):
      assert isinstance(c, prop.Record)
      yield from traverse(c)
  yield from traverse(C.top)

def postorder_records(C):       # ending with top
  assert C.top
  def traverse(x):
    for c in get_inferiors(x):
      yield from traverse(c)
    yield x
  yield from traverse(C.top)

def get_workspace(u):
  return get_source(u).workspace

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
  column = prop.get_column(primary_key_prop, Q)
  count = 0
  for record in prop.get_records(column):
    count += 1
    set_source(record, S)       # not same as EOL "source" column
    resolve_superior(S, record)
  roots = list(get_inferiors(S.top))
  if len(roots) > 1:
    log("-- %s roots" % len(roots))
  # Poor man's cycle detection
  count2 = 0
  for row in preorder_rows(S):
    count2 += 1
  if count != count2-1:
    log("!!! count disagreement: %s rows, %s reachable" % (count, count2))
  return S

# Property scoping context ...

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

def resolve_superior(S, record):
  if record == S.top: return None
  assert isinstance(record, prop.Record)
  assert not get_superior(record, None)
  sup = None

  # If it nontrivially has an accepted record, then it's a synonym
  accepted_key = get_accepted_key(record, None)
  taxonID = get_primary_key(record)
  assert taxonID
  assert accepted_key != taxonID, (accepted_key, taxonID,
                                   get_source_tag(record))  # see start.py
  if accepted_key:
    # set_accepted(record, False)
    accepted = look_up_record(S, accepted_key, record)
    if accepted:                        # Is there an accepted record (parent)?
      status = get_taxonomic_status(record, "synonym")    # default = synonym?
      if status == "equivalent":
        rel = relation(EQ, record, note="input")
        set_equated(accepted, rel)
        sup = relation(EQ, accepted, note="input")
      else:
        if monitor(record):
          log("# checklist: debug: %s %s %s %s" %
              (blurb(record), taxonID, accepted_key, blurb(accepted)))
        sup = relation(synonym_relationship(record, accepted),
                       accepted,
                       note="checklist (synonym)", span=1)
    else:
      log("# Synonym in %s with unresolvable accepted: %s -> %s" %
          (get_tag(S), blurb(record), accepted_key))
  else:
    # Accepted
    assert is_accepted(record), (blurb(record), get_taxonomic_status(record))
    parent_key = get_parent_key(record, None)
    parent = look_up_record(S, parent_key, record)
    if parent:
      sup = relation(HAS_PARENT, parent, note="checklist (child)")
    else:
      if parent_key:
        log("# accepted in %s but has unresolvable parent: %s = %s" %
            (get_tag(S), blurb(record), parent_key))
      # Child of top
      sup = relation(HAS_PARENT, S.top, note="inferior of top")

  if monitor(record):
    "# resolving sup %s = %s" % (blurb(record), blurb(sup))
  if sup:
    set_superior_carefully(record, sup)
    if is_top(sup.record):
      log("# Superior %s = %s" % (blurb(record), blurb(sup)))
  else:
    log("# No superior: %s" % blurb(record))

def synonym_relationship(record, accepted):
  return SYNONYM                # could be LT, LE, EQ

# sup might be top

def set_superior_carefully(x, sup):
  assert isinstance(x, prop.Record)
  assert isinstance(sup, Relation)
  assert sup.record != x
  have = get_superior(x, None)
  if have:
    if (have.relationship != sup.relationship or
        have.record != sup.record):
      log("**** Changing sup of %s from %a to %s" %
          (blurb(x), blurb(have), blurb(sup)))
  assert x != sup.record, (blurb(x), blurb(sup))        # no self-loops

  # If sup itself is a synonym, we're in trouble.
  # Unfortunately GBIF has a bunch of these.
  # We are forced to discard hierarchy information between synonyms.
  if sup.relationship != EQ:
    s = sup
    while not is_accepted(s.record):
      su = get_superior(s.record, None)
      if not su: break          # s.record = top
      s = su
    if s != sup:
      if False:
        # Apparently this happens too often to warrant logging
        log("** parent not accepted: %s -> %s ...-> %s" %
            (blurb(x), blurb(sup.record), blurb(s.record)))
      sup = s

  normalize_ranks(x, sup)

  # OK go for it.  sup might be top
  link_superior(x, sup)

def normalize_ranks(x, sup):
  if is_accepted(x):            # ignore synonym rank?
    r1 = get_rank(x, None)
    r2 = get_rank(sup.record, None)
    if r1 and r2:
      d1 = ranks_dict.get(r1)
      d2 = ranks_dict.get(r2)
      if d1 == None or d2 == None:
        pass
      if d1 < d2:
        pass
      elif d1 == d2:
        log("! Repeated rank: %s %s < %s %s" %
            (blurb(x), r1, blurb(sup.record), r2))
        set_rank(x, None)
      else:
        log("! Ranks out of order: %s %s < %s %s" %
            (blurb(x), r1, blurb(sup.record), r2))
        set_rank(x, None)

# Establish links w -> p and p -> w
# sup is a Relation with record p

def link_superior(w, sup):      # w is inferior Record, sup is Relation
  assert isinstance(w, prop.Record), blurb(w)
  assert isinstance(sup, Relation), blurb(sup)
  assert isinstance(sup.record, prop.Record), blurb(sup)
  # why would this fail?  (for newick)  sup.record is missing
  assert get_source(w) is get_source(sup.record), blurb(sup)

  #assert is_accepted(sup.record), blurb(sup)
  assert get_superior(w, None) == None
  set_superior(w, sup)
  link_inferior(w, sup)

def link_inferior(w, sup):
  if is_accepted(w):                    # based no taxonomicStatus
    ch = get_children(sup.record, None) # list of records
    if ch != None:
      ch.append(w)
    else:
      set_children(sup.record, [w]) # w is a Record
  else:
    ch = get_synonyms(sup.record, None)
    if ch != None:
      ch.append(w)
    else:
      set_synonyms(sup.record, [w])
  if is_top(sup.record):
    log("# Inferior of top: %s -> %s" % (blurb(w), blurb(sup)))

# What is comes_from for?  Potential debugging?

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
  for z in all_records(C):
    count += 1
    if not prop.mep_get(seen, z, False):
      log("** %s not in hierarchy; has to be a cycle" % (blurb(z),))
    sup = get_superior(z, None)
    if not is_accepted_locally(C, z):
      if len(get_children(z, ())) > 0:
        log("** synonym %s has child %s" % (blurb(z), blurb(get_children(z)[0])))
      if len(get_synonyms(z, ())) > 0:
        if monitor(z):
          log("> synonym %s has synonym %s" % (blurb(z), blurb(get_synonyms(z)[0])))
  log("# %s total records" % len(seen))

# -----------------------------------------------------------------------------
# Collect parent/child and accepted/synonym relationships;
# set children and synonyms properties.
# **** Checklist could be either source or coproduct. ****
# Run this AFTER all superior links have been set in other ways (e.g. merge).

# E.g. Glossophaga bakeri (in 1.7) split from Glossophaga commissarisi (in 1.6)

# This seems to be no longer used??
# Replaced by link_inferior?

def ensure_inferiors_indexed(C):
  #if get_children(C.top, None) == None: return

  assert False, "top has no children"

  count = 0

  for record in all_records(C):    # potential inferiors
    if record == C.top: continue
    # top has no superior, all others

    sup = get_superior(record, None)
    assert sup, blurb(record)
    assert_local(record, sup.record)

    if is_accepted_locally(C, record):
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

    else:                       # synonym of some kind
      accepted = sup.record
      # Add record to list of accepted record's synonyms
      ch = get_synonyms(accepted, None)
      if ch != None:
        ch.append(record)
      else:
        set_synonyms(accepted, [record])
      #log("# zzz %s, %s" % (blurb(accepted), blurb(record)))
      # if EQ then also set equivalent ??
      # but the target may be a ghost.
      count += 1

  #log("# %s inferiors indexed.  Top has %s child(ren)" %
  #   (count, len(get_children(C.top, ()))))
  validate(C)

# -----------------------------------------------------------------------------
# More

def assert_local(x, y):
  assert get_source(x) is get_source(y), \
    (blurb(x), get_source_tag(x), blurb(y), get_source_tag(y))

def get_source_tag(x):          # applied to a record
  assert len(x) >= 0  # record
  src = get_source(x)
  return get_tag(src)

def get_tag(check):          # applied to a Source checklist
  return check.meta['tag']

def make_top(C):
  top = make_record(TOP_NAME, C) # primary key is TOP_NAME, source is C
  set_canonical(top, TOP_NAME)   # name = key
  set_rank(top, "top")
  set_taxonomic_status(top, "accepted") # a lie
  C.top = top
  return top

TOP_NAME = "⊤"

# For non-coproduct checklists only
def is_top(x):
  # Why would get_source(x) not be defined? 'missing value'
  sou = get_source(x)
  return x == sou.top

def is_toplike(x):
  return get_canonical(x, None) == TOP_NAME

def get_inferiors(x):
  yield from get_children(x, ())
  yield from get_synonyms(x, ())

# Really this tests for not-synonym, but accepted is nicer to say
# and the doubtful case only occurs in certain sources (such as GBIF)

def is_accepted(x):             # exported
  status = get_taxonomic_status(x, "accepted") # default
  return (status.startswith("accepted") or
          status.startswith("provisionally accepted") or #GBIF
          status.startswith("valid") or
          status.startswith("doubtful"))

def get_accepted_relation(x):
  if is_accepted(x):
    return relation(EQ, x, MISSING)
  else:
    sup = get_superior(x, None)
    if sup:
      assert is_accepted(sup.record), blurb(x)
      return sup
    else:
      return relation(EQ, x, MISSING) # top

def get_accepted(x):            # Doesn't return falsish
  # return get_accepted_relation.record
  if is_accepted(x):
    return x
  else:
    rp = get_superior(x, None)
    if rp:
      return get_accepted(rp.record) # Eeeek!
    else:
      return x

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
  if is_accepted(x):
    sup = get_superior(x, None)
    if sup:
      y = sup.record
      if y:                     # might be top
        return get_primary_key(y, default)
  return default

def recover_accepted_key(x, default=MISSING):
  y = get_accepted(x)
  if y == x:
    return MISSING
  else:
    return get_primary_key(y, default)

def recover_status(x, default=MISSING): # taxonomicStatus
  sup = get_superior(x, None)
  if not sup:
    return default
  status = get_taxonomic_status(x, '?')
  if is_accepted(x):
    if status.startswith("accepted"):
      return status
    else:
      return "accepted " + status
  else:
    if status.startswith("accepted"):
      return "*" + status
    else:
      return status

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
  rel = get_matched(x)
  if rel: return get_primary_key(rel.record)
  else: return default

def recover_basis_of_match(x, default=MISSING):
  m = get_match(x, None)
  return m.note if m else default

def get_matched(x):
  if x:
    rel = get_match(x, None)
    if rel and rel.relationship == EQ:
      return rel
  return None

# -----------------------------------------------------------------------------
# Convert a checklist to csv rows (as an iterable); inverse of rows_to_checklist, above

# This is a basic set of Darwin Core only properties.
# When we serialize a checklist, these are the properties to be selectedd.
# Must contain all the properties to be used by match_records.

usual_props = \
    (primary_key_prop,
     canonical_prop,
     scientific_prop,
     stemmed_prop,
     managed_id_prop,
     rank_prop,
     prop.declare_property("parentNameUsageID",
                           getter=recover_parent_key),
     prop.declare_property("acceptedNameUsageID",
                           getter=recover_accepted_key),
     prop.declare_property("taxonomicStatus",
                           getter=recover_status))

# Alternative ways to order the rows

def checklist_to_rows(C, props=None):
  return records_to_rows(C, all_records(C), props)

def preorder_rows(C, props=None): # includes top
  return records_to_rows(C, preorder_records(C), props)

def records_to_rows(C, records, props=None): # Excludes top
  (header, record_to_row) = begin_table(C, props)
  yield header
  for x in records:
    if not x == C.top:
      assert isinstance(x, prop.Record)
      yield record_to_row(x)

# Filter out unnecessary equivalent A records!

def keep_record_notused(x):
  sup = get_superior(x, None)
  # if not sup: assert False
  if sup.relationship != EQ:
    return True

  # Hmm.  we're an A node that's EQ to some B node.
  # Are they also a record match?
  rel = get_matched(x)
  if (not rel) or rel.record != sup.record:
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
    if get_workspace(r):
      x = get_outject(r, None)
      name = get_source_tag(x) + "." + get_ok_name(x)
    else:
      x = r
      name = get_ok_name(x)
    if not is_accepted(x):
      name = name + "*"
    return name
  elif isinstance(r, Relation):
    if r.note:
      return ("[%s (%s) %s]" %
              (rcc5_symbol(r.relationship), r.note,
               blurb(r.record)))
    else:
      return "[%s %s]" % (rcc5_symbol(r.relationship), blurb(r.record))
  elif isinstance(r, str):
    return "'%s'" % r           # kludge
  elif r:
    return "[not a record]"
  elif r == None:
    return "[none]"
  elif r == False:
    return "[inconsistent]"     # typify.py
  else:
    return "[?]"

def blorb(u):
  if u == None:
    return "None"
  if get_workspace(u):
    x = get_outject(u)
  else:
    x = u
  prefix = get_source_tag(x) + "."
  name = get_scientific(x, "?") or get_canonical(x, "??")
  if not is_accepted(x):
    name = name + "*"
  return prefix + name

def monitor(x):
  if not x: return False
  name = get_ok_name(x)
  return (#name.startswith("Carpitalpa arendsi") or
          #name.startswith("Neophascogale lorentzi") or
          #name == "Rattus satarae" or
          #name == "Tragelaphus typicus" or
          #name == "Pogonomelomys mayeri" or
          #name == "Pseudochirops corinnae" or
          #name == "Pseudohydromys fuscus" or
          #name == "Burramys parvus" or
          #name == "Procyon cancrivorus" or
          #name == "Holochilus vulpinus" or
          #name == "Holochilus brasiliensis" or
          #name == "Holochilus darwini"
          # name == "Gorilla beringei"
          # name == "Sturnira magna"
          # "Mico marcai"
          # "Archaeolemur majori"
          False
          )

# -----------------------------------------------------------------------------
# Load/dump a set of provisional matches (could be either record match
# or taxonomic matches... but basically, record matches).  The matches are stored 
# as Relations under the 'match' property of nodes in AB.

# Obsolete

(get_match_info, set_match_info) = prop.get_set(prop.declare_property("match_info", inherit=False))

def load_matches(row_iterator, AB):
  log("# checklist: loading")

  # match_id,relationship,taxonID,direction,kind,basis_of_match
  header = next(row_iterator)
  plan = prop.make_plan_from_header(header)
  mutual_count = 0
  miss_count = 0

  for row in row_iterator:
    # row = [taxonID (x id), rel, matchID (y id), dir, kind, basis]
    match_record = prop.construct(plan, row)
    A_ids = get_primary_key(match_record, MISSING).split('|')
    xs = map(lambda id: look_up_record(AB.A, id), A_ids)
    us = [AB.in_left(x) for x in xs if x]
    B_ids = get_match_key(match_record, MISSING).split('|')
    ys = map(lambda id: look_up_record(AB.B, id), B_ids)
    vs = [AB.in_right(y) for y in ys if y]
    u = us[0] if len(us) == 1 else None
    v = vs[0] if len(vs) == 1 else None
    assert u or v, (A_ids, B_ids)

    # Put ambiguities into the comment ??

    ship = rcc5_relationship(get_match_relationship(match_record)) # EQ, NOINFO
    # match_direction, match_kind, basis_of_match
    dir = get_match_direction(match_record, '?')
    kind = get_match_kind(match_record, '?')
    basis = get_basis_of_match(match_record, '')
    note = "%s: %s" % (dir, kind)
           
    if monitor(u) or monitor(v):
      log("# match: %s -> %s, %s" %
          (tuple(map(blurb, us)), tuple(map(blurb, vs)), note))

    # u or v might be None with ship=NOINFO ... hope this is OK
    assert dir in ('A->B', 'B->A', 'A<->B')
    if dir == 'A<->B':
      assert ship == EQ
      mutual_count += 1
    else:
      assert ship == NOINFO
      miss_count += 1
    if dir == 'A<->B' or dir == 'A->B':
      assert u
      set_match(u, relation(ship, v, note=note))
      set_match_info(u, (vs, kind, basis))
    if dir == 'A<->B' or dir == 'B->A':
      assert v
      set_match(v, relation(ship, u, note=note))
      set_match_info(v, (us, kind, basis))

  log("# loaded %s match records, %s nonmatches" % (mutual_count, miss_count))

def get_matches(u):
  info = get_match_info(u, None)
  if info:
    (vs, kind, basis) = info
    return vs
  else:
    return ()

def diagnose_match(u):
  def descend(u, prefix, u_return, v_return, count):
    u_info = get_match_info(u, None)
    if u_info:
      (vs, u_kind, u_basis) = u_info
      mess = "%s%s %s %s" % (prefix, xblurb(u), u_kind, u_basis)
      if u == u_return:
        log(mess + " (return)")
      else:
        log(mess)
        if count > 0:
          for v in vs:
            descend(v, ' ' * len(prefix) + "  -> ", v_return, u_return or u, count-1)
  descend(u, "", None, None, 3)

"""
u1 -> v1 -> u2 -> v3
(None, None)
(u1, None)
(v1, u1)
(u1, v1)
...
"""

# Result of parsing = Parts
# Cached on workspace records

(get_parts_cache, set_parts_cache) = prop.get_set(prop.declare_property("parts_cache", inherit=False))

def get_parts(x):
  probe = get_parts_cache(x, None)
  if probe: return probe
  parts = parse.parse_name(get_best_name(x),
                           gn_full = get_gn_full(x, MISSING),
                           gn_stem = get_gn_stem(x, MISSING),
                           gn_auth = get_gn_auth(x, MISSING),
                           canonical = get_canonical(x, MISSING),
                           authorship = get_authorship(x, MISSING))
  set_parts_cache(x, parts)
  return parts

def get_tipe(x, default=None):
  parts = get_parts(x)
  return (parts.epithet, parts.token, parts.year)

def get_best_name(x):
  name = (get_scientific(x, None) or get_ok_name(x))
  assert name
  return name

def get_ok_name(x):
  name = (get_canonical(x, None) or
          get_managed_id(x, None) or get_primary_key(x))
  assert name
  return name

def xblurb(x):
  name = get_best_name(x)
  if not is_accepted(x):
    name = name + "*"
  stuff = [blurb(x)]
  stuff.append(get_primary_key(x))
  if get_nomenclatural_status(x, None):
    stuff.append(get_nomenclatural_status(x))
  return ' '.join(stuff)

if __name__ == '__main__':
  import newick

  def testit(n):
    src = rows_to_checklist(newick.parse_newick(n),
                            {'tag': 'A'})  # meta
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
