#!/usr/bin/env python3

# Source and sum checklists

import sys, csv
import property as prop
from property import mep_get, mep_set
from util import log
from coproduct import *

# -----------------------------------------------------------------------------
# Sum / coproduct / merged checklist

# Returns B side

def make_sum(A, B, rm_sum, kind):
  A_injections = {}
  B_injections = {}
  def join(x, y):
    if x: assert get_checklist(x) == A
    if y: assert get_checklist(y) == B
    z = (x and A_injections.get(x)) or (y and B_injections.get(y))
    if not z:
      z = make_injected(x, y)
      if x: mep_set(A_injections, x.id, z)
      if y: mep_set(B_injections, y.id, z)
    return z
  def split(z):
    return (out_A(z, None), out_B(z, None))
  co = Coproduct(in_A, in_B, out_A, out_B)
  co.rm_sum = rm_sum
  co.kind = kind
  co.index = Index()
  return co.AB

out_A_prop = get_property("out_A"); out_A = getter out_A_prop
out_B_prop = get_property("out_B"); out_B = getter out_B_prop

make_injected = prop.constructor(out_A_prop, out_B_prop)

checklist_prop = make_property("checklist")    # which checklist does this belong to?
get_checklist = prop.getter(checklist_prop)

# Checklist here is supposed to be a Side
# get_checklist is public and has to return Common

class Index:
  self.by_key_dict = {}     # from DwC taxonID column.  maybe others

""
  def record_match(self, x):
    # This doesn't feel right
    y = self.out_left(u)             # y in A
    r = self.out_right(self.rm_sum.in_left(y))
    return self.rm_sum.in_right(r) if r else None
""

# -----------------------------------------------------------------------------

class Source:
  def __init__(self, name):
    self.name = name
    self.index = Index()

# Read and write Darwin Core files
#   Stream of row <-> Source structure

# Darwin Core CSV format hierarchy file ingest

key_prop = prop.Property("taxonID")
get_key = prop.getter(key_prop)

make_usage = prop.constructor(key_prop, checklist_prop)

# Usually table columns

canonical_prop = prop.Property("canonicalName")
get_canonical = prop.getter(canonical_prop)
set_canonical = prop.setter(canonical_prop)

scientific_prop = prop.Property("scientificName")
get_scientific = prop.getter(scientific_prop)
set_scientific = prop.setter(scientific_prop)

rank_prop = prop.Property("taxonRank")
get_rank = prop.getter(rank_prop)
set_rank = prop.setter(rank_prop)

year_prop = prop.Property("year")
get_year = prop.getter(year_prop)
set_year = prop.setter(year_prop)

# Not usually table columns

parent_prop = prop.Property("parent") # Next bigger
get_parent = prop.getter(parent_prop)
set_parent = prop.setter(parent_prop)

accepted_prop = prop.Property("accepted")
get_accepted = prop.getter(accepted_prop)
set_accepted = prop.setter(accepted_prop)

unique_prop = prop.Property("unique")
get_unique = prop.getter(unique_prop)
set_unique = prop.setter(unique_prop)

def load_source(iterator, name):
  source = Source(name)

  header = next(iterator)
  plan = make_plan_from_header(header)
  get_parent_key = getter(get_property("parentNameUsageID"))
  get_accepted_key = getter(get_property("accedptedNameUsageID"))

  fixup = []
  for row in iterator:
    usage = prop.construct(plan, row)
    source.index.by_key_dict[usage.id] = usage
    key = get_key(usage)
    parent_key = get_parent_key(usage, None)
    accepted_key = get_accepted_key(usage, None)
    assert parent_key != key
    assert accepted_key != key

    source.index.by_key_dict[key] = usage
    fixup.append((usage, parent_key, accepted_key))

  seniors = 0
  for (usage, parent_key, accepted_key) in fixup:
    if accepted_key:
      accepted_usage = source.index.by_key_dict.get(accepted_key)
      if accepted_usage:
        set_accepted(usage, accepted_usage)
        # Filter out senior synonyms here
        if seniority(usage, accepted_usage) == "senior synonym":
          seniors += 1
          del source.index.by_key_dict[get_key(usage)] # ????
        else:
          if monitor(usage) or monitor(accepted_usage):
            print("> accepted %s := %s" %
                  (get_blurb(usage), get_blurb(accepted_usage)),
                  file=sys.stderr)
      else:
        print("-- Dangling accepted: %s -> %s" % (get_key(usage), accepted_key),
              file=sys.stderr)
    elif parent_key:
      parent_usage = source.index.by_key_dict.get(parent_key)
      if parent_usage:
        set_parent(usage, parent_usage)
        if monitor(usage) or monitor(parent_usage):
          print("> parent %s := %s" % (get_blurb(usage), get_blurb(parent_usage)),
                file=sys.stderr)
      else:
        print("-- Dangling parent: %s -> %s" % (get_key(usage), parent_key),
              file=sys.stderr)

  if seniors > 0:     # Maybe interesting
    print("Suppressed %s senior synonyms" % seniors,
          file=sys.stderr)

  return init_roots(source)

# -----------------------------------------------------------------------------

if __name__ == '__main__':
  load_source(csv.reader(sys.stdin), "A")

