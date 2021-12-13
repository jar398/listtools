#!/usr/bin/env python3

# Source and sum checklists

import sys, csv
import property as prop

from coproduct import *
from util import log

# -----------------------------------------------------------------------------
# Sum / coproduct / merged checklist

# Returns B side

def make_sum(A, B, kind, rm_sum=None):
  A_injections = {}
  B_injections = {}
  def join(x, y):
    if x: assert get_checklist(x) == A
    if y: assert get_checklist(y) == B
    z = (x and A_injections.get(x, None)) or (y and B_injections.get(y, None))
    if not z:
      z = make_injected(x, y)
      if x: prop.mep_set(A_injections, x.id, z)
      if y: prop.mep_set(B_injections, y.id, z)
    return z
  def split(z):
    return (out_A(z, None), out_B(z, None))
  co = Coproduct(join, split)
  co.rm_sum = rm_sum
  co.kind = kind
  co.index = Index()
  return co.AB

out_A_prop = prop.get_property("out_A"); out_A = prop.getter(out_A_prop)
out_B_prop = prop.get_property("out_B"); out_B = prop.getter(out_B_prop)

make_injected = prop.constructor(out_A_prop, out_B_prop)

checklist_prop = prop.get_property("checklist")    # which checklist does this belong to?
get_checklist = prop.getter(checklist_prop)

# Checklist here is supposed to be a Side
# get_checklist is public and has to return Common

class Index:
  def __init__(self):
    self.by_key_dict = {}     # from DwC taxonID column.  maybe others

"""
  def record_match(self, x):
    # This doesn't feel right
    y = self.out_left(u)             # y in A
    r = self.out_right(self.rm_sum.in_left(y))
    return self.rm_sum.in_right(r) if r else None
"""

# -----------------------------------------------------------------------------

class Source:
  def __init__(self, name):
    self.name = name
    self.index = Index()

# Read and write Darwin Core files
#   Stream of row <-> Source structure

# Darwin Core CSV format hierarchy file ingest

key_prop = prop.get_property("taxonID")
get_key = prop.getter(key_prop)

make_usage = prop.constructor(key_prop, checklist_prop)

# Usually table columns

canonical_prop = prop.get_property("canonicalName")
get_canonical = prop.getter(canonical_prop)
set_canonical = prop.setter(canonical_prop)

scientific_prop = prop.get_property("scientificName")
get_scientific = prop.getter(scientific_prop)
set_scientific = prop.setter(scientific_prop)

rank_prop = prop.get_property("taxonRank")
get_rank = prop.getter(rank_prop)
set_rank = prop.setter(rank_prop)

year_prop = prop.get_property("year")
get_year = prop.getter(year_prop)
set_year = prop.setter(year_prop)

# Not usually table columns

parent_prop = prop.get_property("parent") # Next bigger
get_parent = prop.getter(parent_prop)
set_parent = prop.setter(parent_prop)

accepted_prop = prop.get_property("accepted")
get_accepted = prop.getter(accepted_prop)
set_accepted = prop.setter(accepted_prop)

get_parent_key = prop.getter(prop.get_property("parentNameUsageID"))
get_accepted_key = prop.getter(prop.get_property("accedptedNameUsageID"))

def load_source(iterator, name):
  source = Source(name)

  header = next(iterator)
  plan = prop.make_plan_from_header(header)

  fixup = []
  for row in iterator:
    usage = prop.construct(plan, row)
    key = get_key(usage)
    source.index.by_key_dict[key] = usage
    parent_key = get_parent_key(usage, None)
    accepted_key = get_accepted_key(usage, None)
    assert parent_key != key
    assert accepted_key != key

    assert usage
    fixup.append((usage, parent_key, accepted_key))

  def get_usage(key, usage):
    probe = source.index.by_key_dict.get(key, None)
    if not probe:
      print("-- Dangling taxonID: %s -> %s" % (get_key(usage), key),
            file=sys.stderr)
    return probe

  seniors = 0
  for (usage, parent_key, accepted_key) in fixup:
    assert usage
    if accepted_key:
      accepted_usage = get_usage(accepted_key, usage)
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
    elif parent_key:
      parent_usage = get_usage(parent_key, usage)
      if parent_usage:
        set_parent(usage, parent_usage)
        #if monitor(usage) or monitor(parent_usage):
        print("> parent %s := %s" % (get_blurb(usage), get_blurb(parent_usage)),
                file=sys.stderr)

  if seniors > 0:     # Maybe interesting
    print("-- Suppressed %s senior synonyms" % seniors,
          file=sys.stderr)

  return init_roots(source)

TOP = "[top]"

def init_roots(source):
  # Collect children so we can say children[x]
  roots = collect_inferiors(source.index.by_key_dict.values())

  # We really only want one root (this is so that mrca can work)
  if True or len(roots) > 1:
    top = make_usage(TOP, source)
    # pick_unique(top)
    set_canonical(top, TOP)
    source.index.by_key_dict[TOP] = top
    for root in roots: set_parent(root, top)
    set_children(top, roots)
    roots = [top]
  else:
    top = roots[0]
  source.top = top

  # Prepare for doing within-tree MRCA operations
  cache_levels(roots)

  source.roots = roots
  return source

# Cache every node's level (distance to root)
#   simple recursive descent from roots

# Level is contravariant with RCC5: x < y implies level(x) > level(y)

level_prop = prop.get_property("level")
get_level = prop.getter(level_prop)
set_level = prop.setter(level_prop)

def cache_levels(roots):
  def cache(x, n):
    set_level(x, n)
    for c in get_inferiors(x):
      cache(c, n+1)
  for root in roots:
    cache(root, 1)

# -----------------------------------------------------------------------------
# Collect parent/child and accepted/synonym relationships

# E.g. Glossophaga bakeri (in 1.7) split from Glossophaga commissarisi (in 1.6)

def collect_inferiors(items):
  roots = []
  for item in items:

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

  return roots


(get_children, set_children) = prop.get_set(prop.get_property("children"))
(get_synonyms, set_synonyms) = prop.get_set(prop.get_property("synonyms"))

def get_inferiors(p):
  return (get_children(p, []) +
          get_synonyms(p, []))

def get_superior(x):
  sup = get_parent(x, None) or get_accepted(x, None)
  if sup:
    assert get_checklist(x) is get_checklist(sup)
    if get_level(x, None):
      if get_level(sup, 999) != get_level(x, 999) - 1:
        print("# %s %s" % (get_level(sup), get_level(x)),
              file=sys.stderr)
        assert False
  return sup

# -----------------------------------------------------------------------------

by_unique_dict = {}  # human readable unique name

def pick_unique(x):
  for name in generate_names(x):
    if name in by_unique_dict:
      log("# %s is in by_unique_dict !?? %s" % name)
    else:
      by_unique_dict[name] = name
      return name

# Generate names similar to x's but globally unique

def generate_names(x):
  stem = get_canonical(x, None) or get_key(x)
  yield stem
  check = get_name(get_checklist(x))
  yield "%s sec. %s" % (stem, check)
  while True:
    yield "%s (%x)" % (stem, b)
    yield "%s (%x) sec. %s" % (stem, b, check)

unique_prop = prop.get_property("unique", filler=pick_unique)
get_unique = prop.getter(prop.get_property("unique"))

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


def monitor(x):
  return (x and len(get_canonical(x, None)) == 1) #.startswith("Metachirus"))

# Generate csv rows for a checklist

def generate_rows(cl, props):
  instance_generator = cl.index.by_key_dict.values()
  prop.generate_rows(instance_generator, props)


# -----------------------------------------------------------------------------

if __name__ == '__main__':
  src = load_source(csv.reader(sys.stdin), "A")

  props = (prop.get_property(label) for label in ("taxonID", "canonicalName", "parentNameUsageID"))
  writer = csv.writer(sys.stdout)
  for row in generate_rows(src, props): writer.writerow(row)
