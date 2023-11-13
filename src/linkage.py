#!/usr/bin/env python3

# Record linkage
# This doesn't calculate a match score or probability, but has
# some of that spirit, for which see wikipedia record linkage article...

import util
import parse, rows
import simple
from checklist import *
from workspace import *
from parse import parse_name
from typify import compute_parts_score, homotypic_score, \
  compute_distance, compute_score, pick_better_record, \
  really_get_typification_uf

# For each record, get list of matching records.... hmm
# Future: read exceptions or whole thing from a file

# Definition of 'subproblem':
# If a node pair (u, v) doesn't have u in us and v in vs for some
# 'subproblem' (us, vs) then it is not going to be considered a
# potential match.

# 'pre_estimate' = LUB estimate from first pass, if this is second pass.


# Informational.  Requires a rewrite for exemplars.

def report_on_links(AB, subproblems):
  full_count = 0
  half_count = 0
  amb_count = 0
  for (key, (us, vs)) in subproblems.items():
    for u in us:
      link = get_link(u, None)
      if link == None:
        pass
      elif link == False:
        amb_count += 1
      else:
        link2 = get_link(link, None)
        if link2 is u:
          full_count += 1
        elif link2 == None:
          half_count += 1
    for v in vs:
      link = get_link(u, None)
      if link == None:
        pass
      elif link == False:
        amb_count += 1
      else:
        link2 = get_link(link, None)
        if link2 is v:
          pass                  # already counted
        elif link2 == None:
          half_count += 1
  log("#   Links: %s mutual, %s one-sided, %s ambiguous" %
      (full_count, half_count, amb_count))

# Find blocks/chunks, one per epithet

def find_subproblems(AB):
  (A_index, B_index) = \
    map(lambda CD: \
        index_by_some_key(CD,
                          # should use genus if epithet is missing
                          get_subproblem_key),
        (AB, swap(AB)))
  subprobs = {}
  for (val, us) in A_index.items():
    assert val != MISSING, blurb(us[0])
    vs = B_index.get(val, None)
    if vs != None:
      subprobs[val] = (us, vs)
      # for u in us: set_subproblem(u, subprob)
      # for v in vs: set_subproblem(v, subprob)
  AB.subproblems = subprobs
  return subprobs

# Returns dict value -> key

def index_by_some_key(AB, fn):
  index = {}
  for x in postorder_records(AB.A):
    u = AB.in_left(x)
    key = fn(u)
    #assert key  - MSW has scientificName = ?
    have = index.get(key, None)
    if have:
      have.append(u)
    else:
      index[key] = [u]
  return index

# Each subproblem covers a single epithet (or name, if higher taxon)
# z is in AB

def get_subproblem_key(z):
  parts = get_parts(z)
  ep = parts.epithet
  assert ep != None
  key = ep if ep else parts.genus
  assert key != None  # MSW has scientificName = '?' ...
  return key

# Phase this out?  What's it for?  find_estimate

def get_mutual_link(u, default=-19): # misplaced
  v = get_link(u, None)
  if v:
    u_back = get_link(v, default)
    if u_back == u:
      return v
  return default

def get_link(u, default=-19):
  uf = really_get_typification_uf(u, None)
  if uf:
    (xid, u2, v) = uf.payload()
    return v if separated(u, v) else u2
  return None

# -----------------------------------------------------------------------------
# Plumbing

def generate_linkage_report(AB):
  yield ('from', 'to', 'score')
  for x in all_records(AB.A):
    u = AB.in_left(x)
    v = get_link(u, None)
    if v:
      yield (blurb(u), blurb(v), compute_score(u, v))
    elif v == False:
      yield (blurb(u), 'ambiguous', MISSING)

import argparse
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Generate list of exemplars proposed for two checklists
    """)
  parser.add_argument('--A', help="the A checklist, as path name or -",
                      default='-')
  parser.add_argument('--B', help="the B checklist, as path name or -",
                      default='-')
  parser.add_argument('--test', action='store_true', help="run smoke tests")
  args=parser.parse_args()

  if args.test:
    print(compute_parts_score(parse_name("Sturnira angeli"),
                              parse_name("Sturnira magna"),
                              3))
  else:
    a_name = 'A'; b_name = 'B'
    a_path = args.A
    b_path = args.B
    with rows.open(a_path) as a_rows: # rows object
      with rows.open(b_path) as b_rows:
        # compute name matches afresh
        AB = ingest_workspace(a_rows.rows(), b_rows.rows(),
                              A_name=a_name, B_name=b_name)
        find_links(AB)
        report_gen = generate_linkage_report(AB)
        util.write_rows(report_gen, sys.stdout)
