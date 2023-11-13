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
  compute_distance, compute_score, pick_better_record

# For each record, get list of matching records.... hmm
# Future: read exceptions or whole thing from a file

# Definition of 'subproblem':
# If a node pair (u, v) doesn't have u in us and v in vs for some
# 'subproblem' (us, vs) then it is not going to be considered a
# potential match.

# 'pre_estimate' = LUB estimate from first pass, if this is second pass.

def really_find_links(AB, subproblems, get_pre_estimate):
  # This sets the 'link' property of ... some ... records.

  i = 0
  for (key, (us, vs)) in subproblems.items():
    if i % 1000 == 0:
      log("# Subproblem %s %s %s %s %s" % (i, len(us), len(vs), blurb(us[0]), blurb(vs[0])))
    i += 1
    # classify as A1, A2, A3 (HETEROTYPIC, REVIEW, HOMOTYPIC)
    for u in us:
      for v in vs:
        # **** COMPUTE DISTANCE if 2nd pass ****
        dist = compute_distance(u, v, get_pre_estimate)
        score = compute_score(u, v, dist)
        if homotypic_score(score):
          consider_link(u, v)
          consider_link(v, u)
  #if get_link(AB.in_left(AB.A.top), None) != AB.in_right(AB.B.top):
  #  log("tops don't link")
  report_on_links(AB, subproblems)

# Compare pick_better_record in some_exemplar.py ...
# We have distances on the second pass, which can either remove or create ambiguities.
# So, ignore what was there previously.

def consider_link(u, v):
  set_link(u, pick_better_record(get_link(u, None), v))
  # equate_exemplars(u, v)

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

# -----------------------------------------------------------------------------

# Calibration  -- no longer used

# The most dissimilar things that are similar.  Unite records that are 
# this similar (minimally similar) or more so.
THRESH1 = compute_parts_score(parse_name("Foo bar Jones, 1927"),
                              parse_name("Foo bar"),
                              None)
THRESH2 = compute_parts_score(parse_name("Foo bar Jones, 1927"),
                              parse_name("Quux bar (Jones, 1927)"),
                              5)
THRESH = min(THRESH1, THRESH2)
log("# Match predicted if score >= %s = min(%s,%s)" %
    (THRESH, THRESH1, THRESH2))

THRESH_Q = compute_parts_score(parse_name("Foo bar Jones, 1927"),
                               parse_name("Quux bar Jones, 1927"),
                               5)
log("# For distinct genus starts: %s" % THRESH_Q)
assert THRESH_Q < THRESH, \
  (parse_name("Foo bar Jones, 1927"),
   parse_name("Quux bar Jones, 1927"))

# The most similar things that are dissimilar.  Distinguish records that 
# are this different (minimally different) or more so.
NOTHRESH = compute_parts_score(parse_name("Foo bar Jones, 1927"),
                               parse_name("Foo bar Smith, 1927"))
log("# Unmatch predicted if score <= %s" % NOTHRESH)


# Phase this out?  What's it for?  find_estimate

def get_mutual_link(u, default=-19): # misplaced
  v = get_link(u, None)
  if v:
    u_back = get_link(v, default)
    if u_back == u:
      return v
  return default

link_prop = prop.declare_property("link", inherit=False)
(get_link, set_link) = prop.get_set(link_prop)


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
