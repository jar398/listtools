#!/usr/bin/env python3

# Record linkage - Calculate match score / probability
# See wikipedia record linkage article

import util
import property as prop
import parse, rows
import simple
from checklist import *
from workspace import *
from rcc5 import EQ, NEQ
from parse import parse_name

HOMOTYPIC   = EQ           # Hmm, not sure about this
HETEROTYPIC = NEQ
REVIEW      = HOMOTYPIC | HETEROTYPIC

NEUTRAL = 0

def find_links(AB, m_iter=None, nearness=None):
  return 

# For each record, get list of matching records.... hmm
# Future: read exceptions or whole thing from a file

def really_find_links(AB, distance):
  # This sets the 'link' property of ... some ... records.
  subproblems = find_subproblems(AB)

  for (key, (us, vs)) in subproblems.items():
    # classify as A1, A2, A3 (HETEROTYPIC, REVIEW, HOMOTYPIC)
    for u in us:
      for v in vs:
        have1 = get_link(u, None)
        have2 = get_link(v, None)
        if have1 == None or have2 == None:
          # **** COMPUTE DISTANCE if 2nd pass ****
          dist = distance(u, v)
          score = compute_score(u, v, dist)
          improve_link(u, v, score, key)
          improve_link(v, u, score, key)
          #log("#  Improved: %s %s" % (blurb(u), blurb(get_link(u))))
  assert get_link(AB.in_left(AB.A.top)) == AB.in_right(AB.B.top)
  report_on_links(AB, subproblems)

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
        if link2 == u:
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
        if link2 == v:
          pass                  # already counted
        elif link2 == None:
          half_count += 1
  log("-- Links: %s mutual, %s one-sided, %s ambiguous" %
      (full_count, half_count, amb_count))

def improve_link(u, v, score, key):
  if score_to_ship(score) != HOMOTYPIC:
    return 0
  else:
    have_v = get_link(u, None) # (v0, score) Anything to improve on?
    if have_v == None:                 # not None or False
      set_link(u, v)                   # One side, not nec. reciprocal
      #log("# Set link %s ~ %s" % (blurb(u), blurb(v)))
    elif have_v != False:       # ambiguous
      # ambiguous.  wait until later to figure this out.
      if monitor(u) or monitor(v):
        log("# linkage: Ambiguous (%s): %s -> %s & %s" %
            (key, blurb(u), blurb(have_v), blurb(v)))
      set_link(u, False)

def score_to_ship(score):
  if score >= THRESH:     return HOMOTYPIC    # EQ
  elif score <= NOTHRESH: return HETEROTYPIC  # NEQ
  else: return REVIEW                         # NOINFO

# Assumes both are descended from the same species (or genus?)

def homotypic(u, species):
  score = compute_score(u, species, 10)
  if score >= THRESH2: return HOMOTYPIC
  else: return REVIEW

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

# Each subproblem covers a single associated partial name

def get_subproblem_key(z):
  parts = get_parts(z)
  ep = parts.epithet
  key = ep if ep else parts.genus
  # assert key, parts    MSW has scientificName = '?'
  return key

# -----------------------------------------------------------------------------
# Score potentially contipic taxa.
# If distance is close, give a pass on genus mismatch.

def compute_score(u, v, distance=None):
  score = compute_parts_score(get_parts(u),
                              get_parts(v),
                              distance)
  if (monitor(u) or monitor(v)) and score >= NOTHRESH:
    log("# Score (%s, %s) = %s" % (blurb(u), blurb(v), score))
    #log("# %s" % (get_parts(u),))
    #log("# %s" % (get_parts(v),))
  return score

# Score potentially contipic names from 0 to 100

def compute_parts_score(p, q, distance=None):
  hits = misses = 0

  if p.year != None and q.year != None:
    if p.year == q.year: hits += 2
    else: misses += 2
  if p.token != None and q.token != None:
    if p.token == q.token: hits += 4
    else: misses += 4

  if p.genus != None and q.genus != None:
    #log("# 1 comparing %s, %s (%s, %s)" % (p.genus, q.genus, p.protonymp, q.protonymp))
    if p.genus == q.genus:
      #log("# 2 comparing %s, %s" % (p.genus, q.genus))
      hits += 8
    elif p.protonymp == True and q.protonymp == True:
      #log("# 4 comparing %s, %s -> miss" % (p.genus, q.genus))
      misses += 8
    else:
      #log("# 3 comparing %s, %s -> ?" % (p.genus, q.genus))
      pass

  if distance != None:
    # 15 = too far apart to link
    # 9 = neutral, most distant linkable
    # 0 = equated, exact hit
    if distance <= 9: hits += 16
    elif distance > 15: misses += 16

  if p.epithet != None and q.epithet != None:
    if p.epithet == q.epithet: hits += 32
    else: misses += 32
  if misses > 0: return NEUTRAL - misses
  else: return NEUTRAL + hits

def explain(score):
  def explode(things):
    return ', '.join((z
                      for (i, z)
                      in zip(range(0, 5),
                             (
                              "unused", # 1  = 1 << 0
                              "year",     # 2
                              "token",    # 4
                              "genus",    # 8
                              "vicinity", # 16
                              "epithet",  # 32
                                         ))
                      if (things & (1<<i)) != 0))
  if score == NEUTRAL:
    word = "neutral"
    bits = 0
  if score > NEUTRAL:
    word = "link" if score >= THRESH else "similar, review"
    bits = score - NEUTRAL
  else:
    word = "nolink" if score < NOTHRESH else "dissimilar, review"
    bits = NEUTRAL - score
  return "%s %s(%s)" % (word, explode(bits))
    

# Calibration

# The most dissimilar things that are similar.  Unite records that are 
# this similar (minimally similar) or more so.
THRESH1 = compute_parts_score(parse_name("Foo bar Jones, 1927"),
                              parse_name("Foo bar"),
                              None)
THRESH2 = compute_parts_score(parse_name("Foo bar Jones, 1927"),
                              parse_name("Quux bar (Jones, 1927)"),
                              3)
THRESH = min(THRESH1, THRESH2)
log("# Match predicted if score >= %s = min(%s,%s)" %
    (THRESH, THRESH1, THRESH2))

THRESH_Q = compute_parts_score(parse_name("Foo bar Jones, 1927"),
                               parse_name("Quux bar Jones, 1927"),
                               3)
log("# For distinct genus starts: %s" % THRESH_Q)

# The most similar things that are dissimilar.  Distinguish records that 
# are this different (minimally different) or more so.
NOTHRESH = compute_parts_score(parse_name("Foo bar Jones, 1927"),
                               parse_name("Foo bar Smith, 1927"))
log("# Unmatch predicted if score <= %s" % NOTHRESH)

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
