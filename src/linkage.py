#!/usr/bin/env python3

# Record linkage - Calculate match score / probability
# See wikipedia record linkage article

import util
import property as prop
import parser, rows
import simple
from checklist import *
from workspace import *
from rcc5 import EQ, NEQ
from parse import parse_name

HOMOTYPIC   = EQ           # Hmm, not sure about this
HETEROTYPIC = NEQ
REVIEW      = HOMOTYPIC | HETEROTYPIC

NEUTRAL = 0

# For each record, get list of matching records.... hmm
# Future: read exceptions or whole thing from a file

def find_links(AB, m_iter=None):
  assert not m_iter, 'NYI'
  # This sets the 'link' property of ... some ... records.
  subproblems = find_subproblems(AB)
  count = 0
  for (us, vs) in subproblems.values():
    # classify as A1, A2, A3 (HETEROTYPIC, REVIEW, HOMOTYPIC)
    for u in us:
      for v in vs:
        have1 = get_link(u, None)
        have2 = get_link(v, None)
        if have1 == None or have2 == None:
          score = compute_score(u, v)
          log("# for (%s, %s) have (%s, %s) score %s" %
              (blurb(u), blurb(v),
               blurb(have1) if have1 else '-',
               blurb(have2) if have2 else '-',
               score))
          count += improve_link(u, v, score)
          count += improve_link(v, u, score)
          #log("#  Improved: %s %s" % (blurb(u), blurb(get_link(u))))
  assert get_link(AB.in_left(AB.A.top)) == AB.in_right(AB.B.top)
  log("# %s links forged" % count)

def improve_link(u, v, score):
  log(explain(score))
  if score_to_ship(score) != HOMOTYPIC:
    return 0
  else:
    have_v = get_link(u, None) # (v0, score) Anything to improve on?
    if have_v:
      mrca = simple.mrca(have_v, v)
      # We want the one that is... more protonymical (fewer () )
      if have_v == mrca:
        # v < have_v.  prefer v ??
        set_link(u, v)
        return 0                # Replacing
      elif v == mrca:
        # have_v < v.  Leave it alone.
        return 0
      else:
        # ambiguous.  wait until later to figure this out.
        set_link(u, False)
        return -1               # Retracting
    # Adding
    set_link(u, v)
    return 1

def score_to_ship(score):
  if score >= THRESH:     return HOMOTYPIC
  elif score <= NOTHRESH: return HETEROTYPIC
  else: return REVIEW

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
    assert key
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
  assert ep != None
  return parts.genus if ep == '' else ep

# -----------------------------------------------------------------------------
# Score potentially contipic taxa.
# If distance is close, give a pass on genus mismatch.

def z_compute_score(u, v, nearness=None):
  return compute_score(u, v, nearness)

def compute_score(u, v, nearness=None):
  return compute_parts_score(get_parts(u),
                             get_parts(v),
                             nearness)

# Score potentially contipic names from 0 to 100

def compute_parts_score(p, q, nearness=None):
  # First, normalize moving genera
  g1 = p.genus; g2 = q.genus
  if (g1 != None and g2 != None and g1 != g2):
    if p.moved: g1 = None
    if q.moved: g2 = None

  hits = misses = 0
  if nearness != None:
    # -1 = too far apart to link
    # 0 = neutral, intermediate distance, no info
    # 1 = nearby, hit  ?
    if nearness > 0:
      hits += 1
    elif nearness < 0:
      misses += 1
  if p.year != None and q.year != None:
    if p.year == q.year: hits += 2
    else: misses += 2
  if p.token != None and q.token != None:
    if p.token == q.token: hits += 4
    else: misses += 4
  if g1 != None and g2 != None:
    if g1 == g2: hits += 8
    else: misses += 8
  if p.epithet != None and q.epithet != None:
    if p.epithet != q.epithet: misses += 16
    else: hits += 16
  if misses > 0: return NEUTRAL - misses
  else: return NEUTRAL + hits

def explain(score):
  def explode(things):
    return ', '.join((z
                      for (i, z)
                      in zip(range(0, 5),
                             (
                              "nearness", # 1  = 1 << 0
                              "year",     # 2
                              "token",    # 4
                              "genus",    # 8
                              "epithet",  # 16
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
  return "%s(%s)" % (word, explode(bits))
    

# Calibration

# The most dissimilar things that are similar.  Unite records that are 
# this similar (minimally similar) or more so.
THRESH = compute_parts_score(parse_name("Foo bar"),
                             parse_name("Foo bar"))
log("# Match predicted if score >= %s" % THRESH)

# The most similar things that are dissimilar.  Distinguish records that 
# are this different (minimally different) or more so.
NOTHRESH = compute_parts_score(parse_name("Foo bar Jones, 1927"),
                               parse_name("Foo bar Smith, 1927"))
log("# Unmatch predicted if score <= %s" % NOTHRESH)

THRESH2 = compute_parts_score(parse_name("Foo bar Jones, 1927"),
                              parse_name("Quux bar (Jones, 1927)"),
                              10)
log("# For synonym analysis: %s" % THRESH2)

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
  args=parser.parse_args()

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
