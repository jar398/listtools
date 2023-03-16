#!/usr/bin/env python3

# Record linkage - Calculate match score / probability
# See wikipedia record linkage article

import util
import property as prop
import parser
from parse import parse_name
from checklist import *

NEUTRAL = 0

get_link = get_match
set_link = set_match

# For each record, get list of matching records.... hmm

def find_links(AB):
  # This sets the 'link' property of ... all? ... nodes.
  subproblems = find_subproblems(AB.A, AB.B)
  count = 0
  for (xs, ys) in subproblems.values():
    # classify as A1, A2, A3 (NEQ, NOINFO, EQ)
    for x in xs:
      for y in ys:
        have1 = get_link(x, False)
        have2 = get_link(y, False)
        if have1 == None or have2 == None:
          score = compute_score(x, y)
          count += improve_link(x, y, score)
          count += improve_link(y, x, score)
  return count

def improve_link(x, y, score):
  log(explain(score))
  ship = score_to_ship(score)
  if ship != EQ:
    return 0
  else:
    have_rel = get_link(x, None)
    if have_rel and have_rel.relationship == EQ:
      have_y = have_rel.record
      mrca = simple.mrca(have_y, y)
      if have_y == mrca:
        # y < have_y.  prefer it. ???
        set_link(x, relation(ship, y, "more specific"))
        return 0
      elif y == mrca:
        # have_y < y.  Leave it alone.
        return 0
      else:
        # inconsistent.  wait until later to figure this out.
        set_link(x, False)
        return -1
    elif have_rel != None:
      return 0
    set_link(x, relation(ship, y, "score %s" % score))
    return 1

def score_to_ship(score):
  if score >= THRESH:     return EQ
  elif score <= NOTHRESH: return NEQ
  else: return NOINFO


# Find blocks/chunks, one per epithet

def find_subproblems(A, B):
  (A_index, B_index) = \
    map(lambda C: \
        index_by_some_key(C,
                          # should use genus if epithet is missing
                          get_subproblem_key),
        (A, B))
  subprobs = {}
  for (val, xs) in A_index.items():
    ys = B_index.get(val, None)
    if ys != None:
      subprobs[val] = (xs, ys)
      # for x in xs: set_subproblem(x, subprob)
      # for y in ys: set_subproblem(y, subprob)
  return subprobs

# Each subproblem covers a single associated partial name

def get_subproblem_key(x):
  parts = get_parts(x)
  ep = parts.epithet
  assert ep != None
  return parts.genus if ep == '' else ep

# Returns dict value -> key

def index_by_some_key(A, fn):
  index = {}
  for x in postorder_records(A):
    key = fn(x)
    assert key
    have = index.get(key, None)
    if have:
      have.append(x)
    else:
      index[key] = [x]
  return index

# -----------------------------------------------------------------------------
# Score potentially contipic taxa.
# If distance is close, give a pass on genus mismatch.

def compute_score(x, y, nearness=None):
  return compute_parts_score(get_parts(x),
                             get_parts(y),
                             nearness)

def get_parts(x):
  return parse_name(get_scientific(x, None) or get_canonical(x))

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
                      in zip(range(0, 4),
                                  ("year", "token", "genus", "nearness", "epithet"))
                      if things & 1<<i != 0))
  if score == NEUTRAL:
    return "neutral"
  if score > NEUTRAL:
    return "same(%s)" % explode(score - NEUTRAL)
  else:
    return "different(%s)" % explode(NEUTRAL - score)
    
  

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


import argparse
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("name1")
  parser.add_argument("name2")
  args=parser.parse_args()
  log("Comparing '%s' to '%s'" % (args.name1, args.name2))
  s = compute_parts_score(parse_name(args.name1), parse_name(args.name2))
  log("Score = %s, relationship is %s" % (s, rcc5_symbol(score_to_ship(s))))
