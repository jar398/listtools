#!/usr/bin/env python3

# Align to usage lists, with sensitivity to parent/child relations

import sys, csv
from util import windex, MISSING
import property as prop

# -----------------------------------------------------------------------------
# Supervise the overall process

def align(a_iterator, b_iterator, rm_sum_iterator):
  a_usage_dict = load_usages(a_iterator)
  b_usage_dict = load_usages(b_iterator)

  # Collect children so we can say children[x]
  a_roots = collect_children(a_usage_dict.values())
  b_roots = collect_children(b_usage_dict.values())

  # Prepare for doing within-tree MRCA operations
  cache_levels(a_roots)
  cache_levels(b_roots)

  # Read record matches
  rm_sum = get_sum(rm_sum_iterator, a_usage_dict, b_usage_dict)

  sum = ({}, {}, {})
  find_tipward_record_matches(a_roots, b_roots, rm_sum, sum)
  cache_xmrcas(a_roots, b_roots, sum)

  # annotate_unmatched(a_roots, b_roots, sum)
  build_trees(a_roots, b_roots, sum, rm_sum)
  check_trees(a_roots, b_roots, sum)

  return generate_sum(sum)

# -----------------------------------------------------------------------------
# Hierarchy file ingest

key_prop = prop.Property("primary_key")
get_key = prop.getter(key_prop)
canonical_prop = prop.Property("canonical")
get_canonical = prop.getter(canonical_prop)
make_usage = prop.constructor(key_prop, canonical_prop)

parent_prop = prop.Property("parent")
accepted_prop = prop.Property("accepted")
set_parent = prop.setter(parent_prop)
set_accepted = prop.setter(accepted_prop)
get_parent = prop.getter(parent_prop)
get_accepted = prop.getter(accepted_prop)

def load_usages(iterator):
  header = next(iterator)
  key_pos = windex(header, "taxonID")
  parent_pos = windex(header, "parentNameUsageID")
  accepted_pos = windex(header, "acceptedNameUsageID")
  canonical_pos = windex(header, "canonicalName")

  temp = []
  key_to_usage = {}
  for row in iterator:
    key = row[key_pos]
    usage = make_usage(key, row[canonical_pos])
    key_to_usage[key] = usage
    temp.append((usage,
                 row[parent_pos],
                 row[accepted_pos] if accepted_pos else MISSING))
  for (usage, parent_key, accepted_key) in temp:
    probe = key_to_usage.get(parent_key)
    if probe: set_parent(usage, probe)
    probe = key_to_usage.get(accepted_key)
    if probe: set_accepted(usage, probe)
  return key_to_usage

# -----------------------------------------------------------------------------
# Sum file ingest

# Record-match file ingest

out_a_prop = prop.Property("out_a")
out_a = prop.getter(out_a_prop)
out_b_prop = prop.Property("out_b")
out_b = prop.getter(out_b_prop)
remark_prop = prop.Property("remark")
get_remark = prop.getter(remark_prop)

make_union = prop.constructor(key_prop, out_a_prop, out_b_prop,
                              remark_prop, canonical_prop)

nodefault = []
def mep_get(mep, x, default=nodefault):
  if default is nodefault:
    return mep[prop.get_identity(x)]
  else:
    return mep.get(prop.get_identity(x), default)
def mep_set(mep, x, j):
  mep[prop.get_identity(x)] = j


def get_sum(iterator, a_usage_dict, b_usage_dict):
  header = next(iterator)
  key_pos = windex(header, "taxonID")
  usage_a_pos = windex(header, "taxonID_A")
  usage_b_pos = windex(header, "taxonID_B")
  remark_pos = windex(header, "remark")

  sum = ({}, {}, {})
  for row in iterator:
    key = row[key_pos]
    x = a_usage_dict.get(row[usage_a_pos])
    y = b_usage_dict.get(row[usage_b_pos])
    connect(x, y, row[remark_pos], sum)
  return sum

def connect(x, y, remark, sum):
  (key_to_union, in_a, in_b) = sum
  assert not x or mep_get(in_a, x, None) == None
  assert not y or mep_get(in_b, y, None) == None
  # Invent a key?? and store it in the sum ???
  # To avoid an A/B conflict we'd need the B key->usage map?
  key = get_key(y) if y else get_key(x)
  name = get_canonical(y) if y else get_canonical(x)
  assert name
  u = make_union(key, x, y, remark, name)
  if x: mep_set(in_a, x, u)
  if y: mep_set(in_b, y, u)
  mep_set(key_to_union, key, u)
  return u

def combine_remarks(*remarks):
  "|".join([r for r in remarks if r != MISSING])

# -----------------------------------------------------------------------------
# Get parent/child and accepted/synonym relationships

children_prop = prop.Property("children")
synonyms_prop = prop.Property("synonyms")
get_children = prop.getter(children_prop)
get_synonyms = prop.getter(synonyms_prop)
set_children = prop.setter(children_prop)
set_synonyms = prop.setter(synonyms_prop)

def get_inferiors(y):
  return get_children(y, []) + get_synonyms(y, [])

def get_superior(x):
  return get_parent(x, None) or get_accepted(x, None)

def collect_children(items):
  roots = []
  for item in items:

    assert get_canonical(item, None)

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

    # Add item to list of accepted's synonyms
    accepted_item = get_accepted(item, None)
    if accepted_item:
      ch = get_synonyms(accepted_item, None)
      if ch:
        ch.append(item)
      else:
        set_synonyms(accepted_item, [item])

  return roots

# -----------------------------------------------------------------------------
# 6. mrca

# Cache every node's level (distance to root)
#   simple recursive descent from roots

level_prop = prop.Property("level")
get_level = prop.getter(level_prop)
set_level = prop.setter(level_prop)
levels = {}                     # global

def cache_levels(roots):
  def cache(x, n):
    set_level(x, n)
    for y in get_inferiors(x):
      cache(y, n+1)
  for root in roots:
    cache(root, 1)

def find_peers(x, y):
  if get_level(x) == get_level(y):
    return (x, y)
  elif get_level(x) < get_level(y):
    return find_peers(x, get_superior(y))
  else:
    return find_peers(get_superior(x), y)

def mrca(x, y):
  (x, y) = find_peers(x, y)
  if x == y:
    return x
  else:
    p = get_superior(x, None)
    if p == None:
      return None               # Foo
    else:
      return mrca(p, get_superior(y))

# -----------------------------------------------------------------------------
# Analyze tipward record matches.

# A set of record matches is a 'tipward record match set' if each
# match (x, y) in the set has the property that no ancestor of either
# x or y is also matched in the set.

# The tipward record matches are stored in the given sum object.

def find_tipward_record_matches(a_roots, b_roots, rm_sum, sum):

  def half_find_tipward(roots, inject, outject):
    tipwards = {}               # mep
    def traverse(x_usage):
      saw = False
      for child in get_inferiors(x_usage):
        saw = traverse(child) or saw
      if not saw:
        u = mep_get(inject, x_usage)
        y_usage = outject(u)
        if y_usage:
          mep_set(tipwards, u, u)
          return True
      return saw
    for root in roots:
      traverse(root)
    return tipwards

  (_, rm_in_a, rm_in_b) = rm_sum
  a_tipwards = half_find_tipward(a_roots, rm_in_a, out_a)
  b_tipwards = half_find_tipward(b_roots, rm_in_b, out_b)
  count = 0
  for (u_identity, u) in b_tipwards.items():
    if u_identity in a_tipwards:
      set_xmrca(out_a(u), out_b(u))
      set_xmrca(out_b(u), out_a(u))
      count += 1
  print("-- %s tipward record matches" % count, file=sys.stderr)

# -----------------------------------------------------------------------------
# 7. how-related WITHIN hierarchy
#   checklist.how_related

EQ = 1
LT = 2
GT = 3
DISJOINT = 4
CONFLICT = 5
GT_OR_CONFLICT = 6

def how_related(x, y):
  (x1, y1) = find_peers(x, y)
  if x1 == None:
    return DISJOINT
  elif x1 == y1:
    return EQ
  elif x1 == x:
    return LT
  elif y1 == y:
    return GT
  else:
    return DISJOINT

# -----------------------------------------------------------------------------
# 8. cross-mrcas
#   simplified alignment.infer_partners

xmrca_prop = prop.Property("xmrca")
get_xmrca = prop.getter(xmrca_prop)
set_xmrca = prop.setter(xmrca_prop)

def half_xmrcas(x):
  m = get_xmrca(x, None)
  if m:
    return m                    # Tipward record match
  for c in get_inferiors(x):
    n = half_xmrcas(c)
    if n:
      if m:
        m = mrca(m, n)
        # If m is None here, then we need to merge roots!
      else:
        m = n
  if m:
    set_xmrca(x, m)
  return m

def cache_xmrcas(a_roots, b_roots, sum):
  (_, in_a, in_b) = sum
  for root in a_roots:
    half_xmrcas(root)
  for root in b_roots:
    half_xmrcas(root)

# -----------------------------------------------------------------------------

# Compare two parent candidates, either x's parent p in A,
# vs. its cross-mrca q in B.  Pick the smaller one to be the parent in
# the sum.
# Return values are (rcc5, e, f) where:
#   e is < both p and q
#   f is < p but ! q

# By doing this twice, once in each direction, it's possible to
# distinguish GT from CONFLICT.

def related_how(p, q):
  c = get_xmrca(p)
  d = get_xmrca(q)
  if c == None or d == None:
    return (DISJOINT, None, None)
  elif get_xmrca(c) == d:
    if p == d and q == c:
      return (EQ, None, None)
    else:
      cmp = ((get_level(d) - get_level(p)) -
             (get_level(c) - get_level(q)))
      if cmp < 0:
        return (LT, None, None)
      elif cmp == 0:
        # TBD: Show more finesse
        # especially: make use of record matches between monotype chains
        return (EQ, None, None)
      else:
        return (GT, None, None)
  else:
    (e, f) = seek_conflict(p, q)
    assert e
    if f:
      (e2, f2) = seek_conflict(q, p)
      assert e2
      if f2:
        return (CONFLICT, e, f)    # or (e, f, f2)?
      else:
        return (GT, e, f)
    else:
      return (LT, e, None)

# p < xmrca(q); but is p < q?

def seek_conflict(p, q):
  c = get_xmrca(p)
  if c == None: return (None, None)
  rel = how_related(c, q)
  if rel == LT or rel == EQ:
    return (p, None)
  elif rel == DISJOINT:
    return (None, q)
  p_and_q = p_not_q = None                # in A
  for ch in get_inferiors(p):
    (ch_and_q, ch_not_q) = seek_conflict(ch, q)
    if ch_and_q and ch_not_q:
      return (ch_and_q, ch_not_q)
    if ch_and_q: p_and_q = ch
    if ch_not_q: p_not_q = ch
    if p_and_q and p_not_q:             # hack: cut it short
      return (p_and_q, p_not_q)
  return (p_and_q, p_not_q)

# -----------------------------------------------------------------------------

conflict_prop = prop.Property("conflict")
get_conflict = prop.getter(conflict_prop)
set_conflict = prop.setter(conflict_prop)

# Determine whether x in A is consistent with the B hierarchy.
# Let y = xmrca(x).
# x is consistent with B if it's consistent with every child of y.
# Do we know that x > c for all children c of y??
# Yes, because if x <= c, then c would be the xmrca, not y.

# This is merge.test_compatibility  in cldiff, sort of

def annotate_unmatched(a_roots, b_roots, sum):
  def traverse(x):
    (_, in_a, _) = sum
    if prop.get_identity(x) in in_a:
      # Already congruent, or something
      pass
    else:
      # Check for incompatibility of x with B
      y = get_xmrca(x)
      if y == None:
        connect(x, None, ".add peripherally", sum)
      else:
        for d in get_inferiors(y):
          (rel, e, f) = related_how(d, x)
          if rel == CONFLICT:
            set_conflict(x, (e, f))
            connect(x, None, "conflict: %s | %s" % (e, f), sum)
          else:
            connect(x, None, ".add to hierarchy", sum)

    # Now do the same for all descendants
    for c in get_inferiors(x):
      traverse(c)

  for x in a_roots: traverse(x)

  def traverse_b(y):
    (_, _, in_b) = sum
    if prop.get_identity(y) in in_b:
      pass
    else:
      connect(None, y, MISSING, sum)

  for y in b_roots: traverse_b(y)

# -----------------------------------------------------------------------------

# 13. set parent pointers
#   - build up new A/B match set by recursive descent
#   - retract nullified terminal nodes
#   - retract matched internal nodes (keep only matched)
#   - add new internal matches (singletons)

"""
  p  <= h =>  q   parents outside the chains
  v     |     v
  x     k     y   rootward nodes in chains
  |           |
  |           m   record match to s
  s  =  n  =  |
  |           t
  |  \     /  |
  |     k     |
  |     |     |
  |     j     |
  |  /     \  |
  u  =  g  =  v   tipward nodes in chains
"""

def build_trees(a_roots, b_roots, sum, rm_sum):
  (_, in_a, in_b) = sum
  (_, rm_in_a, rm_in_b) = rm_sum

  def weave(x, u, v):
    # v and u need to be linked.
    g = connect(u, v, "join point", sum)

    s = get_superior(u)
    t = get_superior(v)
    k = g
    while True:
      s_done = (not s or get_xmrca(s) != v)
      t_done = (not t or get_xmrca(t) != u)
      if s_done and t_done:
        break
      elif s_done:
        k = connect(None, t, "upper t", sum)
        t = get_superior(t)
      elif t_done:
        k = connect(s, None, "upper s", sum)
        s = get_superior(s)
      else:
        n = mep_get(rm_in_a, s)
        if n and out_b(n) == t:
          k = connect(s, t, "record match in cluster", sum)
          s = get_superior(s)
          t = get_superior(t)
        elif ((lambda m: m and get_xmrca(m) == u and get_level(m) <= get_level(t)) \
              (out_b(n))):
          k = connect(None, t, "t unmatched in cluster", sum)
          t = get_superior(t)
        else:    # no record match
          k = connect(s, None, "s unmatched in cluster", sum)
          s = get_superior(s)
    return (k, s, t)

  # Need to do this bottom up.

  def traverse(x):
    p = get_parent(x, None)
    v = get_xmrca(x)
    if v == None:
      k = connect(x, None, "unmatched", sum)
      set_superior(k, get_mep(in_a, p))
    else:
      u = get_xmrca(v)          # in A
      if how_related(x, u) == LT:
        k = connect(x, None, "unmatched", sum)
        choose_parent(k, p, v, sum)
      else:
        (k, s, t) = weave(x, u, v)
        assert s == p
        assert not t or mrca(t, v) == t
        choose_parent(k, s, t, sum)
        x = u
    # Deal with descendants of x, and inject x
    for c in get_inferiors(x): traverse(c)
  for x in a_roots: traverse(x)

# Returns either in_a(p) or in_b(q)
# We are given
#  (a) p in A and q in B such that p is the least node in A
#   greater than x and q is the least node in B greater than x
# We want to choose
#  (b) h = the smaller of inj_a(p) and inj_b(q) in A+B

def choose_parent(k, p, q, sum):
  (_, in_a, in_b) = sum
  # No monotype chains.  Ordering is by tipward-match sets.
  if p == None:
    h = mep_get(in_b, q) if q else None
  elif q == None:
    h = mep_get(in_a, p)
  else:
    (rel, _, _) = related_how(p, q)
    if rel == LT or rel == EQ:
      h = mep_get(in_a, p)
    elif rel == GT:
      h = mep_get(in_b, q)
    else:  # rel == CONFLICT
      h = choose_parent(p, get_parent(q), sum)
  set_superior(k, h)

# j is to be either a child or synonym of k

def set_superior(j, k):
  ja = out_a(j)
  jb = out_b(j)

  # If j has any children, then j is a child, not a synonym
  if ((jb and get_children(jb, None)) or
      (ja and get_children(ja, None))):
    set_parent(j, k)

  # If j is considered a synonym on both sides, then j is a synonym
  elif ((not ja or get_accepted(ja, None)) and
        (not jb or get_accepted(jb, None))):
    set_accepted(j, k)

  else: set_parent(j, k)

def check_trees(a_roots, b_roots, sum):
  (_, in_a, in_b) = sum
  def traverse(x, in_a, out_a):
    j = mep_get(in_a, x, None)
    assert j
    assert out_a(j) == x
    for c in get_inferiors(x): traverse(c)
  for root in a_roots: traverse(root, in_a, out_a)
  for root in b_roots: traverse(root, in_b, out_b)

  unions = {}
  for union in in_a.values():
    assert out_a(union, None)
    mep_set(unions, union, union)
  for union in in_b.values():
    assert out_b(union, None)
    mep_set(unions, union, union)
  roots = collect_children(unions.values())
  print("%s roots in sum" % len(roots))


# -----------------------------------------------------------------------------
# 14. emit the new sum with additional parent column

# Returns a row generator

def generate_sum(sum):
  (key_to_union, _, _) = sum

  yield ["taxonID", "taxonID_A", "taxonID_B",
         "parentNameUsageID", "acceptedNameUsageID",
         "canonicalName", "remark"]

  for (key, union) in key_to_union.items():
    a_usage = out_a(union, None)
    b_usage = out_b(union, None)
    name = (get_canonical(union, None) or get_canonical(b_usage, None)
            or get_canonical(a_usage, MISSING))
    p = get_parent(union, None)
    a = get_accepted(union, None)
    yield [key,
           get_key(a_usage) if a_usage else MISSING,
           get_key(b_usage) if b_usage else MISSING,
           get_key(p) if p else MISSING,
           get_key(a) if a else MISSING,
           name, get_remark(union)]

def write_generated(gen, outfile):
  writer = csv.writer(outfile)
  for row in gen:
    writer.writerow(row)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    TBD
    """)
  parser.add_argument('--target', help="B hierarchy")
  parser.add_argument('--matches', help="record matches")
  args=parser.parse_args()

  a_file = sys.stdin
  b_path = args.target
  rm_sum_path = args.matches

  with open(b_path) as b_file:
    with open(sum_path) as sum_file:
      sum = align(csv.reader(a_file),
                  csv.reader(b_file),
                  csv.reader(rm_sum_file))
      write_generated(sum)
