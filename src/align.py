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

  roots = set_congruences(a_roots, b_roots, sum, rm_sum)
  roots = build_tree(a_roots, b_roots, sum, rm_sum)

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
    note_match(x, y, row[remark_pos], sum)
  return sum

def note_match(x, y, remark, sum):
  assert remark
  (key_to_union, in_a, in_b) = sum
  if x and y:
    j = mep_get(in_a, x, None)
    assert mep_get(in_b, y, None) == j
  elif x:
    j = mep_get(in_a, x, None)
  elif y:
    j = mep_get(in_b, y, None)
  else:
    j = None
  if j:
    assert out_a(j) == x
    assert out_b(j) == y
    return j
  # Invent a key?? and store it in the sum ???
  # To avoid an A/B conflict we'd need the B key->usage map?
  if y:
    key = get_key(y)
    if key in key_to_union: key = "B.%s" % key
  else:
    key = get_key(x)
    if key in key_to_union: key = "A.%s" % key
  name = get_canonical(y) if y else get_canonical(x)
  u = make_union(key, x, y, remark, name or MISSING)
  if x: mep_set(in_a, x, u)
  if y: mep_set(in_b, y, u)
  mep_set(key_to_union, key, u)
  return u

def combine_remarks(*remarks):
  return ";".join([r for r in remarks if r != MISSING])

# -----------------------------------------------------------------------------
# Get parent/child and accepted/synonym relationships

children_prop = prop.Property("children")
synonyms_prop = prop.Property("synonyms")
get_children = prop.getter(children_prop)
get_synonyms = prop.getter(synonyms_prop)
set_children = prop.setter(children_prop)
set_synonyms = prop.setter(synonyms_prop)

def get_inferiors(p):
  return get_children(p, []) + get_synonyms(p, [])

def get_superior(x):
  return get_parent(x, None) or get_accepted(x, None)

def collect_children(items):
  roots = []
  for item in items:

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
    for c in get_inferiors(x):
      cache(c, n+1)
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
    p = get_superior(x)
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
  a_tipwards = half_find_tipward(a_roots, rm_in_a, out_b)
  b_tipwards = half_find_tipward(b_roots, rm_in_b, out_a)
  count = 0
  for u in b_tipwards.values():
    if mep_get(a_tipwards, u, None):
      set_xmrca(out_a(u), out_b(u))
      set_xmrca(out_b(u), out_a(u))
      count += 1
  if len(a_tipwards) != count:
    print("-- align: %s tipward record matches in A" % len(a_tipwards),
          file=sys.stderr)
  if len(b_tipwards) != count:
    print("-- align: %s tipward record matches in B" % len(b_tipwards),
          file=sys.stderr)
  print("-- align: %s tipward record matches" % count, file=sys.stderr)

# -----------------------------------------------------------------------------
# 7. how-related WITHIN hierarchy
#   checklist.how_related

EQ = 1
LT = 2
GT = 3
DISJOINT = 4
CONFLICT = 5

def how_related(x, y):
  (x1, y1) = find_peers(x, y)
  if x1 == y1:
    if x == y:
      return EQ
    elif x1 == x:
      return GT
    else:
      return LT
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
  o = get_xmrca(p)
  if o == None: return (None, None)
  rel = how_related(o, q)
  if rel == LT or rel == EQ:
    return (p, None)
  elif rel == DISJOINT:
    return (None, q)
  assert rel == GT
  p_and_q = p_not_q = None                # in A
  for x in get_inferiors(p):
    (x_and_q, x_not_q) = seek_conflict(x, q)
    if x_and_q and x_not_q:
      return (x_and_q, x_not_q)
    if x_and_q: p_and_q = x
    if x_not_q: p_not_q = x
    if p_and_q and p_not_q:             # hack: cut it short
      return (p_and_q, p_not_q)
  return (p_and_q, p_not_q)

# -----------------------------------------------------------------------------

# Detect and record A/B congruences.
# In the process of doing this, also set parent pointers in the sum
# within "clusters".  A cluster is a set of nodes, in both trees, 
# that all subtend the same set of tipward record matches.

# At the end:
# 1. congruences have been established as needed
# 2. every node in every cluster has an assigned node in the union (injection)
# 3. every assigned node for every node in every cluster - other than the
#    rootward one - has its parent pointer set to another node in the
#    cluster, forming a single linear chain.

"""
Diagram explaining the variables used in "weave".  Rootward is toward
the top.

  p  <= h =>  q   parents outside the cluster
  v     |     v
  x  ?  k  ?  y   rootward nodes in cluster
  |           |
  |           m   record match to s
  s  = (r) =  |
  |           t
  |  \     /  |
  |     k     |   s, t, and k are the iteration variables
  |     |     |
  |           |
  |  /     \  |
  u  =  g  =  v   tipward nodes in cluster
"""

def connector(rm_sum, sum):
  (_, rm_in_a, rm_in_b) = rm_sum
  def connect(x, y, remark):
    assert remark
    remarks = [remark]
    if x:
      j = mep_get(rm_in_a, x, None)
      assert j
      if j: remarks.append(get_remark(j))
    if y:
      j = mep_get(rm_in_b, y, None)
      assert j
      if j: remarks.append(get_remark(j))
    return note_match(x, y, combine_remarks(*remarks), sum)
  return connect

def set_congruences(a_roots, b_roots, sum, rm_sum):
  (_, in_a, in_b) = sum
  (_, rm_in_a, _) = rm_sum

  connect = connector(rm_sum, sum)

  def weave(x, u, v):
    # v and u need to be linked.
    g = connect(u, v, "join point")

    s = get_superior(u)
    t = get_superior(v)
    k = g
    while True:
      s_done = (not s or get_xmrca(s) != v)
      t_done = (not t or get_xmrca(t) != u)
      if s_done and t_done:
        break
      elif s_done:
        k = connect(None, t, "upper t")
        t = get_superior(t)
      elif t_done:
        k = connect(s, None, "upper s")
        s = get_superior(s)
      else:
        # Use record matching to connect cluster nodes when possible
        r = mep_get(rm_in_a, s)
        if r:
          m = out_b(r)
          if m == t:
            k = connect(s, m, "record match in cluster")
            s = get_superior(s)
            t = get_superior(t)
          elif m and get_xmrca(m) == u and get_level(m) < get_level(t):
            k = connect(None, t, "t inferior to record match")
            t = get_superior(t)
          else:
            k = connect(s, None, "near miss")
            s = get_superior(s)
        else:    # no record match
          k = connect(s, None, "no record matched in cluster")
          s = get_superior(s)
      set_parent(g, k)
      g = k
    return (k, s, t)

  def traverse(x):
    p = get_superior(x)
    v = get_xmrca(x, None)
    if v != None:
      u = get_xmrca(v)          # in A
      if how_related(x, u) != LT:
        (k, s, t) = weave(x, u, v)
        # k is the rootward cluster node in sum
    # Deal with descendants of x, and inject x
    for c in get_inferiors(x): traverse(c)

  for x in a_roots: traverse(x)

# -----------------------------------------------------------------------------
# Now set the parent pointers for the merged tree

def build_tree(a_roots, b_roots, sum, rm_sum):
  (_, in_a, in_b) = sum
  connect = connector(rm_sum, sum)
  def fasten_a(x):
    return connect(x, None, "peripheral in A")
  def fasten_b(y):
    return connect(None, y, "peripheral in B")
  finish(b_roots,
         in_b, in_a, fasten_b, fasten_a, True)
  finish(a_roots,
         in_a, in_b, fasten_a, fasten_b, False)
  unions = all_unions(sum)
  print("align: %s nodes in sum" % len(unions), file=sys.stderr)

  roots = collect_children(unions)
  print("align: %s roots in sum" % len(roots), file=sys.stderr)
  return roots

def cap(x, in_a, fasten):
  j = mep_get(in_a, x, None)
  if not j:
    j = fasten(x)
    assert mep_get(in_a, x) == j
  return j

def finish(roots, in_a, in_b, fasten_a, fasten_b, priority):
  def traverse(x):
    j = cap(x, in_a, fasten_a)
    for c in get_inferiors(x): traverse(c)
    if not get_superior(j):
      (pq, ab) = determine_superior_in_sum(x, in_a, in_b, priority)
      if pq:
        if ab:
          h = cap(pq, in_a, fasten_a)
        else:
          h = cap(pq, in_b, fasten_b)
        set_superior(j, h)
  for root in roots: traverse(root)

def determine_superior_in_sum(x, in_a, in_b, priority):
  v = get_xmrca(x, None)
  p = get_superior(x)
  if v == None:
    return (p, True)
  u = get_xmrca(v)
  q = v
  # Ascend until we find a q that's definitely bigger than x
  if how_related(x, u) != LT:
    # In cluster
    while True:
      q = get_superior(q)
      if not q or get_xmrca(q) != u:
        break
  def relax(q):
    # No monotype chains.  Ordering is by tipward-match sets.
    if p == None and q == None:
      return (None, None)
    else:
      if p == None:
        return (q, False)
      elif q == None:
        return (p, True)
      else:
        (rel, _, _) = related_how(p, q)
        if rel == GT:
          return (q, False)
        elif rel == LT or rel == EQ:
          return (p, True)
        elif rel == CONFLICT:
          if priority:          # is this is the priority tree,
            return (p, True)
          else:
            return relax(get_superior(q))
        else:
          assert False
  return relax(q)

# j is to be either a child or synonym of k

def set_superior(j, k):
  assert k
  assert get_superior(j) == None
  x = out_a(j)
  y = out_b(j)

  # If j has any children or synonyms, then j is a child, not a synonym
  if ((x and len(get_inferiors(x)) > 0) or
      (y and len(get_inferiors(y)) > 0)):
    set_parent(j, k)

  # If j is considered a synonym on both sides, then j is a synonym
  elif y and get_accepted(y, None):
    set_accepted(j, k)

  elif not y and get_accepted(x, None):
    set_accepted(j, k)

  else: set_parent(j, k)

# -----------------------------------------------------------------------------
# 14. emit the new sum with additional parent column

# Returns a row generator

def generate_sum(sum):
  yield ["taxonID", "taxonID_A", "taxonID_B",
         "parentNameUsageID", "acceptedNameUsageID",
         "canonicalName", "remark"]

  for union in all_unions(sum):
    a_usage = out_a(union, None)
    b_usage = out_b(union, None)
    p = get_parent(union, None)
    a = get_accepted(union, None)
    yield [get_key(union),
           get_key(a_usage) if a_usage else MISSING,
           get_key(b_usage) if b_usage else MISSING,
           get_key(p) if p else MISSING,
           get_key(a) if a else MISSING,
           get_canonical(union, MISSING),
           get_remark(union)]

def all_unions(sum):
  (key_to_union, _, _) = sum
  return key_to_union.values()

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
