#!/usr/bin/env python3

# Align two trees (given as DwC files i.e. usage lists), with
# sensitivity to parent/child relations.

import sys, csv, argparse
import util
from util import windex, MISSING
import property as prop
from property import mep_get, mep_set

troublemaker = "0000000000"

# -----------------------------------------------------------------------------
# Supervise the overall process

def align(a_iterator, b_iterator, rm_sum_iterator):
  (a_usage_dict, a_roots) = load_usages(a_iterator)
  (b_usage_dict, b_roots) = load_usages(b_iterator)

  # Read record matches
  rm_sum = load_sum(rm_sum_iterator, a_usage_dict, b_usage_dict)
  find_tipward_record_matches(a_roots, b_roots, rm_sum)
  cache_xmrcas(a_roots, b_roots)

  # Merge the two trees
  sum = ({}, {}, {})
  set_congruences(a_roots, b_roots, sum, rm_sum)
  roots = build_tree(a_roots, b_roots, sum, rm_sum)

  # Prepare to assign names
  assign_canonicals(sum)

  # Emit tabular version of merged tree
  return generate_sum(sum, roots, rm_sum)

def assign_canonicals(sum):
  (key_to_union, in_a, in_b) = sum
  b_index_by_name = {}
  for u in key_to_union.values():
    y = out_b(u)
    if y:
      name = get_canonical(y, None)
      if name:
        b_index_by_name[name] = y
  print("# %s B canonicals" % len(b_index_by_name),
        file=sys.stderr)
  count = 0
  losers = 0
  for u in key_to_union.values():
    y = out_b(u)
    if y:
      name = get_canonical(y, None)
    else:
      x = out_a(u)
      name = get_canonical(x, None)
      if name in b_index_by_name:
        y = b_index_by_name[name]
        name = name + " sec. A"
        if False:
          (rcc5, _, _, _) = related_how(x, y)
          print("  %s %s %s" %
                (name,
                 rcc5_symbols[rcc5],
                 get_canonical(y, "[no canonical]")),
                file=sys.stderr)
    if name:
      set_canonical(u, name)
      count += 1
    else:
      losers += 1
  print("# %s A+B canonicals (%s without)" % (count, losers),
        file=sys.stderr)

# -----------------------------------------------------------------------------
# Hierarchy file ingest

key_prop = prop.Property("primary_key")
get_key = prop.getter(key_prop)

canonical_prop = prop.Property("canonical")
get_canonical = prop.getter(canonical_prop)
set_canonical = prop.setter(canonical_prop)

scientific_prop = prop.Property("scientific")
get_scientific = prop.getter(scientific_prop)
set_scientific = prop.setter(scientific_prop)

rank_prop = prop.Property("rank")
get_rank = prop.getter(rank_prop)
set_rank = prop.setter(rank_prop)

make_usage = prop.constructor(key_prop)

parent_prop = prop.Property("parent")
get_parent = prop.getter(parent_prop)
set_parent = prop.setter(parent_prop)

accepted_prop = prop.Property("accepted")
get_accepted = prop.getter(accepted_prop)
set_accepted = prop.setter(accepted_prop)

def load_usages(iterator):
  header = next(iterator)
  key_pos = windex(header, "taxonID")
  canonical_pos = windex(header, "canonicalName")
  scientific_pos = windex(header, "scientificName")
  rank_pos = windex(header, "taxonRank")
  parent_pos = windex(header, "parentNameUsageID")
  accepted_pos = windex(header, "acceptedNameUsageID")

  fixup = []
  key_to_usage = {}
  for row in iterator:
    key = row[key_pos]
    usage = make_usage(key)
    if canonical_pos != None:
      name = row[canonical_pos]
      if name != MISSING: set_canonical(usage, name)
    if scientific_pos != None:
      sci = row[scientific_pos]
      if sci != MISSING: set_scientific(usage, sci)
    if rank_pos != None:
      rank = row[rank_pos]
      if rank != MISSING: set_rank(usage, rank)
    key_to_usage[key] = usage

    accepted_key = row[accepted_pos] if accepted_pos else MISSING
    if accepted_key == key: accepted_key = MISSING
    parent_key = row[parent_pos]
    if accepted_key != MISSING: parent_key = MISSING
    fixup.append((usage, parent_key, accepted_key))

  for (usage, parent_key, accepted_key) in fixup:
    if accepted_key != MISSING:
      probe2 = key_to_usage.get(accepted_key)
      if probe2: set_accepted(usage, probe2)
    elif parent_key != MISSING:
      probe1 = key_to_usage.get(parent_key)
      if probe1: set_parent(usage, probe1)

  # Collect children so we can say children[x]
  roots = collect_children(key_to_usage.values())

  # Prepare for doing within-tree MRCA operations
  cache_levels(roots)

  return (key_to_usage, roots)

# -----------------------------------------------------------------------------
# Get parent/child and accepted/synonym relationships

children_prop = prop.Property("children")
synonyms_prop = prop.Property("synonyms")
get_children = prop.getter(children_prop)
set_children = prop.setter(children_prop)
get_synonyms = prop.getter(synonyms_prop)
set_synonyms = prop.setter(synonyms_prop)

def get_inferiors(p):
  return (get_children(p, []) +
          get_synonyms(p, []))

def get_superior(x):
  sup = get_parent(x, None) or get_accepted(x, None)
  if sup:
    if get_level(sup, 0) != get_level(x, 1) - 1:
      print("# %s %s" % (get_level(sup), get_level(x)),
            file=sys.stderr)
      assert False
  return sup

def collect_children(items):
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

# -----------------------------------------------------------------------------
# MRCA

# Cache every node's level (distance to root)
#   simple recursive descent from roots

level_prop = prop.Property("level")
get_level = prop.getter(level_prop)
set_level = prop.setter(level_prop)

def cache_levels(roots):
  def cache(x, n):
    set_level(x, n)
    for c in get_inferiors(x):
      cache(c, n+1)
  for root in roots:
    cache(root, 1)

def find_peers(x, y):
  while get_level(x) < get_level(y):
    y = get_superior(y)
  while get_level(x) > get_level(y):
    x = get_superior(x)
  return (x, y)

def mrca(x, y):
  (x, y) = find_peers(x, y)
  while x and not (x is y):
    x = get_superior(x)
    y = get_superior(y)
  return x

# -----------------------------------------------------------------------------
# Sum file ingest

# Record-match file ingest

out_a_prop = prop.Property("out_a")
out_a = prop.getter(out_a_prop)
out_b_prop = prop.Property("out_b")
out_b = prop.getter(out_b_prop)
remark_prop = prop.Property("remark")
get_remark = prop.getter(remark_prop)
set_remark = prop.setter(remark_prop)

make_union = prop.constructor(key_prop, out_a_prop, out_b_prop,
                              remark_prop)


# Load record mathes from a file or whatever

def load_sum(iterator, a_usage_dict, b_usage_dict):
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

# Sadly, rows are repeated sometimes 

def note_match(x, y, remark, sum):
  global union_count
  assert isinstance(remark, str)
  (key_to_union, in_a, in_b) = sum
  if x:
    j = mep_get(in_a, x, None)
  elif y:
    j = mep_get(in_b, y, None)
  else:
    j = None
  if j:
    assert out_a(j) == x
    if out_b(j) != y:
      print("!! %s -> %s (%s)\n!! but := %s (%s)" %
            (get_key(x), get_key(out_b(j)), get_remark(j),
             get_key(y), remark),
            file=sys.stderr)
    assert out_b(j) == y
    assert mep_get(in_a, x, None) == j
    assert mep_get(in_b, y, None) == j
    return j
  # Invent a key?? and store it in the sum ???
  # To avoid an A/B conflict we'd need the B key->usage map?
  if y:
    key = get_key(y)
    if key in key_to_union: key = "B.%s" % key
  else:
    key = get_key(x)
    if key in key_to_union: key = "A.%s" % key
  j = make_union(key, x, y, remark)
  if x: mep_set(in_a, x, j)
  if y: mep_set(in_b, y, j)
  key_to_union[key] = j
  union_count += 1
  return j

union_count = 0

def combine_remarks(*remarks):
  return ";".join([r for r in remarks if r != MISSING])

# -----------------------------------------------------------------------------
# Analyze tipward record matches.

# A set of record matches is a 'tipward record match set' if each
# match (x, y) in the set has the property that no ancestor of either
# x or y is also matched in the set.

# The tipward record matches are stored in the given sum object.

def find_tipward_record_matches(a_roots, b_roots, rm_sum):

  def half_find_tipward(roots, inject, outject):
    tipwards = {}               # mep
    def traverse(x_usage):
      saw = False
      for child in get_inferiors(x_usage):
        saw = traverse(child) or saw
      if not saw:
        u = mep_get(inject, x_usage)
        if outject(u):
          if get_key(u).startswith(troublemaker):
            print("!! tipward: %s" % (get_key(u),), file=sys.stderr)
          mep_set(tipwards, u, u)
          return True
        else:
          return False
      return saw
    for root in roots:
      traverse(root)
    return tipwards

  (_, rm_in_a, rm_in_b) = rm_sum
  a_tipwards = half_find_tipward(a_roots, rm_in_a, out_b)
  b_tipwards = half_find_tipward(b_roots, rm_in_b, out_a)
  count = 0
  for u in a_tipwards.values():
    if mep_get(b_tipwards, u, None):
      if get_key(out_a(u)).startswith(troublemaker):
        print("!! %s = %s" % (get_key(out_a(u)), get_key(out_b(u))),
              file=sys.stderr)
      set_xmrca(out_a(u), out_b(u))
      set_xmrca(out_b(u), out_a(u))
      count += 1
  if len(a_tipwards) != count:
    print("-- align: %s tipward record matches in A" % len(a_tipwards),
          file=sys.stderr)
  if len(b_tipwards) != count:
    print("-- align: %s tipward record matches in B" % len(b_tipwards),
          file=sys.stderr)
  print("-- align: %s net tipward record matches" % count, file=sys.stderr)

# -----------------------------------------------------------------------------
# 7. how-related WITHIN hierarchy

EQ = 1
LT = 2
GT = 3
DISJOINT = 4
CONFLICT = 5
UNCLEAR = 6     # co-synonyms of the same accepted
rcc5_symbols = ['---', '=', '<', '>', '!', '><', '?']

def how_related(x, y):
  (x1, y1) = find_peers(x, y)
  if x1 == y1:
    if x == y:
      return EQ
    elif x1 == x:
      return GT
    else:
      return LT
  elif (get_parent(x1, 123) == get_parent(y1, 456) and
        get_accepted(x1, None) or get_accepted(y1, None)):
    return UNCLEAR
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
    assert get_level(m) > 0
    return m                    # Tipward record match
  for c in get_inferiors(x):
    n = half_xmrcas(c)
    if n:
      assert get_level(n) > 0
      if m:
        m = mrca(m, n)
      else:
        m = n
  if m:
    if get_key(x).startswith(troublemaker):
        print("!! xmrca(%s) = %s" % (get_key(x), get_key(m)),
              file=sys.stderr)
    set_xmrca(x, m)
  return m

def cache_xmrcas(a_roots, b_roots):
  for root in a_roots:
    half_xmrcas(root)
  for root in b_roots:
    half_xmrcas(root)
  if True:
    def check(x):
      for c in get_inferiors(x): check(c)
      y = get_xmrca(x, None)
      if y and not get_xmrca(y, None):
        print("!! x %s %s y %s %s" %
              (get_key(x), get_canonical(x), get_key(y), get_canonical(x),),
              file=sys.stderr)
        assert False
    for root in a_roots: check(root)
    for root in b_roots: check(root)

# -----------------------------------------------------------------------------

# Compare an A node with a B node, obtaining the RCC5 relationship
# that holds between them.
# Return value is (rcc5, e, f, g) where:
#   e is <= both p and q
#   f is <= p but ! q
#   g is <= q but ! p

def related_how(p, q):
  pa = get_accepted(p, p)
  qa = get_accepted(q, q)
  result = related_how_accepted(pa, qa)
  (rcc5, _, _, _) = result
  if (rcc5 == EQ and
      (pa != p or qa != q) and
      not (get_xmrca(p, None) == q and
           get_xmrca(q, None) == p)):
    return (UNCLEAR, _, _, _)
  else:
    return result

def related_how_accepted(p, q):
  c = get_xmrca(p, None)
  d = get_xmrca(q, None)
  if c == None or d == None:
    return (DISJOINT, None, p, q)
  elif get_xmrca(c, None) == d:      # c <= d?
    # They belong to a 'monotype' cluster.
    # Not consistent with what 'weave' does.
    # I hope the result doesn't matter a whole lot.
    cmp = ((get_level(d) - get_level(p)) -
           (get_level(c) - get_level(q)))
    if cmp < 0:
      return (LT, p, None, True)
    elif cmp == 0:
      return (EQ, p, None, None)
    else:
      return (GT, p, True, None)
  else:
    #  e   f   g
    # yes yes yes  CONFLICT
    # yes yes no   GT
    # yes no  no   EQ
    # yes no  yes  LT
    # no  yes yes  DISJOINT
    #  the other 3 cases can't occur
    (e, f) = seek_conflict(p, q)       # Both in A
    if not e:
      assert f
      return (DISJOINT, None, f, q)    # g would be nice
    elif f:
      (e2, g) = seek_conflict(q, p)   # Both in B
      assert e2
      if g:
        return (CONFLICT, e, f, g)    # or (e, f, g)?
      else:
        return (GT, e, f, None)
    else:
      # EQ case is covered above
      return (LT, e, None, True)    # g would be nice

# p < xmrca(q); but is p < q?

def seek_conflict(p, q):
  o = get_xmrca(p, None)
  if o == None: return (None, None)
  rel = how_related(o, q)
  if rel == GT:
    p_and_q = p_not_q = None                # in A
    for x in get_inferiors(p):
      (x_and_q, x_not_q) = seek_conflict(x, q)
      if x_and_q and x_not_q:
        return (x_and_q, x_not_q)
      if x_and_q: 
        p_and_q = x_and_q
      if x_not_q:
        p_not_q = x_not_q
      if p_and_q and p_not_q:             # hack: cut it short
        return (p_and_q, p_not_q)
    return (p_and_q, p_not_q)
  elif rel == DISJOINT:
    return (None, q)
  else:  # LT, EQ, UNCLEAR
    return (p, None)

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

def set_congruences(a_roots, b_roots, sum, rm_sum):
  (_, in_a, in_b) = sum
  (_, rm_in_a, rm_in_b) = rm_sum

  connect = connector(rm_sum, sum)

  def weave(u, v):
    # v and u need to be linked.
    if mep_get(in_b, v, None) != None:
      return

    # Set the comment
    
    j = mep_get(rm_in_a, u)
    if j == mep_get(rm_in_b, v):
      g = connect(u, v, get_remark(j))
    else:
      g = connect(u, v, "join point")

    s = get_superior(u)
    t = get_superior(v)
    k = g
    while True:
      s_done = (not s or get_xmrca(s, None) != v)
      t_done = (not t or get_xmrca(t, None) != u)
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
          elif m and get_xmrca(m, None) == u and get_level(m) < get_level(t):
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
    #return (k, s, t)

  def traverse(x):
    # Deal with descendants of x, and inject x, bottom up
    for c in get_inferiors(x): traverse(c)
    # Skip if already processed
    if mep_get(in_a, x, None) == None:
      v = get_xmrca(x, None)      # in B
      if v != None:
        u = get_xmrca(v, None)          # in A
        if get_key(x).startswith(troublemaker):
          print("# x %s v %s u %s" % (get_key(x), get_key(v), get_key(u)),
                file=sys.stderr)
        if u == x:  # and how_related(x, u) != LT:
          # u and v will be connected
          weave(u, v)
          # k is the rootward cluster node in sum

  for x in a_roots: traverse(x)

# Get remarks from rm_sum ... ?

def connector(rm_sum, sum):
  (_, rm_in_a, rm_in_b) = rm_sum
  def connect(x, y, remark):
    # Get remark from prior record match
    if x:
      j = mep_get(rm_in_a, x, None)
    elif y:
      j = mep_get(rm_in_b, y, None)
    assert j
    assert remark
    remarks = [remark]
    if out_a(j) == x and out_b(j) == y:
      rem = get_remark(j)
      if not rem in remarks:
        remarks.append(rem)
    return note_match(x, y, combine_remarks(*remarks), sum)
  return connect

# -----------------------------------------------------------------------------
# Now set the parent pointers for the merged tree

def build_tree(a_roots, b_roots, sum, rm_sum):
  (key_to_union, in_a, in_b) = sum
  (_, rm_in_a, rm_in_b) = rm_sum
  connect = connector(rm_sum, sum)
  def fasten_a(x):
    return connect(x, None, "peripheral in A")
  def fasten_b(y):
    return connect(None, y, "peripheral in B")
  finish_sum(b_roots,
             in_b, in_a, fasten_b, fasten_a, True, rm_in_b)
  finish_sum(a_roots,
             in_a, in_b, fasten_a, fasten_b, False, rm_in_a)
  print("# %s in a, %s in b, %s union keys" %
        (len(in_a), len(in_b), len(key_to_union)),
        file=sys.stderr)

  unions = key_to_union.values()
  roots = collect_children(unions)
  print("-- align: %s nodes, %s roots in sum" % (len(unions), len(roots)),
        file=sys.stderr)
  return roots

# Set parents of nodes whose parent wasn't set by 'weave'

def finish_sum(a_roots, in_a, in_b, fasten_a, fasten_b, priority, rm_in_a):
  stats = [0, 0, 0, 0, 0]
  def cap(x, in_a, fasten_a):
    j = mep_get(in_a, x, None)
    if not j:
      j = fasten_a(x)
      assert mep_get(in_a, x) == j
      stats[2] += 1
    return j

  def traverse(x):
    stats[0] += 1
    j = cap(x, in_a, fasten_a)
    for c in get_inferiors(x):
      if get_xmrca(c, None): traverse(c)
    for c in get_inferiors(x):
      if not get_xmrca(c, None): traverse(c)
    # All joint nodes derived from nodes descended from x
    # have their superior node set at this point.
    if not get_superior(j):
      (pq, ab) = determine_superior_in_sum(x, priority)
      if pq:
        if ab:
          h = cap(pq, in_a, fasten_a)
        else:
          h = cap(pq, in_b, fasten_b)
        set_superior(j, h)
        stats[1] += 1
      else:
        stats[3] += 1
    else:
      stats[4] += 1

  def determine_superior_in_sum(x, priority):
    # Two candidates for superior: superior in a, and xmrca in b
    p = get_superior(x)         # p is in a
    q = get_xmrca(x, None)      # q is in b
    if q == None:
      # No tipward matches in sight under x.
      # If there is a record match to x, assume x is splitting it, and
      # put x under it.
      m = out_b(mep_get(rm_in_a, x), None)
      if (m and p and get_xmrca(m, None) and
          # Make sure p [parent of x in a] ≈ parent of m in b.
          mep_get(in_a, p) == mep_get(in_b, get_superior(m))):
        print("# Putting %s under %s" % (get_blurb(x), get_blurb(m)),
              file=sys.stderr)
        return (m, False)
      else:
        return (p, True)
    u = get_xmrca(q, None)
    assert u
    # Ascend until we find a B-side ancestor that's definitely bigger than x
    # This isn't right for monotype chains but those are handled separately
    while q and get_xmrca(q, None) == u:
      q = get_superior(q)
    # Ascend until we find a priority-side ancestor that doesn't conflict with p
    if p == None:
      return (q, False)
    while q:
      # No monotype chains.  Ordering is by tipward-match sets.
      if q == None:
        return (p, True)
      else:
        (rel, _, _, _) = related_how(p, q)
        if rel == GT:
          return (q, False)
        elif rel == LT or rel == EQ:
          return (p, True)
        elif rel == CONFLICT:
          if False:
            print("!! %s ? %s" % (get_key(p), get_key(q)),
                file=sys.stderr)
          if priority:
            # Ignore conflicting low-priority candidate q
            return (p, True)
          else:
            # Seek a nonconflicting high-priority candidate
            pass
        else:                   # DISJOINT or UNCLEAR
          assert False
      q = get_superior(q)
    return (None, True)

  for root in a_roots: traverse(root)
  print("# Finish: touched %s, set sup %s, capped %s, orphans %s, pass %s" % tuple(stats),
        file=sys.stderr)
  # end finish_sum


# j is to be either a child or synonym of k.  Figure out which.
# Invariant: A synonym must have neither children nor synonyms.

def set_superior(j, k):
  assert k
  assert not get_superior(j)
  assert not get_accepted(k, None)
  x = out_a(j)
  y = out_b(j)
  if True:
    if y:
      if get_accepted(y, None):
        assert not get_children(j, None)
        assert not get_synonyms(j, None)
        set_accepted(j, k)
      else:
        set_parent(j, k)
    else:  # x
      if get_accepted(x, None):
        assert not get_children(j, None)
        assert not get_synonyms(j, None)
        set_accepted(j, k)
      else:
        set_parent(j, k)
  else:
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
# Report on differences between record matches and hierarchy matches

def report(rm_sum, sum):
  (rm_map, rm_in_a, rm_in_b) = rm_sum
  (a_map, in_a, in_b) = sum
  drop_a = 0
  for r in rm_map.values():     # Record match
    x = out_a(r)
    if x:
      zu = mep_get(in_a, x)     # Bridge
      y = out_b(r)              # x record matches y
      z = out_b(zu)             # x aligns to z
      if y != z:
        # Record match != aligned
        if y and z:
          # related_how(x in A, y in B)
          (rcc5, _, _, _) = related_how(x, y)
          print("  %s %s %s" %
                (get_blurb(mep_get(in_b, y)),
                 rcc5_symbols[rcc5],
                 get_blurb(zu)),
                file=sys.stderr)
        else:
          drop_a += 1
  if drop_a > 0:
    print("-- align: %s broken in A" % (drop_a,), file=sys.stderr)
    
def get_blurb(r):
  if r:
    name = get_canonical(r, None) or get_scientific(r, None)
    if name != None:
      if get_accepted(r, None):
        return name + " (synonym)"
      else:
        return name
    return "[no canonical]"
  else:
    return "[no match]"

# -----------------------------------------------------------------------------
# 14. emit the new sum with additional parent column

# Returns a row generator

def generate_sum(sum, roots, rm_sum):
  yield ["taxonID", "taxonID_A", "taxonID_B",
         "parentNameUsageID", "acceptedNameUsageID",
         "canonicalName", "previousID", "change", "remark"]
  (_, in_a, in_b) = sum
  (_, rm_in_a, rm_in_b) = rm_sum

  # was: for union in all_unions(sum): ...

  def compose_row(union):
    a_usage = out_a(union, None)
    b_usage = out_b(union, None)
    p = get_parent(union, None)
    a = get_accepted(union, None)
    z = m = None
    sym = MISSING
    if b_usage:
      z = out_a(mep_get(rm_in_b, b_usage), None)
      if z:
        (rcc5, _, _, _) = related_how(z, b_usage)
        sym = rcc5_symbols[rcc5]
        m = mep_get(in_a, z, None)
      # if a_usage and b_usage then "renamed"
      # if b_usage then if get_xmrca(b_usage) then "collected"
    else:
      z = out_b(mep_get(rm_in_a, a_usage), None)
      if z:
        (rcc5, _, _, _) = related_how(a_usage, z)
        sym = rcc5_symbols[rcc5]
        m = mep_get(in_b, z, None)
      # if a_usage then if get_xmrca(a_usage) then "dissolved"
    return [get_key(union),
            get_key(a_usage) if a_usage else MISSING,
            get_key(b_usage) if b_usage else MISSING,
            get_key(p) if p else MISSING,
            get_key(a) if a else MISSING,
            get_canonical(union, MISSING),
            get_key(m) if m else MISSING,
            sym,
            get_remark(union, MISSING)]
  def traverse(union):
    yield compose_row(union)
    for inf in get_inferiors(union):
      for row in traverse(inf): yield row
  for root in roots:
    for row in traverse(root): yield row

def all_unions(sum):
  (key_to_union, _, _) = sum
  return key_to_union.values()

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    TBD.  stdin = A hierarchy
    """)
  parser.add_argument('--target', help="B hierarchy")
  parser.add_argument('--matches', help="record matches")
  args=parser.parse_args()

  a_file = sys.stdin
  b_path = args.target
  rm_sum_path = args.matches

  with open(b_path) as b_file:
    # TBD: Compute sum if not provided ?
    with open(rm_sum_path) as rm_sum_file:
      sum = align(csv.reader(a_file),
                  csv.reader(b_file),
                  csv.reader(rm_sum_file))
      util.write_generated(sum, sys.stdout)
