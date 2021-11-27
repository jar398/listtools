#!/usr/bin/env python3

# Align two trees (given as DwC files i.e. usage lists), with
# sensitivity to parent/child relations.

import sys, csv, argparse
import util
from util import windex, MISSING
import property as prop
from property import mep_get, mep_set

TOP = "[top]"
debug = False

def is_top(x): return get_key(x) == TOP

def monitor(x):
  # (x and (get_canonical(x, None) == "Myomyscus angolensis"))
  return False

# -----------------------------------------------------------------------------
# Supervise the overall process

def align(a_iterator, b_iterator, rm_sum_iterator):
  (a_usage_dict, a_roots) = load_usages(a_iterator)
  (b_usage_dict, b_roots) = load_usages(b_iterator)

  # Read record matches
  rm_sum = load_sum(rm_sum_iterator, a_usage_dict, b_usage_dict)

  calculate_xmrcas(a_roots, b_roots, rm_sum)

  # Merge the two trees
  sum = ({}, {}, {})
  set_congruences(a_roots, b_roots, sum, rm_sum)
  roots = build_tree(a_roots, b_roots, sum, rm_sum)

  # Prepare to assign names
  assign_canonicals(sum)

  # Emit tabular version of merged tree
  return generate_alignment(sum, roots, rm_sum)

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

year_prop = prop.Property("year")
get_year = prop.getter(year_prop)
set_year = prop.setter(year_prop)

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
  year_pos = windex(header, "year")
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
    if year_pos != None:
      year = row[year_pos]
      if year != MISSING: set_year(usage, year)
    key_to_usage[key] = usage

    accepted_key = row[accepted_pos] if accepted_pos else MISSING
    if accepted_key == key: accepted_key = MISSING
    parent_key = row[parent_pos]
    if accepted_key != MISSING: parent_key = MISSING
    assert parent_key != key
    fixup.append((usage, parent_key, accepted_key))

  for (usage, parent_key, accepted_key) in fixup:
    if accepted_key != MISSING:
      probe2 = key_to_usage.get(accepted_key)
      if probe2:
        set_accepted(usage, probe2)
        if monitor(usage) or monitor(probe1):
          print("> accepted %s := %s" % (get_blurb(usage), get_blurb(probe1)),
                file=sys.stderr)
    elif parent_key != MISSING:
      probe1 = key_to_usage.get(parent_key)
      if probe1:
        set_parent(usage, probe1)
        if monitor(usage) or monitor(probe1):
          print("> parent %s := %s" % (get_blurb(usage), get_blurb(probe1)),
                file=sys.stderr)

  # Collect children so we can say children[x]
  roots = collect_inferiors(key_to_usage.values())

  # We really only want one root (this is so that mrca can work)
  if True or len(roots) > 1:
    top = make_usage(TOP)
    set_canonical(top, TOP)
    key_to_usage[TOP] = top
    for root in roots: set_parent(root, top)
    set_children(top, roots)
    roots = [top]

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

# E.g. Glossophaga bakeri (in 1.7) split from Glossophaga commissarisi (in 1.6)

def collect_inferiors(items):
  roots = []
  seniors = 0
  for item in items:

    # Add item to list of accepted's synonyms
    accepted_item = get_accepted(item, None)
    if accepted_item:
      # Filter out senior synonyms here
      if seniority(item, accepted_item) == "senior synonym":
        seniors += 1
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

  if False and seniors > 0:     # Maybe interesting
    print("-- Saw %s senior synonyms" % seniors,
          file=sys.stderr)
  return roots

# Is given synonym usage a senior synonym of its accepted usage?
# In the case of splitting, we expect the synonym to be a senior
# synonym of the item.

# We could also look to see if taxonomicStatus is 'senior synonym'.

# Hmm, if there's exactly one senior synonym, we should remember it as
# being the 'split from' taxon.

def seniority(item, accepted_item):
  year1 = get_year(item, None)  # the synonym
  year2 = get_year(accepted_item, None)
  if year1 and year2:
    if year1 <= year2:
      return "senior synonym"
    else:
      return "junior synonym"
  else:
    return "synonym"

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
  assert x and y
  while get_level(x) < get_level(y):
    y = get_superior(y)
  while get_level(x) > get_level(y):
    x = get_superior(x)
  return (x, y)

def mrca(x0, y0):
  (x, y) = find_peers(x0, y0)
  (x0, y0) = (x, y)    # save for error messages
  while not (x is y):
    x = get_superior(x)
    if not x:
      print("!! mrca off top: x0 %s y0 %s" % (get_key(x0), get_key(y0)),
            file=sys.stderr)
      break
    y = get_superior(y)
    if not y:
      print("!! mrca top off: x0 %s y0 %s" % (get_key(x0), get_key(y0)),
            file=sys.stderr)
      break
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

  # These pass somehow
  a_top = a_usage_dict.get(TOP)
  b_top = b_usage_dict.get(TOP)
  assert is_top(a_top)
  assert is_top(b_top)

  note_match(a_usage_dict.get(TOP),
             b_usage_dict.get(TOP),
             'top', sum)

  # More sanity checks
  (_, rm_in_a, rm_in_b) = sum
  top = mep_get(rm_in_a, a_top)
  assert mep_get(rm_in_b, b_top) == top
  assert out_a(top) == a_top
  assert out_b(top) == b_top

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

  if monitor(x) or monitor(y):
    print("> note match %s to %s; %s" %
          (get_blurb(x), get_blurb(y), remark),
          file=sys.stderr)

  return j

union_count = 0

def combine_remarks(*remarks):
  return ";".join([r for r in remarks if r != MISSING])

# -----------------------------------------------------------------------------
# 8. Cross-mrcas

xmrca_prop = prop.Property("xmrca")
get_xmrca = prop.getter(xmrca_prop)
set_xmrca = prop.setter(xmrca_prop)

def calculate_xmrcas(a_roots, b_roots, rm_sum):
  (_, rm_in_a, rm_in_b) = rm_sum
  def half_xmrcas(roots, record_match):
    def traverse(x):
      count = 0
      y = None
      for c in get_inferiors(x):
        n = traverse(c)
        if n:
          if y:
            y = mrca(y, n)
          else:
            y = n
          count += 1
      if y == None:
        y = record_match(x)     # See if this is a tipward match
      if y:
        if monitor(x):
            print("> xmrca(%s) := %s" % (get_blurb(x), get_blurb(y)),
                  file=sys.stderr)
        set_xmrca(x, y)
      return y
    for root in roots: traverse(root)
  half_xmrcas(a_roots, lambda x:out_b(mep_get(rm_in_a, x)))
  half_xmrcas(b_roots, lambda y:out_a(mep_get(rm_in_b, y)))

  fix_xmrcas(a_roots, b_roots, rm_sum)
  check_xmrcas(a_roots, b_roots)

# Use record matches to resolve correspondence within clusters

def fix_xmrcas(a_roots, b_roots, rm_sum):
  (_, rm_in_a, rm_in_b) = rm_sum
  def half_fix(roots, record_match):
    def traverse(x):
      y = record_match(x)
      if y:
        y0 = get_xmrca(x)
        x0 = get_xmrca(y0)
        if le(x0, x):
          # assert get_xmrca(x0) == y0 !! why not?
          if get_xmrca(y) == x0:
            # y is in the x0/y0 cluster
            if y != y0:
              z = x
              while z and get_xmrca(z) == y0:
                if monitor(x):
                 print("# Fixing A %s: %s -> %s" %
                      (get_blurb(z), get_blurb(get_xmrca(z)), get_blurb(y)),
                      file=sys.stderr)
                set_xmrca(z, y)
                z = get_superior(z)
            if x != x0:
              z = y
              while z and get_xmrca(z) == x0:
                if monitor(x):
                 print("# Fixing B %s: %s -> %s" %
                      (get_blurb(z), get_blurb(get_xmrca(z)), get_blurb(x)),
                      file=sys.stderr)
                set_xmrca(z, x)
                z = get_superior(z)
      for c in get_inferiors(x):
        traverse(c)
    for root in roots: traverse(root)
  half_fix(a_roots, lambda x:out_b(mep_get(rm_in_a, x)))
  half_fix(b_roots, lambda y:out_a(mep_get(rm_in_b, y)))
  a_top = a_roots[0]
  b_top = b_roots[0]
  assert get_xmrca(a_top) == b_top
  assert get_xmrca(b_top) == a_top # FAILS

def check_xmrcas(a_roots, b_roots):
  if True:
    def check(x):
      for c in get_inferiors(x): check(c)
      y = get_xmrca(x, None)
      if y and not get_xmrca(y, None):
        print("!! non-mutual xmrca: x %s %s y %s %s" %
              (get_key(x), get_canonical(x), get_key(y), get_canonical(x),),
              file=sys.stderr)
        assert False
    for root in a_roots: check(root)
    for root in b_roots: check(root)

# -----------------------------------------------------------------------------
# 7. how-related WITHIN hierarchy

EQ = 1
LT = 2
GT = 3
DISJOINT = 4
CONFLICT = 5
UNCLEAR = 6     # co-synonyms of the same accepted
LT_OR_CONFLICT = 7
GT_OR_CONFLICT = 8
rcc5_symbols = ['---', '=', '<', '>', '!', '><', '?', '<?', '>?']

def how_related(x, y):
  (x1, y1) = find_peers(x, y)
  if x1 == y1:
    if x == y:
      return EQ
    elif x1 == x:
      return GT
    else:
      return LT
  elif ((get_accepted(x1, None) or get_accepted(y1, None)) and
        get_superior(x1) == get_superior(y1)):
    return UNCLEAR
  else:
    return DISJOINT

# -----------------------------------------------------------------------------

# Compare an A node with a B node, obtaining the RCC5 relationship
# that holds between them.
# Return value is (rcc5, e, f, g) where:
#   e is <= both p and q
#   f is <= p but ! q
#   g is <= q but ! p

def related_how(p, q):
  result = really_related_how(p, q)
  if False:
    result2 = really_related_how(q, p)
    want = result2[0]
    if want == LT: want = GT
    elif want == GT: want = LT
    if want != result[0]:
      print("!! Comparison is asymmetric: %s %s %s because %s\n   vs. %s because %s" %
            (get_blurb(p), rcc5_symbols[result[0]], get_blurb(q),
             result[1], rcc5_symbols[result2[0]], result2[1]),
            file=sys.stderr)
  if (result[0] == DISJOINT and
      (get_accepted(p, None) or get_accepted(q, None)) and
      get_superior(p) == get_superior(q)):
    result = (UNCLEAR, "synonym", None, None)
  if monitor(p) or monitor(q):
    tag = result[1] if isinstance(result[1], str) else "..."
    print("> %s %s %s because %s" %
          (get_blurb(p), rcc5_symbols[result[0]], get_blurb(q), tag),
          file=sys.stderr)
  return result

def really_related_how(p, q):

  d = get_xmrca(p, None)        # p ? d
  c = get_xmrca(q, None)        # q ? c

  if c == None: return (DISJOINT, "q is peripheral", p, q)
  if d == None: return (DISJOINT, "p is peripheral", p, q)

  # p0 is the xmrca round-trip via d
  # q0 is the xmrca round-trip via c
  p0 = get_xmrca(d, None)             # d <= p0
  q0 = get_xmrca(c, None)             # c <= q0

  if p0 == None: return (DISJOINT, "q is peripheral?", p, q)
  if q0 == None: return (DISJOINT, "p is peripheral?", p, q)

  # We can now test relative order of p, c, p0 and q, d, q0
  # subject to d <= p0 and c <= q0

  # p ? c <= q0 ? q      and vice versa, p ? p0 >= d ? q

  p0_le_p = le(p0, p)
  q0_le_q = le(q0, q)
  p_lesseq_q = (le(p, c) and q0_le_q)
  q_lesseq_p = (le(q, d) and p0_le_p)

  if p_lesseq_q:
    if p != c:
      return (LT, "p < c <= q0 <= q", None, None)
    elif q0 != q:
      return (LT, "p <= c <= q0 < q", None, None)
    elif q_lesseq_p:
      return (EQ, "p <= c <= q0 <= q <= d <= p0 <= p", None, None)
  if q_lesseq_p:
    if q != d:
      return (GT, "q < d <= p0 <= p", None, None)
    elif p0 != p:
      return (GT, "q <= d <= p0 < p", None, None)

  if p0_le_p and p0 == c and q0_le_q and q0 == d:
    return compare_in_cluster(p, p0, q, q0)

  return compare_carefully(p, q)

def compare_carefully(p, q):

  # Brute force
  # Check if get_xmrca(p) > q ... ?
  (f1, e) = seek_conflict(p, q)
  (f2, g) = seek_conflict(q, p)

  if f1 or f2:
    if e:
      if g:
        return (CONFLICT, "on careful examination: %s %s %s %s" %
                (get_blurb(e), get_blurb(f1), get_blurb(f2),
                 get_blurb(g)), e, g)
      else:
        return (GT, "checked", None, True)
    elif g:
      return (LT, "checked", True, None)
    else:
      if False and not (f1 and f2):
        print("!! Missing f1/f2. %s ? %s" % (get_blurb(p), get_blurb(q)),
              file=sys.stderr)
      return (EQ, "checked ??", True, None) # shouldn't happen?
  elif e or g:
    return (DISJOINT, "checked", p, q)
  else:
    if monitor(p) or monitor(q):
      print("!! Unclear. %s ? %s" % (get_blurb(p), get_blurb(q)),
            file=sys.stderr)
    return (UNCLEAR, "checked", None, None)

# Returns two sibling(?) descendants of p
# Assumes initially that either p > q or p >< q

def seek_conflict(p, q):
  o = get_xmrca(p, None)
  if o == None: return (None, None)

  rel = how_related(o, q)
  if rel == LT:
    return (p, None)
  elif rel == DISJOINT:
    return (None, p)
  elif rel == UNCLEAR:
    if monitor(p) or monitor(q):
      print("!! Equivocal: %s ? %s" %
            (get_blurb(p), get_blurb(q)),
            file=sys.stderr)
    return (None, None)
  else:    # GT EQ
    p_and_q = p_not_q = None                # in A
    for x in get_children(p, []):           # Disjoint
      (x_and_q, x_not_q) = seek_conflict(x, q)
      if x_and_q: p_and_q = x_and_q
      if x_not_q: p_not_q = x_not_q
      if p_and_q and p_not_q:             # hack: cut it short
        return (p_and_q, p_not_q)
    if p_and_q or p_not_q:
      return (p_and_q, p_not_q)
    # No information from children.  Look at the synonyms.
    for x in get_synonyms(p, []):
      (x_and_q, x_not_q) = seek_conflict(x, q)
      if x_and_q: p_and_q = x_and_q
      if x_not_q: p_not_q = x_not_q
      if p_and_q and p_not_q:             # hack: cut it short
        if monitor(p) or monitor(q):
          print("!! Inconclusive evidence for %s vs. %s\n   %s ? %s" %
                (get_blurb(p), get_blurb(q),
                 get_blurb(p_and_q), get_blurb(p_not_q)),
                file=sys.stderr)
        return (None, None)  # was: return (None, p_not_q)
    return (p_and_q, p_not_q)


def compare_in_cluster(p, p0, q, q0):
  # if x > y then level x < level y
  cmp = ((get_level(p0) - get_level(p)) -
         (get_level(q0) - get_level(q)))
  print("# Kludge: %s - %s = %s" % (get_blurb(p), get_blurb(q), cmp),
        file=sys.stderr)
  if cmp == 0:
    return (EQ, "Kludge: = in cluster", None, None)
  elif cmp > 0:
    return (GT, "Kludge: > in cluster", True, None)
  else:
    return (LT, "Kludge: < in cluster", None, True)

def compare_in_cluster_2(p, p0, q, q0):
  if p == p0 and q == q0:
    return (EQ, "p <= p0 <= d = q <= q0 <= c = p", True, None)
  if p == p0:
    return (LT, "< in cluster", None, True)
  elif q == q0:                 # q > q0
    return (GT, "> in cluster", True, None)
  elif is_top(p):
    if is_top(q):
      return (EQ, "top = top", True, None)
    else:
      return (GT, "top > q", True, None)
  elif is_top(q): return (LT, "p < top", None, True)
  elif get_canonical(p) == get_canonical(q):
    # TEMP KLUDGE
    return (EQ, "name= kludge in cluster", None, None)
  else:
    print("!! Don't know how to compare %s to %s\n   p0=%s q0=%s" %
          (get_blurb(p), get_blurb(q), get_blurb(p0), get_blurb(q0)),
          file=sys.stderr)
    return ((LT, "name< kludge in cluster", None, True)
            if get_canonical(p) < get_canonical(q)
            else (GT, "name> kludge in cluster", True, None))

def le(x, y): return mrca(x, y) == y
def lt(x, y): return le(x, y) and x != y

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

  u  =  g  =  v   tipward nodes in cluster.  u <= x
  |  \     /  |
  |           |
  |     |     |
  |     k     |   s, t, and k are the iteration variables
  |  /     \  |    they progress downwards
  |           t
  s  = (r) =  |
  |           m   record match to s
  |           |
  x  ?  k  ?  y   rootward nodes in cluster
  ^     |     ^
  p  <= h =>  q   parents outside the cluster
"""

def set_congruences(a_roots, b_roots, sum, rm_sum):
  (_, in_a, in_b) = sum
  (_, rm_in_a, rm_in_b) = rm_sum

  connect = connector(rm_sum, sum)
  assert is_top(a_roots[0])
  assert is_top(b_roots[0])

  def weave(u, v):
    if debug:
      print("# weaving %s with %s" % (get_blurb(u), get_blurb(v)), file=sys.stderr)
    # v and u need to be linked.
    if mep_get(in_b, v, None):
      return

    # Set the comment based on the record match
    j = mep_get(rm_in_a, u)
    if j == mep_get(rm_in_b, v):
      g = connect(u, v, get_remark(j))
    else:
      g = connect(u, v, "join point")

    s = get_superior(u)
    t = get_superior(v)
    k = g
    while True:
      s_done = (not s or get_xmrca(s, None) != v) # x
      t_done = (not t or get_xmrca(t, None) != u) # y
      if s_done and t_done:
        break
      elif s_done:              # t not done, deal with t
        k = connect(None, t, "right chain")
        t = get_superior(t)
      elif t_done:              # s not done, deal with s
        k = connect(s, None, "left chain")
        s = get_superior(s)
      else:
        # Use record matching to connect cluster nodes when possible
        r = mep_get(rm_in_a, s)
        if r:
          m = out_b(r)
          if m == t:
            assert False        # Obsolete logic
            k = connect(s, m, "record match in cluster")
            s = get_superior(s)
            t = get_superior(t)
          elif m and get_xmrca(m, None) == u and get_level(m) < get_level(t):
            k = connect(None, t, "pursuing record match in cluster")
            t = get_superior(t)
          else:
            k = connect(s, None, "arbitrary superior choice")
            # Could choose t, I think
            s = get_superior(s)
        else:    # no record match
          k = connect(s, None, "no record matched in cluster")
          # Could choose t, I think
          s = get_superior(s)
      set_superior(g, k)
      if debug:
        print("# weave super %s = %s" % (get_blurb(j), get_blurb(k)), file=sys.stderr)
      g = k

  def traverse(x):
    # Skip if already processed
    if mep_get(in_a, x, None) == None:

      y = get_xmrca(x, None)      # in B
      if y != None and x == get_xmrca(y, None):
        # Mutual xmrca with xmrca... that means they're equivalent
        if monitor(x):
          print("> weave(%s, %s)" % (get_blurb(x), get_blurb(y)),
                file=sys.stderr)
        weave(x, y)
      else:
        if monitor(x):
          print("> Not weaving(%s, %s)" % (get_blurb(x), get_blurb(y)),
                file=sys.stderr)

    # Deal with descendants of x, and inject x, bottom up
    for c in get_inferiors(x): traverse(c)

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
  roots = collect_inferiors(unions)
  print("-- align: %s usages, %s roots in sum" % (len(unions), len(roots)),
        file=sys.stderr)
  return roots

# Set parents of nodes whose parent wasn't set by 'weave'

def finish_sum(a_roots, in_a, in_b, fasten_a, fasten_b, priority, rm_in_a):
  stats = [0, 0, 0, 0, 0]
  def cap(x, in_a, fasten_a):
    j = mep_get(in_a, x, None)
    if not j:
      j = fasten_a(x)
      stats[2] += 1
    assert mep_get(in_a, x) == j
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
      if monitor(x):
        print("> superior of %s will be %s in %s checklist" %
              (get_blurb(x), get_blurb(pq), "same" if ab else "other"),
              file=sys.stderr)

      if pq:                    # parent is either p or q
        if ab:                  # parent is in A, not in B
          h = cap(pq, in_a, fasten_a)
        else:
          h = cap(pq, in_b, fasten_b)
        set_superior(j, h)
        if debug:
          print("# finish super %s = %s" % (get_blurb(j), get_blurb(h)), file=sys.stderr)
        stats[1] += 1
      else:
        if debug:
          print("# hmm pq %s no superior why not?" % (get_blurb(j)), file=sys.stderr)
        stats[3] += 1
    else:
      stats[4] += 1

    if not (get_superior(j) or is_top(j)):
      print("!! Missing parent: %s = %s + %s" %
            (get_blurb(out_a(j)), get_blurb(out_b(j)), get_blurb(j)),
            file=sys.stderr)
      assert False

  def ad_hoc_split(x, p):
    if True: return False
    # If there is a record match to x, assume x is splitting it, and
    # put x under it.
    m = out_b(mep_get(rm_in_a, x), None)
    if (m and p and get_xmrca(m, None) and
        # Make sure p [parent of x in a] ≈ parent of m in b.
        mep_get(in_a, p, 123) == mep_get(in_b, get_superior(m), 456)):
      print("# Inferring that %s split off, sibling to %s" %
            (get_blurb(get_inferiors(m)[0]), get_blurb(m),),
            file=sys.stderr)
      if debug:
        print("# A-parent of %s is %s" % (get_blurb(x), get_blurb(m)), file=sys.stderr)
      return m
    else:
      return False

  def determine_superior_in_sum(x, priority):

    if is_top(x): return (None, True)

    # Find minimal p, q such that either x < p <= q or x < q <= p.
    # If not priority, also require q does not conflict with p.
    p = get_superior(x)         # p is in a (possibly top)
    assert lt(x, p)
    m = ad_hoc_split(x, p)
    if m: return (m, False)

    # This method can probably be made much more efficient, but for
    # now at least it seems to work

    q = get_xmrca(x, None)      # q is in b
    if q == None: return (p, True)    # peripheral
    # q = get_accepted(q, q)     # parent can't be a synonym (Rhinolophus hirsutus)

    # Take q up the b-lineage until x < q.
    # If not priority, further ascend b-lineage until not (p >< q).
    while related_how(x, q)[0] != LT:
      # wait what if x is a synonym?
      # if priority, then stop at conflict.
      if is_top(q):    # don't go up from top
        print("!! %s %s %s" % (get_blurb(x),
                               rcc5_symbols[related_how(x, q)[0]],
                               get_blurb(q)),
              file=sys.stderr)
        assert False
      q = get_superior(q)
    if not priority:
      # a tree (p, x) is low-priority.  Climb up to avoid conflict.
      while related_how(p, q)[0] == CONFLICT:
        # TBD: All these conflicting nodes need to turned into synonyms
        # of the final (rootward) unconflicting node.  Their nonconflicting
        # children should become children of the nonconflicting node.
        assert not is_top(q)    # don't go up from top
        q = get_superior(q)
    if related_how(p, q)[0] == GT:
      assert related_how(x, q)[0] == LT
      return (q, False)
    else:                     # LT EQ DISJOINT CONFLICT UNCLEAR
      return (p, True)

  for root in a_roots: traverse(root)
  if debug:
   print("# Finish: touched %s, set sup %s, capped %s, orphans %s, pass %s" % tuple(stats),
         file=sys.stderr)
  # end finish_sum


# j is to be either a child or synonym of k.  Figure out which.
# Invariant: A synonym must have neither children nor synonyms.

def set_superior(j, k):
  assert k
  assert not get_superior(j)
  assert not is_top(j)
  if j == k:
    print("!! self-loop %s" % get_blurb(j), file=sys.stderr)
    print("!! j = (%s, %s)" % (get_blurb(out_a(j)),
                               get_blurb(out_b(j))),
          file=sys.stderr)
    assert j != k
  k = get_accepted(k, k)
  if j == k:
    print("!! self-cycle %s" % get_blurb(j), file=sys.stderr)
    assert j != k
  x = out_a(j)
  y = out_b(j)
  # x and y might be synonym/accepted or accepted/synonym
  if ((y and get_accepted(y, None)) or
      (x and get_accepted(x, None))):
    # y a synonym; convert x from accepted to synonym, perhaps
    assert not get_children(j, None)
    assert not get_synonyms(j, None)
    set_accepted(j, k)
  else:
    set_parent(j, k)
  if monitor(x) or monitor(y):
    print("> superior of %s := %s" % (get_blurb(j), get_blurb(k)),
          file=sys.stderr)
          

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
    return (get_subblurb(r) or
            get_subblurb(out_b(r, None)) or
            get_subblurb(out_a(r, None)) or
            "[no canonical %s]" % get_key(r))
  else:
    return "[no match]"

def get_subblurb(r):
  if r:
    name = get_canonical(r, None) or get_scientific(r, None)
    if name:
      if get_accepted(r, None):
        return "*" + name
      else:
        return name
  else:
    return None

# -----------------------------------------------------------------------------
# 14. emit the new sum with additional parent column

# Returns a row generator

def generate_alignment(sum, roots, rm_sum):
  yield ["taxonID", "taxonID_A", "taxonID_B",
         "parentNameUsageID", "acceptedNameUsageID",
         "taxonomicStatus",
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
    change = MISSING
    if b_usage:
      z = out_a(mep_get(rm_in_b, b_usage), None)
      if z:
        (rcc5, _, _, _) = related_how(z, b_usage)
        change = rcc5_symbols[rcc5]
        m = mep_get(in_a, z, None)
      # if a_usage and b_usage then "renamed"
      # if b_usage then if get_xmrca(b_usage) then "collected"
    else:
      z = out_b(mep_get(rm_in_a, a_usage), None)
      if z:
        (rcc5, _, _, _) = related_how(a_usage, z)
        change = rcc5_symbols[rcc5]
        m = mep_get(in_b, z, None)
      # if a_usage then if get_xmrca(a_usage) then "dissolved"
    return [get_key(union),
            get_key(a_usage) if a_usage else MISSING,
            get_key(b_usage) if b_usage else MISSING,
            get_proper_key(p) if p else MISSING,
            get_proper_key(a) if a else MISSING,
            figure_taxonomic_status(union),
            get_canonical(union, MISSING),
            # m is record match
            get_key(m) if m else MISSING, # = previousID for name
            change,
            get_remark(union, MISSING)]
  def traverse(union):
    if get_key(union) != TOP:
      yield compose_row(union)
    for inf in get_inferiors(union):
      for row in traverse(inf): yield row
  for root in roots:
    for row in traverse(root): yield row

def figure_taxonomic_status(u):
  u_basis = out_a(u, None) or out_b(u, None)
  acc = get_accepted(u, None)
  if acc:
    return seniority(u_basis, acc)
  else:
    return "accepted"

def get_proper_key(u):
  if is_top(u): return MISSING
  else: return get_key(u)

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
