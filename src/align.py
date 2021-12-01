#!/usr/bin/env python3

# Align two trees (given as DwC files i.e. usage lists), with
# sensitivity to parent/child relations.

import sys, csv, argparse
import util
from util import windex, MISSING
import property as prop
from property import mep_get, mep_set
import match_records

TOP = "[top]"
debug = False

def is_top(x): return get_key(x) == TOP

def monitor(x):
  return (x and get_canonical(x, None).startswith("Metachirus"))

# -----------------------------------------------------------------------------
# Supervise the overall process

def align(a_iterator, b_iterator, rm_sum_iterator=None):
  if rm_sum_iterator == None:
    # Need to copy the iterators!
    (a_iterator, a_iterator_copy) = dup_iterator(a_iterator)
    (b_iterator, b_iterator_copy) = dup_iterator(b_iterator)
    rm_sum_iterator = match_records.match_records(a_iterator_copy, b_iterator_copy)
  # Load record matches
  (a_usage_dict, a_roots) = load_usages(a_iterator)
  (b_usage_dict, b_roots) = load_usages(b_iterator)
  rm_sum = load_sum(rm_sum_iterator, a_usage_dict, b_usage_dict)

  calculate_xmrcas(a_roots, b_roots, rm_sum)

  # Connect the two trees
  sum = set_equivalences(a_roots, b_roots, rm_sum)

  # Assign names
  assign_canonicals(sum)

  # Create a single merged hierarchy
  roots = set_superiors(a_roots, b_roots, sum, rm_sum)

  # Emit tabular version of merged tree
  return generate_alignment(sum, roots, rm_sum)

def dup_iterator(iter):
  clunk = list(iter)
  return ((x for x in clunk), (x for x in clunk))

def assign_canonicals(sum):
  (key_to_union, in_a, in_b) = sum
  b_index_by_name = {}
  for u in key_to_union.values():
    y = out_b(u)
    if y:
      name = get_canonical(y, None)
      if name:
        b_index_by_name[name] = y
  if debug:
    print("# %s B canonicals" % len(b_index_by_name),
          file=sys.stderr)
  count = sci_count = 0
  losers = 0
  for u in key_to_union.values():
    y = out_b(u)
    if y:
      name = get_canonical(y, None)
      sci_name = get_scientific(y, None)
    else:
      x = out_a(u)
      name = get_canonical(x, None)
      sci_name = get_scientific(x, None)
      if name in b_index_by_name:
        name = name + " sec. A"
    if name:
      set_canonical(u, name)
      count += 1
    if sci_name:
      set_scientific(u, sci_name)
      sci_count += 1
    if not name and not sci_name:
      losers += 1
  if debug:
    print("# %s canonicals, %s scientifics, %s nameless" %
          (count, sci_count, losers),
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

  seniors = 0
  for (usage, parent_key, accepted_key) in fixup:
    if accepted_key != MISSING:
      accepted_usage = key_to_usage.get(accepted_key)
      if accepted_usage:
        set_accepted(usage, accepted_usage)
        # Filter out senior synonyms here
        if seniority(usage, accepted_usage) == "senior synonym":
          seniors += 1
          del key_to_usage[get_key(usage)] # ????
        else:
          if monitor(usage) or monitor(accepted_usage):
            print("> accepted %s := %s" %
                  (get_blurb(usage), get_blurb(accepted_usage)),
                  file=sys.stderr)
      else:
        print("-- Dangling accepted: %s -> %s" % (get_key(usage), accepted_key),
              file=sys.stderr)
    elif parent_key != MISSING:
      parent_usage = key_to_usage.get(parent_key)
      if parent_usage:
        set_parent(usage, parent_usage)
        if monitor(usage) or monitor(parent_usage):
          print("> parent %s := %s" % (get_blurb(usage), get_blurb(parent_usage)),
                file=sys.stderr)
      else:
        print("-- Dangling parent: %s -> %s" % (get_key(usage), parent_key),
              file=sys.stderr)

  if seniors > 0:     # Maybe interesting
    print("Suppressed %s senior synonyms" % seniors,
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
  while get_level(x) < get_level(y):
    y = get_superior(y)
  while get_level(x) > get_level(y):
    x = get_superior(x)
  return (x, y)

def le(x, y):
  while get_level(x) > get_level(y):
    x = get_superior(x)
  return x == y

def lt(x, y): return le(x, y) and x != y

def mrca(x0, y0):
  (x, y) = find_peers(x0, y0)
  while not (x is y):
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
remarks_prop = prop.Property("remarks")
get_remarks = prop.getter(remarks_prop)
set_remarks = prop.setter(remarks_prop)

make_union = prop.constructor(key_prop, out_a_prop, out_b_prop,
                              remarks_prop)


# Load record mathes from a file or whatever

def load_sum(iterator, a_usage_dict, b_usage_dict):
  header = next(iterator)
  usage_a_pos = windex(header, "taxonID_A")
  usage_b_pos = windex(header, "taxonID_B")
  remarks_pos = windex(header, "remarks")

  sum = ({}, {}, {})
  for row in iterator:
    x = a_usage_dict.get(row[usage_a_pos])
    y = b_usage_dict.get(row[usage_b_pos])
    if x or y:
      note_match(x, y, row[remarks_pos], sum)

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
# This function applies to both record match sums and to alignments

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
            (get_key(x), get_key(out_b(j)), get_remarks(j),
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

# -----------------------------------------------------------------------------
# 8. Cross-mrcas

center_prop = prop.Property("center")
get_center = prop.getter(center_prop)
set_center = prop.setter(center_prop)

xmrca_prop = prop.Property("xmrca")
get_xmrca = prop.getter(xmrca_prop)
set_xmrca = prop.setter(xmrca_prop)

only_prop = prop.Property("only") # or 'monotype'
get_only = prop.getter(only_prop)
set_only = prop.setter(only_prop)

def calculate_xmrcas(a_roots, b_roots, rm_sum):
  calculate_centers(a_roots, b_roots, rm_sum)

def calculate_centers(a_roots, b_roots, rm_sum):
  (_, rm_in_a, rm_in_b) = rm_sum
  def half(roots, other, record_match):
    def traverse(x):
      only = None
      y = None                  # mrca of the central children
      for c in get_inferiors(x):
        y_c = traverse(c)         # in b
        if y_c:
          if y:
            y = mrca(y, y_c)
            only = None
          else:
            y = y_c
            only = y_c
      if not y:
        # Match tip to tip or to internal node (asymmetric center)
        y = record_match(x)
      if y:
        set_center(x, y)
        if monitor(x):
          print("> center(%s) := %s %s" % (get_blurb(x), get_blurb(y), not not only),
                file=sys.stderr)
        if only and not is_top(x):        
          assert y == only
          if monitor(x):
            print("> only(%s) := %s" % (get_blurb(x), get_blurb(only)),
                  file=sys.stderr)
          set_only(x, only)
      return y
    for root in roots:
      traverse(root)
      set_center(root, other)
  half(a_roots, b_roots[0], lambda x:out_b(mep_get(rm_in_a, x)))
  half(b_roots, a_roots[0], lambda y:out_a(mep_get(rm_in_b, y)))

  set_xmrcas(a_roots, b_roots, rm_sum)
  check_xmrcas(a_roots, b_roots)

# Use record matches to improve xmrcas within clusters

def set_xmrcas(a_roots, b_roots, rm_sum):
  (_, rm_in_a, rm_in_b) = rm_sum
  def half(roots, top, record_match):
    def traverse(x, y_upper):
      d = get_center(x, None)
      assert y_upper
      if is_top(x) or not d or is_top(d):
        y = top
      elif get_only(x, None):
        # Go just rootward of center
        y = get_superior(d)
        assert y

        # See if record match (even more rootward) can override
        r = record_match(x)    # on the other side
        if r:
          # Match has to be to another only taxon...
          # and it has to be rootward of x's center ('only')...
          # and it has to have the same center as x's center
          if (get_only(r, None) and le(y, r) and 
              get_center(r) == get_center(d)):
            if monitor(x):
              # assert lt(s, y_upper)  etc etc
              print("> Record match for %s" %
                    (get_blurb(x),),
                    file=sys.stderr)
            # Cleave the cluster!
            y = r
      else:
        # Not linear, so no change to xmrca, but propagate.
        y = d or y_upper
      assert not d or le(d, y)
      if d != y:
        if monitor(x):
          print("> %s: center %s -> xmrca %s, only=%s" %
                (get_blurb(x), get_blurb(d),
                 get_blurb(y), not not get_only(x, None)),
                file=sys.stderr)
      set_xmrca(x, y)
      for c in get_inferiors(x): traverse(c, y)
    for x in roots:
      traverse(x, top)

  def equivalent(x, y):
    return x == get_center(y) and get_center(x) == y

  half(a_roots, b_roots[0], lambda x:out_b(mep_get(rm_in_a, x)))
  half(b_roots, a_roots[0], lambda y:out_a(mep_get(rm_in_b, y)))

def check_xmrcas(a_roots, b_roots):
  def check(x):
    for c in get_inferiors(x): check(c)
    y = get_xmrca(x, None)
    if y and not get_xmrca(y, None):
      print("!! Non-mutual xmrca: x %s y %s" %
            (get_blurb(x), get_blurb(x)),
            file=sys.stderr)
      assert False
  for root in a_roots: check(root)
  for root in b_roots: check(root)

# -----------------------------------------------------------------------------
# Detect and record A/B equivalences.  At the end, every usage will be
# part of the sum.

def set_equivalences(a_roots, b_roots, rm_sum):
  sum = ({}, {}, {})
  (_, in_a, in_b) = sum
  (_, rm_in_a, rm_in_b) = rm_sum

  connect = connector(rm_sum, sum)

  def traverse(x):
    mem = False
    for c in get_inferiors(x):
      mem = traverse(c) or mem
    y = get_xmrca(x, None)      # in B
    if y != None and x == get_xmrca(y, None):
      # Mutual xmrca with xmrca... that means they're equivalent
      j = connect(x, y)

      # If there was a record match, preserve reason
      r = mep_get(rm_in_a, x, None)
      if r and out_b(r) == y:
        add_remark(j, get_remarks(r))

      # If internal node match, make a note of it
      if mem: add_remark(j, "consistent membership")
      return True
    else:
      connect(x, None)
      if monitor(x):
        print("> Not equivalent(%s, %s)" % (get_blurb(x), get_blurb(y)),
              file=sys.stderr)
      return False
  for x in a_roots: traverse(x)

  def cotraverse(y):
    for c in get_inferiors(y): cotraverse(c)
    if not mep_get(in_b, y, None):
      connect(None, y)
  for y in b_roots: cotraverse(y)

  return sum

def connector(rm_sum, sum):
  (_, rm_in_a, rm_in_b) = rm_sum
  def connect(x, y):
    return note_match(x, y, MISSING, sum)
  return connect

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
  assert p and q
  result = really_related_how(p, q)
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
  if d:
    p0 = get_xmrca(d, None)     # d <= p0
  else:
    p0 = p

  if monitor(p) or monitor(p0):
   if not le(p, p0):
    print("%s not<= %s <= %s" % (get_blurb(p), get_blurb(d), get_blurb(p0)),
          file=sys.stderr)

  c = get_xmrca(q, None)        # q ? c
  if c:
    q0 = get_xmrca(c, None)     # c <= q0
  else:
    q0 = q

  if not c and not d:
    return (DISJOINT, "both peripheral", None, None)

  # We can now test relative order of p, c, p0 and q, d, q0
  # subject to d <= p0 and c <= q0

  # p ? c <= q0 ? q      and vice versa, p ? p0 >= d ? q

  p0_le_p = le(p0, p)
  q0_le_q = le(q0, q)
  p_lesseq_q = (c and le(p, c) and q0_le_q)
  q_lesseq_p = (q and le(q, d) and p0_le_p)

  if p_lesseq_q and q_lesseq_p:
    return (EQ, "p <= c <= q0 <= q <= d <= p0 <= p", None, None)
  if p_lesseq_q and (p != c or q0 != q):
    return (LT, "p <= c <= q0 <= q", None, None)
  if q_lesseq_p and (q != d or p0 != p):
    return (GT, "q <= d <= p0 <= p", None, None)

  # Hmm, tricky.  A peripheral node is either disjoint from or less
  # than each node in the other tree.  Both cases are handled above.
  # Therefore neither p nor q is peripheral at this point.
  assert c and d

  # Brute force
  (f1, e) = seek_conflict(p, q)
  (f2, g) = seek_conflict(q, p)

  if f1 or f2:
    if e:
      if g:
        return (CONFLICT, e, f1 or f2, g)
      else:
        return (GT, "checked", None, True)
    elif g:
      return (LT, "checked", True, None)
    else:
      # Not enough information to be able to decide
      if monitor(p) or monitor(q):      
        print("!! Unclear - %s ? %s\n   %s, %s" %
              (get_blurb(p), get_blurb(q), get_blurb(f1), get_blurb(f2)),
              file=sys.stderr)
      return (UNCLEAR, "checked ??", True, True)
  elif e or g:
    return (DISJOINT, "checked", p, q)
  else:                         # Within a cluster
    if monitor(p) or monitor(q):
      print("!! Unclear. %s ? %s" % (get_blurb(p), get_blurb(q)),
            file=sys.stderr)
    return (UNCLEAR, "checked", None, None)

# Returns two sibling(?) descendants of p
# Assumes initially that either p > q or p >< q

def seek_conflict(p, q):
  o = get_center(p, None)
  if o == None:
    # Peripheral - no way to know
    return (None, None) #  ????????????
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
    return (p, p)     # (None, None)
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
        # return (None, None)  # was: return (None, p_not_q)
        return (p_and_q, p_not_q)
    return (p_and_q, p_not_q)

# -----------------------------------------------------------------------------
# Now set the parent pointers for the spanning tree

def set_superiors(a_roots, b_roots, sum, rm_sum):
  (key_to_union, in_a, in_b) = sum
  (_, rm_in_a, rm_in_b) = rm_sum
  connect = connector(rm_sum, sum)

  half_set_superiors(b_roots, in_b, in_a, True)    # first, B = priority
  half_set_superiors(a_roots, in_a, in_b, False)
  if debug:
    print("# %s in a, %s in b, %s union keys" %
          (len(in_a), len(in_b), len(key_to_union)),
          file=sys.stderr)

  unions = key_to_union.values()
  roots = collect_inferiors(unions)
  print("-- align: %s usages, %s roots in sum" % (len(unions), len(roots)),
        file=sys.stderr)
  return roots

# Set parents in the 'a' hierarchy, which might be high priority (B)
# or low priority (A).

def half_set_superiors(a_roots, in_a, in_b, priority):

  elisions = {}

  def traverse(x):
    choose_superior(x, priority)
    for c in get_inferiors(x): traverse(c)

  def choose_superior(x, priority):
    if is_top(x): return None
    p = get_superior(x)         # p is in a (possibly top)
    u = mep_get(in_a, x)
    answer = get_superior(u)
    if not answer:
      q = get_xmrca(x, None)      # q is in b
      if q == None:
        answer = mep_get(in_a, p)    # peripheral
      else:
        rel = related_how(x, q)[0]
        assert rel != GT
        if rel == EQ:    # GT case should never happen
          q = get_superior(q)
        answer = choose_nearest(p, q, priority)
      set_superior(u, answer)
    return answer

  # Pick the most rootward of the two

  def choose_nearest(p, q, priority):
    assert p and q
    result = related_how(p, q)
    rel = result[0]
    if rel == LT or rel == EQ:
      answer = mep_get(in_a, p)
    elif rel == GT:
      answer = mep_get(in_b, q)
    elif priority:
      answer = choose_nearest(p, get_superior(q), priority)
    else:
      # a = A = low priority
      answer = choose_nearest(get_superior(p), q, priority)
      mep_set(elisions, p, (p, result, q, answer))
    return answer

  for root in a_roots: traverse(root)

  # Report on elisions (well wait, what are these semantically?)

  def report_on_elisions(elisions):
    print("* %s elisions.  Report follows." % len(elisions),
          file=sys.stderr)
    yield(["p in A", "rcc5", "q in B", "⊆ p-q", "⊆ p∩q",
           "⊆ q-p", "p∨q", "multi-parent"])
    for (p_elide, result, q, safe) in elisions.values():
      # Need to remove q_eject from the merged hierarchy.  Turn it into a synonym
      # of safe and move its inferiors there.
      eliding = mep_get(in_a, p_elide)
      (rel2, g, f, e) = result
      set_parent(eliding, None)   # Detach it from the A hierarchy
      # Move the p-side children to a safe place.
      degraded = 0
      for c in get_inferiors(p_elide):
        j = mep_get(in_a, c)
        if get_superior(j) == eliding:
          # set_alt_parent(j, eliding)
          change_superior(j, safe)
          if get_accepted(j, False):
            add_remark(j, "synonym of elided %s" % get_blurb(eliding))
          else:
            add_remark(j, "child of elided %s" % get_blurb(eliding))
          degraded += 1

      # Demote the ejecting node to a synonym
      add_remark(eliding, ("elided because %s %s; parent was %s" %
                           (rcc5_symbols[rel2], get_blurb(q),
                            get_blurb(get_superior(p_elide)))))
      # set_alt_parent(eliding, get_parent(eliding))
      set_accepted(eliding, safe)

      if rel2 == CONFLICT:
        writer.writerow((get_blurb(p_elide), rcc5_symbols[rel2], get_blurb(q),
                         get_blurb(e), get_blurb(f), get_blurb(g),
                         get_blurb(safe), str(degraded)))
      else:   # UNCLEAR or DISJOINT
        yield((get_blurb(p_elide), rcc5_symbols[rel2], get_blurb(q),
               MISSING, MISSING, MISSING,
               get_blurb(safe), str(degraded)))

  if len(elisions) > 0:
    writer = csv.writer(sys.stderr)
    for row in report_on_elisions(elisions):
      writer.writerow(row)

def add_remark(x, rem):
  have = get_remarks(x, None)
  if have and have != MISSING:
    set_remarks(x, have + '|' + rem)
  else:
    set_remarks(x, rem)

# j is to be either a child or synonym of k.  Figure out which.
# Invariant: A synonym must have neither children nor synonyms.

def set_superior(j, k):
  assert not get_superior(j)
  change_superior(j, k)

def change_superior(j, k):
  assert k
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
          rcc5 = related_how(x, y)[0]
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
  if isinstance(r, dict):
    s = get_subblurb(r)
    if s: return s
    s1 = get_subblurb(out_b(r, None))
    s2 = get_subblurb(out_a(r, None))
    if s1 and (s1 == s2): return s1
    if s1: return "B." + s1
    if s2: return "A." + s2
    return "[no name %s]" % get_key(r)
  elif r:
    return "[not a dict]"
  else:
    return "[no match]"

def get_subblurb(r):
  if r:
    name = get_canonical(r, None) or get_scientific(r, None)
    if name:
      if get_accepted(r, None):
        return name + "*"     # attend to sort order
      else:
        return name
  return None

# -----------------------------------------------------------------------------
# 14. emit the new sum with additional parent column

# Returns a row generator

def generate_alignment(sum, roots, rm_sum):
  yield ["taxonID", "taxonID_A", "taxonID_B",
         "parentNameUsageID", "acceptedNameUsageID",
         "taxonomicStatus",
         "canonicalName", "scientificName", "recordMatch", "change", "remarks"]
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
        rcc5 = related_how(z, b_usage)[0]
        change = rcc5_symbols[rcc5]
        m = mep_get(in_a, z, None)
      # if a_usage and b_usage then "renamed"
      # if b_usage then if get_xmrca(b_usage) then "collected"
    else:
      z = out_b(mep_get(rm_in_a, a_usage), None)
      if z:
        rcc5 = related_how(a_usage, z)[0]
        change = rcc5_symbols[rcc5]
        m = out_a(mep_get(in_b, z, None))
      # if a_usage then if get_xmrca(a_usage) then "dissolved"
    return [get_key(union),
            get_key(a_usage) if a_usage else MISSING,
            get_key(b_usage) if b_usage else MISSING,
            get_proper_key(p) if p else MISSING,
            get_proper_key(a) if a else MISSING,
            figure_taxonomic_status(union),
            get_canonical(union, MISSING),
            get_scientific(union, MISSING),
            # m is record match
            get_key(m) if m else MISSING, # = recordMatch for name
            change,
            get_remarks(union, MISSING)]
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

  a_iter = csv.reader(a_file)
  with open(b_path) as b_file:
    b_iter = csv.reader(b_file)
    # TBD: Compute sum if not provided ?
    if rm_sum_path:
      with open(rm_sum_path) as rm_sum_file:
        rm_sum_iter = csv.reader(rm_sum_file)
        sum = align(a_iter, b_iter, rm_sum_iter)
    else:
      sum = align(a_iter, b_iter)
    util.write_generated(sum, sys.stdout)
