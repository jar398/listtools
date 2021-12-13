#!/usr/bin/env python3

NEW_STYLE = True

# Align two trees (given as DwC files i.e. usage lists), with
# sensitivity to parent/child relations.

import sys, csv, argparse
import util
from util import windex, MISSING, log
import match_records
from rcc5 import *

import property as prop
from property import mep, mep_get, mep_set
from property import get_property

import checklist
from checklist import Source, get_checklist

TOP = "[top]"
debug = False

def is_top(x): return get_key(x) == TOP

alt_parent_prop = prop.get_property("alt_parent") # Next bigger
get_alt_parent = prop.getter(alt_parent_prop)
set_alt_parent = prop.setter(alt_parent_prop)

# -----------------------------------------------------------------------------
# Supervise the overall process

def align(a_iterator, b_iterator, rm_sum_iterator=None):
  global rm_sum  # foo
  if rm_sum_iterator == None:
    # Need to copy the iterators!
    (a_iterator, a_iterator_copy) = dup_iterator(a_iterator)
    (b_iterator, b_iterator_copy) = dup_iterator(b_iterator)
    rm_sum_iterator = match_records.match_records(a_iterator_copy, b_iterator_copy)
  A = checklist.load_source(a_iterator, "A")
  B = checklist.load_source(b_iterator, "B")
  # Load record matches
  rm_sum = load_sum(rm_sum_iterator, A, B)

  sum = checklist.make_sum(A, B, "", rm_sum=rm_sum)
  forge_links(sum, rm_sum)   # Find links to other checklist for all nodes (lower bounds)

  # Create sum and connect the two trees
  analyze_order(sum, rm_sum)

  # Create a single merged hierarchy
  roots = spanning_tree(A, B, sum)

  # Emit tabular version of merged tree
  return emit_spanning_tree(sum, roots, rm_sum)

def dup_iterator(iter):
  clunk = list(iter)
  return ((x for x in clunk), (x for x in clunk))


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
# MRCA within the same checklist

# E.g. level(x) >= level(y)  should imply  x <= y  if not disjoint

def le(x, y):
  y1 = x     # proceeds down to y
  # if {level(x) >= level(y1) >= level(y)}, then x <= y1 <= y or disjoint
  stop = get_level(y)
  while get_level(y1) > stop:
    y1 = get_superior(y1)    # y1 > previously, level(y1) < previously
  # level(y) = level(y1) <= level(x)
  return y1 == y    # Not > and not disjoint

def lt(x, y): return le(x, y) and x != y

def find_peers(x, y):
  while get_level(x) < get_level(y):
    y = get_superior(y)
  while get_level(x) > get_level(y):
    x = get_superior(x)
  return (x, y)

def mrca(x, y):
  if x == BOTTOM: return y
  if y == BOTTOM: return x
  (x, y) = find_peers(x, y)
  while not (x is y):
    x = get_superior(x)
    y = get_superior(y)
  return x

def relation(x, y):             # Within a single tree
  (x1, y1) = find_peers(x, y)    # Decrease levels as needed
  assert not x1.checklist is y1.checklist
  if x1 == y1:
    if x == y:
      return EQ
    elif x1 == x:
      return LT     # y > y1 = x
    elif y1 == y:
      return GT     # x > x1 = y
    else:
      assert False
  else:
    x2 = get_accepted(x1, x1)
    y2 = get_accepted(y1, x1)
    if x2 == y2:                                  # Same accepted
      if x1 != x2 and y1 != y2: return NOINFO    # synonym ? synonym
      if x1 != x2: return LE                      # synonym <= accepted
      else:        return GE                      # accepted >= synonym
    else:
      return DISJOINT

# -----------------------------------------------------------------------------
# Sum file ingest

# Record-match file ingest

remarks_prop = prop.get_property("remarks")
get_remarks = prop.getter(remarks_prop)
set_remarks = prop.setter(remarks_prop)


# Load record mathes from a file or whatever

def load_sum(iterator, A, B):
  header = next(iterator)
  usage_a_pos = windex(header, "taxonID_A")
  usage_b_pos = windex(header, "taxonID_B")
  remarks_pos = windex(header, "remarks")

  rm_sum = checklist.make_sum(B, A, "match")
  for row in iterator:
    x = A.index.by_key_dict.get(row[usage_a_pos])
    y = B.index.by_key_dict.get(row[usage_b_pos])
    if x or y:
      rm_sum.join(x, y, row[remarks_pos])
    # To do: save other fields ??

  return rm_sum

# Sadly, rows are repeated sometimes 
# This function applies to both record match sums and to alignments

def really_connect(sum, x, y, remark):
  global union_count
  if x and y:
    assert is_top(x) == is_top(y), \
      "connecting %s to %s" % (get_blurb(x), get_blurb(y))
  j1 = sum.in_left(x, None) if x else None
  j2 = sum.in_right(y, None) if y else None
  if j1 and j2: assert j1 == j2
  j = j1 or j2
  if j:
    assert sum.out_left(j) == x and sum.out_right(j) == y, \
      ("want <%s,%s>\n  but already have <%s,%s>  ;%s" %
       (get_unique(x), get_unique(y), 
        get_unique(sum.out_left(j)), get_unique(sum.out_right(j)),
        get_remarks(j)))
    if x: assert sum.in_left(x) == j
    if y: assert sum.in_right(y) == j
    return j
  # Invent a key?? and store it in the sum ???
  # To avoid an A/B conflict we'd need the B key->usage map?
  if y:
    key = get_key(y)
    if key in sum.index.by_key_dict: key = "B.%s" % key
  else:
    key = get_key(x)
    if key in sum.index.by_key_dict: key = "A.%s" % key
  
  if y:
    if x:
      assert get_unique(x, None), get_blurb(x)
      assert get_unique(y, None), get_blurb(y)
      (p, q) = (get_unique(x), get_unique(y))
      if p == q:
        unique = "<%s>" % p
      else:
        unique = "<%s,%s>" % (p, q)
    else:
      unique = "<,%s>" % pick_unique(y)
  else: unique = "<%s,>" % pick_unique(x)

  j = make_union(x, y, sum)
  if x:
    mep_set(sum.in_left_dict, x, j)
    mep_set(sum.out_left_dict, j, x)
    set_out_A(j, x)
    add_sourced_remark(j, x)
  if y:
    mep_set(sum.in_right_dict, y, j)
    mep_set(sum.out_right_dict, j, y)
    set_out_B(j, y)
    add_sourced_remark(j, y)
  sum.index.by_key_dict[key] = j
  sum.index.by_unique_dict[unique] = j
  union_count += 1

  if monitor(x) or monitor(y):
    print("> note: connecting %s to %s; %s" %
          (get_blurb(x), get_blurb(y), remark),
          file=sys.stderr)

  return j

union_count = 0

def add_sourced_remark(j, z):
  rem = get_remarks(z, None)
  if rem:
    add_remark(j, "%s:%s" % (z.source.name, rem))

# -----------------------------------------------------------------------------
# Links between checklists

# Explanation of 'linkage' (lower bound).
# We're seeking an interpretation [·] of A+B that is consistent with
# whatever the interpretations A and B are.
# To do this we need to construct a child/parent relationship in A+B
# satisfying:
#   parent[x] = the least z in A+B such that z > x (i.e. [x] ⊂ [z]).
# This is done bottom up: z is replaced by its parent for as long as
# [x] ⊄ [z].
# At each stage we know the relation between z and x, although it
# may not be an RCC5 relation.  The pair (relation, candidate) is
# called x's 'link'.  Initially the link is just (LINKED y)
# where y is computed by bounds analysis.

# The 'link' is always ? to z, and is cranked up to the parent if
# necessary (so that we can capture EQ relationships along the way).

# The purpose of LINKED is to speed up comparisons in the
# inner loop of cross_not_le - you never need to descend below the
# link, because we know it's contained entirely in both x and y.

link_prop = get_property("link")
get_link = prop.getter(link_prop)
set_link = prop.setter(link_prop)

BOTTOM = None

def half_forge_links(AB, rm_AB):

  def traverse(u):              # u in A
    for c in A_inferiors(u): traverse(c)
    ry = calculate_link(c_j)
    if ry:
      (rel, y_j) = ry
      if monitor(u):
        log("> link(%s) := (%s, %s)" %
            (get_blurb(u), rcc5_symbol(rel), get_blurb(y)))

      set_bound(x, ry)
    return ry

  def calculate_link(x, AB):   # x in AB
    # y = functools.reduce(mrca, A_inferiors(u), BOTTOM)
    y = BOTTOM                  # mrca of opposite mrcas
    for c in A_inferiors(u):
      (y_c, more_c) = traverse(c)
      y = mrca(y, y_c)

    if y != BOTTOM:
      return (LINKED, y)
    else:
      r = record_match(y)
      # Necessarily [x] <= [r], so link is (LE, r)
      return (EQ, AB.in_left(r)) if r else (BOTTOM, NOINFO)    # hmm.  helps caching

  def record_match(u):          # u in A+B
    y = out_left(u)             # y in A
    r = out_right(rm_AB.in_left(y))
    return AB.in_right(r) if r else None

  for root in AB.A.roots: traverse(root)

def forge_links(AB, rm_AB):
  half_forge_links(AB.flip(), rm_AB.flip())
  half_forge_links(AB, rm_AB)

def A_inferiors(x, AB):
  return (AB.inject(c) for c in get_inferiors(out_left(x)))

def get_xmrca(x, AB):           # compatibility kludge
  return link_to_xmrca(get_bound(x), AB)
def link_to_xmrca(ry, AB):
  (rel, y) = ry                # x rel y
  if y == BOTTOM: return None
  if rel == EQ: return (y, 0)
  else: return (y, -1)

# -----------------------------------------------------------------------------
# Detect and record A/B equivalences

#  p  ?  q  ?  x0        - parents
#  ∨     ∨     ∥
#  x  →  y  →  x0  ?  x

def half_analyze_order(sum, rm_sum):
  assert sum.A is rm_sum.A and sum.B is rm_sum.B

  sum.rm_sum = sum              # ??

  # Start with [x] ⊃ [y]

  # Eventual link y_final will satisfy y rel y_final
  # rel starts out being LINKED, or LE for tips

  # Find smallest z such that x < z

  def peer_or_superior(x):
    ry = get_link(x)            # Start here and go rootward
    (rel, y) = ry
    if rel & DISJOINT != 0: return ry
    while True:
      z = peer_or_superior_candidate(x, y)
      if is_peer_or_superior(z): break
      y = get_superior_same_side(z)
    assert sum_le(x, y)
    return ry

    # Different clusters
      # y is in a bigger cluster - do nothing, keep going rootward
      # return (BOTTOM, NOINFO)

  def sum_superior(x):
    ry = peer_or_superior(x)    # LT, LE, or EQ
    if ry[0] == LT: return 
    if ry[0] == EQ: return 
    assert ry[0] == LE

  # Assumes ry is a peer or superior candidate
  def is_peer_or_superior(ry):
    (rel, y) = ry
    return rel & ~GE == 0

  # Find least z among {p, y, q} subject to x <= z

  def superior_candidate(x, y):
    p = get_superior_same_side(x)
    q = get_superior_same_side(y)
    # Consider candidates p, y, and q
    def good_enough(z):     # winner must be smallest of the three
      return (sum_le(x, z) and
              sum_le(z, p) and
              sum_le(z, y) and
              sum_le(z, q))
    z = None
    if   good_enough(p): z = p
    elif good_enough(y): z = y
    elif good_enough(q): z = q
    else: return get_bound(p)    # ?
    return (LE, z)

  def brute_eq(x, y):
    return sum_ge(x, y) and sum_ge(y, x)

  def in_same_cluster(y, x):    # is y in same cluster as x
    (rel2, x0) = sum.get_bound(y)
    return le(x0, x)

  # TDB: Write a really simple one that works
  def brute_propagate(x):
    for c in get_inferiors(x):
      super = peer_or_superior(c)
    if super: set_superior(x, super)

  def propagate(x, y = None):
    if y == None:
      (rel1, y) =  sum.get_bound(x) # almost always LINKED, LE at tips
      # or, y = parent(d) ??
    (rel2, x0) = sum.get_bound(y)
    if le(x0, x):
      # SAME CLUSTER - maybe room for improvement
      r = record_match(x)
      if rel2 == LE:
        # Tip matching something - take least option
        link = (EQ, y)
      elif get_mono(x, None) == y0 or get_superior_same_side(x) == y0:
        # Ambiguous on the A side
        ...
      else:
        # Internal to internal or to tip
        q = get_superior_same_side(y)
        q0 = sum.get_bound(q)[1]
        if q0 == x0:
          # AMBIGUOUS - look for record match
          r = record_match(x)
          if r == y:
            bound = (EQ, y) # Success!
          elif r == q:
            connect(None, y)    # Skip this one
            y = q               # Go toward root
          else:
            # Choice between y and q must be resolved by brute force, I think
            # q to y could still be any of < > = ><.  Not = I suppose.
            # Just call it a reciprocal conflict though - this is a race
            # since the outcome will depend on which polarity is processed
            # first?
            # CONFLICT if cross_not_le(x, y, sum) else LT
            bound = (UNRESOLVED, y)
        else:
          # Unambiguous in B, answer is y by process of elimination
          bound = (EQ, y)
    else:
      # Different clusters - x < y - y is parent of x, if x's parent isn't
      bound = (LT, y)

    if lt(x, x0): return False
    # Ambiguous iff child or parent is in same cluster.
    return get_mono(x, None) == y0 or get_superior_same_side(x) == y0

  def traverse(x):
    # TBD: sum.connect(x, None)

    mem = False
    for c in get_inferiors(x):
      mem = traverse(c) or mem
    y = get_equivalent(x, AB)
    if y != None:
      j = sum.join(y, x)

      # If there was a record match, preserve reason
      r = record_match(x)
      if r and r == y:
        add_remark(j, get_remarks(r, None))
      elif mem:
        # If internal node match, make a note of it
        add_remark(j, "consistent membership")
      mem = True
    else:
      sum.join(None, x)
      if monitor(x):
        log("> %s has no equivalent in %s" %
            (get_blurb(x), sum.A.name))
    assert sum.inject(x)
    return mem

  def record_match(x):
    r = rm_sum.in_right(x)
    return sum.out_left(r) if r else None

  for x in get_inferiors(sum.B.top): traverse(x)

def analyze_order(sum, rm_sum):
  sum.top = sum.join(sum.A.top, sum.B.top)
  half_analyze_order(sum, rm_sum)
  half_analyze_order(sum.flip(), rm_sum.flip())

  # Assign names - for Newick - too early
  assign_canonicals(sum)

  return sum

def sum_eq(x, y):
  if x == y: return True
  if get_checklist(x) == get_checklist(y): return False
  else:
    assert get_checklist(x) == get_checklist(y).flip()
    (x_rel, y2) == get_bound(x)
    (y_rel, x2) == get_bound(x)
    return ((x_rel == EQ and y == y2) or
            (y_rel == EQ and x == x2) or
            quick_eq(x, y))

# wait, this is not exact

def sum_le(x, y):
  AB = get_checklist(x)
  u = out_left(x)
  v = out_left(y)
  if same_side(x, y):
    return le(u, v)
  else:
    (rel, y2) == get_bound(x)
    # BLAH BLAH
    if ((rel == EQ or rel == LE) and le(out_left(y2), y)):
      return True
    return quick_le(u, v, AB)

def same_side(x, y):
  if get_checklist(x) == get_checklist(y):
    return True # eg AB==AB
  else:
    assert get_checklist(x).flip() == get_checklist(y)
    return False

# -----------------------------------------------------------------------------
# Some approximate methods for cross-checklist relationships

def quick_eq(x, y):
  return (quick_le(x, y) and
          quick_le(y, x))

# Find y, if any, in other checklist, such that x = y

def get_equivalent(x, AB):
  y = quick_bound(x, AB)
  return y if x is quick_bound(y, AB) else None

def quick_le(x, y, AB):
  y1 = quick_bound(x, AB)
  return y1 and le(y1, y)

# Find least y in other checklist such that x <= y, if any

def quick_bound(x, AB):
  # get_proper(get_xmrca(x))
  yt = get_xmrca(x, AB) or get_xmrca(attachment_point(x, AB), AB)
  assert yt
  y1 = get_proper(yt, AB)
  if monitor(x):
    log("> link for %s is %s" % (get_blurb(x), get_blurb(y1)))
  return y1

def attachment_point(x, AB):        # For peripheral nodes
  x1 = get_superior_same_side(x)
  if get_xmrca(x1, AB):
    return x1
  else:
    return attachment_point(x1)

# Most rootward node in SAME checklist having same shared-tipward set as x,
# i.e. an ancestor that is guaranteed to be non-conflicting.
# This gets us to the 'top' of 'monotypic' chains.
# Cf. find_xmrca(x), above, which goes to the 'bottom' of the chain.
# (conflict pairs can also form chains ... think about this.
#   p1 < p2 < p3 < p4, q1 < q2 < q3 < q4, 
#   p1 >< q1, p2 >< q2, p3 ≈ q3, p4 ≈ q4 )

def get_proper(xt, AB):      # Maybe cache this?
  if xt == None: return None    # foo
  # while superior(x1) == xmrca(x): x1 = superior(x1)
  (x, tweakx) = xt
  assert x
  if tweakx > 0:
    x2 = get_superior_same_side(x)
    log("> using tweak %s to force %s to superior %s" %
        (tweakx, get_blurb(x), get_blurb(x2)))
    if x2: x = x2
  yt = get_xmrca(x, AB)    # proper must have same xmrca
  if yt:
    (y0, tweak0) = yt
    x1 = x
    while x1:
      if is_top(x1): break
      x2 = get_superior(x1)
      x2t = get_xmrca(x2, AB)
      if not x2t: break
      (y2, tweak2) = x2t
      if y2 != y0: break    # want: y2 >= y0
      # tweak2 < tweak0 means tweak2 is rootward of tweak0 (like levels)
      if tweak2 < tweak0:
        log("> using tweaks to rule out %s as more proper for %s (@ %s)" %
             (get_blurb(x2), get_blurb(x), get_blurb(y0)))
        break
      if monitor(x):
        log("> proper to superior: %s -> %s, %s -> %s" %
            (get_blurb(x1), get_blurb(y0), get_blurb(x2), get_blurb(y2)))
      x1 = x2
    if monitor(x):
      log("> proper(%s) = %s" % (get_blurb(x), get_blurb(x1)))
    return x1
  else:  
    return x

# Estimate RCC5 relationship ACROSS trees

def quick_relation(x, y, AB):
  y1 = quick_bound(x, AB)
  x1 = quick_bound(y, AB)
  if le(y1, y):
    return EQ if le(x1, x) else LT
  if le(x1, x):
    return GT
  rel1 = relation(x1, x)
  rel2 = relation(y1, y)
  return rel1 if rel1 == reverse_relation(rel2) else CONFLICT

# -----------------------------------------------------------------------------
# Create a spanning tree by setting parent pointers.

def spanning_tree(A, B, sum):
  # First time, B is high priority, A is low (per usual)
  # Second time, B is low priority, A is high (reversed)
  half_spanning_tree(sum)
  half_spanning_tree(sum.flip())
  if debug:
    print("# %s in a, %s in b, %s union keys" %
          (len(sum.in_left_dict), len(sum.in_b_dict), len(sum.index.by_key_dict)),
          file=sys.stderr)

  unions = sum.index.by_key_dict.values()
  roots = collect_inferiors(unions)
  print("-- align: %s usages, %s roots in sum" % (len(unions), len(roots)),
        file=sys.stderr)
  return roots

# Set parents in the 'a' hierarchy, which might be high priority (B)
# or low priority (A).

def half_spanning_tree(sum):
  a = sum.A
  b = sum.B
  b_priority = sum.b_priority
  elisions = {}

  # For each node a, set a's superior's spanning tree node, if unset.

  def traverse(x):
    j = sum.in_left(x)
    if not get_superior(j):
      sup = spanning_superior(j)
      if sup:
        assert get_checklist(j) is get_checklist(sup), \
          "checklists differ: %s %s" % (get_blurb(j), get_blurb(sup))
        set_superior(j, sup, sum)
    for c in get_inferiors(x): traverse(c)

  def spanning_superior(j):
    x = sum.out_left(j)
    y = sum.out_right(j)
    if x: p = get_superior(x)
    if y: q = get_superior(y)
    if not x or not p:
      return sum.in_right(q) if q else None
    if not y or not q:
      return sum.in_left(p) if p else None
    return cross_min(p, q, sum)

  def cross_min(p, q, sum):
    if not p: return sum.in_right(q) if q else None
    if not q: return sum.in_left(p)
    if cross_le(p, q, sum): return sum.in_left(p)
    if cross_le(q, p, sum): return sum.in_right(q)
    if b_priority:
      p1 = get_superior(p)      # Skip p, which conflicts
      return cross_min(p1, q, sum) if p1 else sum.in_right(q)
    else:
      q1 = get_superior(q)
      return cross_min(p, q1, sum) if q1 else sum.in_left(p)

  for root in sum.A.roots: traverse(root)

  if len(elisions) > 0:
    writer = csv.writer(sys.stderr)
    for row in report_on_elisions(elisions, sum):
      writer.writerow(row)

# Report on elisions (well wait, what are these semantically?)

def report_on_elisions(elisions, sum):
  print("* %s elisions.  Report follows." % len(elisions),
        file=sys.stderr)
  yield(["p in A", "rcc5", "q in B", "⊆ p-q", "⊆ p∩q #1", "⊆ p∩q #2",
         "⊆ q-p", "p∨q", "multi-parent"])
  for (p_elide, q, safe) in elisions.values():
    safe = get_superior(sum.in_right(q))
    # Need to remove q_eject from the spanning hierarchy.  Turn it into a synonym
    # of safe and move its inferiors there.
    eliding = sum.in_left(p_elide)
    (rel2, e, f1, f2, g) = get_conflict_proof(p, q, sum)
    set_parent(eliding, None)   # Detach it from the A hierarchy
    # Move the p-side children to a safe place.
    degraded = 0
    for c in get_inferiors(p_elide):
      j = sum.in_left(c)
      if get_superior(j) == eliding:
        set_alt_parent(j, eliding)
        change_superior(j, safe)
        if get_accepted(j, False):
          add_remark(j, "synonym of elided %s" % get_blurb(eliding))
        else:
          add_remark(j, "child of elided %s" % get_blurb(eliding))
        degraded += 1

    # Demote the ejecting node to a synonym
    add_remark(eliding, ("elided because %s %s; parent was %s" %
                         (rcc5_symbol(rel2), get_blurb(q),
                          get_blurb(get_superior(p_elide)))))
    set_alt_parent(eliding, get_parent(eliding))
    set_accepted(eliding, safe)

    if rel2 == CONFLICT:
      writer.writerow((get_blurb(p_elide), rcc5_symbol(rel2), get_blurb(q),
                       get_blurb(e), get_blurb(f1), get_blurb(f1), get_blurb(g),
                       get_blurb(safe), str(degraded)))
    else:   # NOINFO or DISJOINT
      yield((get_blurb(p_elide), rcc5_symbol(rel2), get_blurb(q),
             MISSING, MISSING, MISSING,
             get_blurb(safe), str(degraded)))

def add_remark(x, rem):
  if rem:
    have = get_remarks(x, None)
    if have and have != MISSING:
      set_remarks(x, have + '|' + rem)
    else:
      set_remarks(x, rem)

# j is to be either a child or synonym of k.  Figure out which.
# Invariant: A synonym must have neither children nor synonyms.

def set_superior(j, k, sum):
  assert k
  assert get_checklist(j) is get_checklist(k)
  assert not is_top(j)
  if j == k:
    print("!! self-loop %s" % get_blurb(j), file=sys.stderr)
    print("!! j = (%s, %s)" % (get_blurb(sum.out_left(j)),
                               get_blurb(sum.out_right(j))),
          file=sys.stderr)
    assert j != k
  k = get_accepted(k, k)
  if j == k:
    print("!! self-cycle %s" % get_blurb(j), file=sys.stderr)
    assert j != k
  x = sum.out_left(j)
  y = sum.out_right(j)
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

# Complete implementation of RCC5 decisions

# Now, exact methods for cross-checklist relationships

def cross_relation(x, y, AB):
  return quick_relation(x, y, AB) or brute_relation(x, y, AB)
  
# Exact computation of relationship across trees
# Assumes they're neither = nor !

def brute_relation(x, y, AB):
  if cross_not_le(x, y, AB):
    # relation is GT, CONFLICT, or DISJOINT
    if cross_not_le(y, x, AB):
      # CONFLICT or DISJOINT
      return CONFLICT   # DISJOINT should have been caught by quick_relation
    else:
      # LT or EQ
      return LT                 # EQ is caught by quick_relation
  elif cross_not_le(y, x, AB):
    # CONFLICT or DISJOINT
    return CONFLICT

def cross_le(x, y, AB):
  return cross_not_le(y, x, AB) and not cross_not_le(x, y, AB)

# Return descendant of p that is not in q, if any

def cross_not_le(p, q, AB):
  def search(p1, q):
    if quick_le(p1, q, AB):
      return None
    else:
      infs = get_inferiors(p1)
      if len(infs) > 0:
        for c in infs:
          if search(c, q):
            return c
        return None
      else:
        return p1    #Tip not found by quick_le ... ?
  c = search(p, q)
  if monitor(c) or monitor(p) or monitor(q):
    if c:
      log("> %s seems to be <= %s" % (get_blurb(p), get_blurb(q)))
    else:
      log("> %s proves that %s is not <= %s" %
          (get_blurb(c), get_blurb(p), get_blurb(q)))

  return c

# For debugging

def get_conflict_proof(p, q, AB):
  def le_or_conflict_proof(p, q):
    qt = get_xmrca(p, AB)
    if not qt: return (None, None)
    (tweak, q1) = qt
    rel = relation(q1, q)
    if rel == LT:
      return (None, p)          # p < q
    elif rel == DISJOINT:
      return (p, None)          # p ! q
    e = f = None
    for c in get_children(p, []):
      (ec, fc) = get_conflict_proof(x, q)
      if e == None or get_accepted(ec): e = ec
      if f == None or get_accepted(fc): f = fc
      return (e, f)
  (e, f1) = le_or_conflict_proof(p, q)
  (g, f2) = le_or_conflict_proof(q, p)
  return (CONFLICT, e, f1, f2, g)

def assign_canonicals(sum):
  in_left_dict = sum.in_left_dict
  in_b_dict = sum.in_b_dict
  b_index_by_name = {}
  for u in sum.index.by_key_dict.values():
    y = sum.out_right(u)
    if y:
      name = get_canonical(y, None)
      if name:
        b_index_by_name[name] = y
  if debug:
    print("# %s B canonicals" % len(b_index_by_name),
          file=sys.stderr)
  count = sci_count = 0
  losers = 0
  for u in sum.index.by_key_dict.values():
    y = sum.out_right(u)
    if y:
      name = get_canonical(y, None)
      sci_name = get_scientific(y, None)
    else:
      x = sum.out_left(u)
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
# 14. emit the new sum with additional parent column

# Returns a row generator

def emit_spanning_tree(sum, roots, rm_sum):
  if new_style:
    def generate_usages():
      yield union
    generate_rows(generate_usages(),
                  [key_prop,
                   A_key_prop, B_key_prop,
                   parent_id_prop,
                   accepted_id_prop,
                   taxonomic_status_prop,
                   scientific_name_prop,
                   ])
  else:
    yield ["taxonID", "taxonID_A", "taxonID_B",
           "parentNameUsageID", "acceptedNameUsageID",
           "taxonomicStatus",
           "canonicalName", "scientificName", "recordMatch", "change", "remarks"]
    in_left_dict = sum.in_left_dict
    in_b_dict = sum.in_b_dict
    rm_in_left_dict = rm_sum.in_left_dict
    rm_in_b_dict = rm_sum.in_b_dict

  def compose_row(union):
    a_usage = sum.out_left(union)
    b_usage = sum.out_right(union)
    p = get_parent(union, None)
    a = get_accepted(union, None)
    z = m = None
    change = MISSING
    if b_usage:
      z = sum.out_left(mep_get(rm_in_b_dict, b_usage))
      if z:
        rcc5 = cross_relation(z, b_usage, sum)
        change = rcc5_symbol(rcc5)
        m = mep_get(in_a_dict, z, None)
      # if a_usage and b_usage then "renamed"
    else:
      z = sum.out_right(rm_sum.in_left(a_usage))
      if z:
        rcc5 = cross_relation(a_usage, z, sum)
        change = rcc5_symbol(rcc5)
        m = sum.out_left(mep_get(in_b_dict, z, None))
    return [get_key(union),
            get_key(a_usage) if a_usage else MISSING,
            get_key(b_usage) if b_usage else MISSING,
            get_proper_key(p) if p else MISSING,
            get_proper_key(a) if a else MISSING,
            figure_taxonomic_status(union, sum),
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

A_key_prop = prop.getter(prop.get_property("taxonID_A"))  # getter
B_key_prop = prop.getter(prop.get_property("taxonID_B"))

def figure_taxonomic_status(u, sum):
  u_basis = sum.out_right(u) or sum.out_left(u)
  acc = get_accepted(u, None)
  if acc:
    return seniority(u_basis, acc)
  else:
    return "accepted"

def get_proper_key(u):
  if is_top(u): return MISSING
  else: return get_key(u)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    TBD.  stdin = the A checklist
    """)
  parser.add_argument('--target', help="the B checklist")
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
