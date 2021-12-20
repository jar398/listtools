#!/usr/bin/env python3
""" ≥ ≤ ≳ ≲ ≁ ⋧ ⋦ ≵ ≷ ≂ ≌ ≠ ⪝ ⪞ ⊥⊤ ∧∨ ⊔⊓ ∪∩ → ∎ 

Five objectives here.

  #0: Specify how to find 'mutual tipward record matches' or MTRMs,
  which provide a possible starting point for decision and spanning
  tree method.

  #1: RCC-5 decision procedure for the (possibly incomplete) combined
  theory (2 trees linked by MTRMs, defined ff.).

  #2: complete the theory by adding axioms relating within a block.
  The two sides need to be threaded together into a linear sequence.

  #3: breaks some parent/child links in order to eliminate ><.
  Some nodes will want to have 2 parents and can't if we're going to
  make a spanning tree

  #4: Explore optimization options to make comparison adequately fast.

1. Decision procedure

Since we'll be using *some* estimation scheme (at the very least, to
disregard 'peripheral' nodes - those not subtending any MTRMs), there
will be an easily decidable equivalence relation ~ that holds between
nodes having the same estimate.  Write [x] for the equivalence class
that contains x.

First, a 'brute force' decision procedure:

We can consider the RCC5 relations on blocks [x] (equivalence classes)
induced by the RCC5 relations in the original theory.

    [x] = [y]   iff  x ~ y
    [x] < [y]   iff  x ⋦ y, i.e. x < y and not x ~ y
    [x] > [y]   iff  x ⋧ y, i.e. x > y and not x ~ y
    [x] >< [y]  iff  x >< y  and not x ~ y
    [x] ! [y]   iff  x ! y   and not x ~ y

These five relations ~ ⋦ ⋧ >< ! are mutually exclusive.

Because it's an equivalence relation, [x] = [y] (or x ~ y) could mean
that any of the five RCC5 relations hold between x and y - a priori we
have no information.  The finest (most precise) estimation scheme that
has any justification is based on the MTRM sets, each one of which is
the set of MTRMs subtended by some node.  If we do this we can drop
the 'not x ~ y' condition on ! :

    [x] ! [y]   iff x ! y    - could never be justified

That is, if x ~ y, then x and y overlap.

We can easily compute the RCC5 relations on blocks.  Code for
'outlier' and 'overlap' is below.

    outlier(x, y)  iff  some MTRM is subtended by [y] but not [x]
    overlap(x, y)  iff  some MTRM is subtended by [y] and [x]

that is,

    outlier(x, y)  iff  x {⋦ or >< or !} y
    overlap(x, y)  iff  x {⋦ or ⋧ or ~ or ><} y   i.e. x ∩ y

So...

We can decide x ~ y i.e. [x] = [y] with

def similar(x, y):
  return ( not outlier(x, y) and 
           not outlier(y, x) and 
           overlap(x, y) )

or perhaps something simpler than that depending on the estimation scheme.

We can decide {x < y and not x ~ y}, i.e. x ⋦ y, i.e. [x] < [y], with

def lt_not_similar(x, y):
  return ( not outlier(x, y) and 
           outlier(y, x) and 
           overlap(x, y) )

and similarly for the reverse [x] > [y].

We can decide {x >< y and not x ~ y}, i.e. [x] >< [y], with

def conflict_not_similar(x, y):
  return ( outlier(x, y) and 
           outlier(y, x) and 
           overlap(x, y) )

(except for sibling synonyms???)

So now we can decide RCC5 in the blocks theory (with ~ for = and so
on), for which see part 4.  ∎

We might want to substitute, for performance reasons, certain kinds of
estimates for the exact blocks in the process, and revert to MTRM set
comparison when necessary.  It has to be possible to detect cases
where the estimates are not conclusive, and fall back to the 'brute
force' method in that situation.

2. Completion (decision within blocks).

Complete the theory by adding relationships within a block.  This is
prerequisite to creating a spanning tree.

Define the infimum of x and y, written x ∧ y, to be the the largest
element (according to >) that's smaller (according to <) than all
others that are smaller, if there is one.  

  If x = y, then x ∧ y = either x or y (it doesn't matter).
  If x < y, then x ∧ y = x.
  If x > y, then x ∧ y = y.
  If x >< y, then there is no x ∧ y, although we might postulate the
    existence of one (usually called 'bottom') so that we can write 
    x ∧ y = ⊥.
  Similarly, x ! y = ⊥.

The spanning tree problem is to find a parent w = up(z) in the
combined system (both checklists combined).  Algorithms for finding a
spanning tree (whether top down or bottom up) always have two nodes p
and q, one from each checklist, on hand.  These are chosen so as to be
the two best candidates for up(z), one from each checklist, and we aim
to establish up(z) = p ∧ q ≠ ⊥.  This means we require 

  either z < p ≤ q or z < q ≤ p.

(Peripheral nodes are handled separately and trivially.)

By induction, x and y are known to overlap (it's not the case that x !
y).

To synthesize a tree, `up` has to be total, which requires that we
have a theory (as close as possible to the combination of the two
input checklists) that is complete and contains no ><
relationships; this is practically the definition of a tree.  ><
removal we delay until the next section.

For completion, we need to form the members of each ~ block into a
single tree, by adding and perhaps removing axioms according to
heuristics.

First, if ~ contains members with different subtended MTRM sets, take
the quotient of the block over MTRM sets, and treat each sub-block
differently.

Here are some trees to think about, then:

    ((f,((a,b),c)d1,g)e1)h1
    ((f,(a,(b,c))d2,g)e2)h2

The MTRMs are unsubscripted.  e1, e2, h1, and h2 all have the same
MTRM set, namely {f, a, b, c, g}.

For each sub-block, we perform the following procedure.  Start with x,
y = the 'floors' of the block i.e. the minimal (lowest, most tipward)
nodes in the block.

  . Start with x = the least (most tipward) element of the A chain and
    y = the least element of the B chain.  We'll be going up one or
    the other chain (or both chains) (rootward) at each step.  

  . At each step x and y are given.

  . Usually the RCC5 relationship between x and y will not be
    determined by the existing theory-in-progress.  If it is, our
    choice of x, y relationship is constrained by established
    constraint.  If that constraint does not include any of > < = then
    >< is demanded and we have to do conflict resolution (as for e1,
    e2 above).  [There is a conflict resolution method in the next
    section; linking this part of the writeup to that is TBD.]

  . We can assume it's not the case that x >< y, because we went from
    carved up the blocks into sub-blocks if we had to.  Therefore each
    sub-block consists of two linear < (or >) chains, one in each
    checklist.  We will add axioms to totally order the chains,
    compatibly with existing ordering of each side.

  . Put p = x's parent in A, and q = y's parent in B.  We seek to know
    whether x < y or y < x.

  . If x and y are record-matched to one another, posit that x = y, so
    that x ∧ y = x = y.  Blend the peripheral children of x and y.
    (Next: find up(x ∧ y) = p ∧ q.)

  . If p ≁ y and q ≁ x (the candidates are in different sub-blocks),
    then posit that x = y (the "compatible extensions" case, which is
    sort of nice) or that x < y or y < x.  (Next: find up(x ∧ y) = p
    ∧ q.)

  . If p ≁ y but q ~ x, then necessarily x ∧ y = x.  Similarly, if q ≁
    x but p ~ y, then x ∧ y = y.  (This is because we have to exhaust
    all the ~ options before going higher.)
    (Next: find up(x ∧ y) = p ∧ y or x ∧ q).

  . If x is record-matched to y' with y' > y, and y isn't
    record-matched to any x' with x' > x (or vice versa), take the
    matched one to be the inf (so that we will aspirationally maximize
    our chances of a match further up).  To be sure, make sure that
    the match is in same ~ equivalence class as y.  
    (Next: find up(x ∧ y) = p ∧ y or x ∧ q).

  . In all other cases (both unmatched, with or without matches to
    other nodes; parents both in block), choose between <, >, and =
    arbitrarily (y by priority rule, or something about the record
    matches, or any rule really).  (Next: find up(x ∧ y) = p ∧ y, x ∧
    q, or p ∧ q.)

3. Conflict (><) resolution

If x >< y, and we want to find a spanning tree, there is no principled
way to do it.  We need to retract axioms (cut parent/child links in A
and/or B) in order to form a spanning tree.  We choose x ∧ y between
<, >, and = using any ad hoc rule (similar to the ~ case but with more
serious consequences), positing y < x (or x < y or x = y or even x !
y) and retracting some set of links that yield x >< y.  (Next: find
up(x ∧ y) = p ∧ y, x ∧ q, or p ∧ q.)

As above, if there is a record match it will advocate for x = y.

Hmm, this just seems to happen by itself.  We get to a point where a
node z (such that z = x = y) yearns to have to have two parents, we
can't figure out how to choose between them (neither is < the other,
i.e. they are ><; which is symmetric relation), and we have to drop
one of them (the low priority node).

  1. We factor z as z = x ∧ y, and 'split' it so we treat x and y
     separately.  If there is no such equation then we have no problem;
     and none would have triggered by the resolution process.

  2. We make q the parent of y (from priority checklist) and call p the
     'shadow' (or 'alternative') parent of x.  x becomes a junior
     synonym of y (x <= y).

  3. y (in q) takes on all the combined children of x and y ('winner
     takes all'), after which x has no children, which is good because
     x has become a synonym.

4. Efficiency via estimation schemes

Given checklists as disjoint sets A and B (each equipped with parent
function), define an 'estimation scheme' to be

  . a T set of 'estimates',
  . a function η : MTRM -> T (basis case),
  . a binary (associative, commutative) operator ⊓: T x T -> T
  . a function E: A ∪ B -> T defined by:
      . E(m) = η(m)
      . E(x) = E(x1) ⊔ ... ⊔ E(xn), where x1 ... xn are x's
        children
  . an equivalence relation x ~ y defined by: x ~ y iff E(x) = E(y).

Any of these equivalence relation ~ is suitable for use in parts 1-3.

Here are three estimation schemes that I've seen used:

1. Exact, or set-based (Holder and Redelings):
     . T = sets of MTRMs
     . η(m) = {m} if m is an MTRM,
     . e ⊔ f = e ∪ f

2. Range-based (Smasher):
     . Impose a total ordering on the MTRMs
     . T = {[m, n] | m and n MTRMs with m ≤ n}
     . η(m) = [m, m],
     . [m1, n1] ⊔ [m2, n2] = [min(m1, m2), max(n1, n2)]
       with min and max according to the ordering.

3. MRCA-based (smasher conflict analyzer, listtools):
     . T = A ∪ B
     . η(m) = m
     . x ⊔ y = MRCA(E(x), E(y))

Other possibilities:

4. Sets of ranges (i.e. run-length encoded sets).

5. Kalman filters.

(Maybe write: - do we need these? -

    e ≤ f iff e = e ⊔ f
    e ≥ f iff f = e ⊔ f. )

 
[... implementation note ... we need in Python: a new kind of instance,
type 'estimate', so that we can put properties on it (yes? or maybe no
need? there is a function that goes from A-estimates to B-estimates
and back, yes? but maybe that doesn't need to be cached.)

(or maybe we can do without an instance and use a direct
representation.  try it without first.)]

[For typographic fun: ≌ ⪝ ⪞.  x ≌ y : 'x and y have the same estimate'
... but ... I will continue to use the other notation (same notation
for MTRM set comparison and for coarser equivalences)]

"""



# C is for checklist, either sum or source but usually sum
# S is for source, used in that role but conceivably a sum


import sys, csv, argparse
import functools
import util
import property as prop
import checklist

from rcc5 import *

# Reading and writing articulation sets, and doing inference on them

# Proposing alignments that create spanning trees (maximally
# consistent, ><-free theories)

# Read a theory from 4 files (A, B, seed, record matches)

def load_theory(a_iterator, b_iterator,
                seed_iterator=None, rm_sum_iterator=None):
  if rm_sum_iterator == None:
    # Need to copy the iterators!
    (a_iterator, a_iterator_copy) = dup_iterator(a_iterator)
    (b_iterator, b_iterator_copy) = dup_iterator(b_iterator)
    rm_sum_iterator = match_records.match_records(a_iterator_copy, b_iterator_copy)
  A = checklist.load_source(a_iterator, "A")
  B = checklist.load_source(b_iterator, "B")
  # Load record matches
  rm_sum = load_sum(rm_sum_iterator, A, B)
  # Create 'theory object'
  AB = checklist.make_sum(A, B, "", rm_sum=rm_sum)
  init_summand(A)
  init_summand(B)


# S is usually a source checklist, but could be a sum

def init_summand(S):
  (S.get_level, S.set_level) = \
    prop.get_set(prop.get_property("level"), context=A)
  cache_levels(S, S.roots)


def dup_iterator(iter):
  clunk = list(iter)
  return ((x for x in clunk), (x for x in clunk))


# Read/write an articulation set

# Finally need to sort this out...

"""
a way to identify a record in either summand;
hmmmm...

need a nice way to get unique names for records.
this will be hard, but can use .... ?

canonical / scientific / dbid ...

hmm.  maybe do it at the taxon id level?  need tooling to create these
things, but much more reliable and concise that way.

so, taxonID_A and taxonID_B columns as before.

Oh, the merged taxonomy is a very different artifact!!!

So, at least 2 different outputs.

SUGGESTED:
  links inherited from sources  ?  goes without saying
  every = articulation  ('peers')
  consequences of =.  equations between children
  child/parent links in monotypic blocks, for completeness
    ('altParent', 'altAccepted')
  dropped children, to avoid conflict and change >< to > or < or =
    'splice out'  or 'synonymize'

INFERRED:
  'estimates'  (individual or as [a, z] estimates)
  child/parent and synonym/accepted links, many of them

So... I think this is not too bad.  An alignment has
articulations
  - at most one = per node (I think)
  - at most parent or accepted link per node
  - either can override ... how??
  - also, deleted or deprecated conflicting nodes
    (conversion to synonym)
  - mutex combination

OK a list with pairs of TNUs and a single relation seems about right.
But the relationship will always be = or < (or <= for synonyms).
a = b
a < c
x <= c
with reasons:
  tipward match on canonicalName,
  similar extension and match on canonicalName,
  coextensional,
  conflict resolution,

"""

# -----------------------------------------------------------------------------
# RCC5 decision procedure
#   1. for nodes in same tree
#   2. for nodes in different trees

# This seems wrong somehow
# Things would be simpler if only source nodes had levels

def common_summand(C, x, y):    # C is AB or BA
  (x1, y1) = C.split(x)
  (x2, y2) = C.split(y)
  if y1 and y2: return C.B
  if x1 and x2: return C.A
  return None

def le(C, x, y):
  return tree_le(C, x, y) if common_summand(C, x, y) else cross_le(C, x, y)

def lt(C, x, y):
  return lt(C, x, y) and not eq(C, x, y) # ?

def eq(C, x, y):
  return tree_eq(C, x, y) if common_summand(C, x, y) else cross_le(C, x, y)

# relation(C, x, y) and so on

# -----------------------------------------------------------------------------
# Decide about a single source tree

# E.g. level(x) >= level(y)  should imply  x <= y  if not disjoint

def tree_eq(C, x, y):
  return x == y

def tree_le(C, x, y):
  y1 = x     # proceeds down to y
  # if {level(x) >= level(y1) >= level(y)}, then x <= y1 <= y or disjoint
  stop = get_level(y)
  while get_level(y1) > stop:
    y1 = get_superior(y1)    # y1 > previously, level(y1) < previously
  # level(y) = level(y1) <= level(x)
  return y1 == y    # Not > and not disjoint

def tree_lt(x, y): return tree_le(x, y) and x != y

def find_peers(C, x, y):
  while C.get_level(x) < C.get_level(y):
    y = get_superior(y)
  while C.get_level(x) > C.get_level(y):
    x = get_superior(x)
  return (x, y)

# MRCA within the same tree

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
# Decide about relations between two trees

# Returns a pair (overlap, outlier)
#   overlap is in y and x {< = > ><} y (but not !)
#   outlier is in y and x {< >< !} y (but not {> =}) x ⋧ y
#   either can be None

def analyze(AB, x, y):
  (over, out) = (None, None)
  for d in AB.get_children(y):
    (over2, out2) = analyze(x, d)
    over = over or over2; out = out or out2
    if over and out: return (over, out)
  for d in AB.get_synonyms(y):
    (over2, out2) = analyze(x, d)
    over = over or over2; out = out or out2
    if over and out: return (over, out) # (over, None) !!!??? conflict
  if over or out:      # Exclude peripherals
    return (over, out)
  m = AB.record_match(y)
  j = AB.join(m, y)
  return (j, None) if m else (None, j)

def outlier(AB, x, y):
  (over, _) = analyze(x, y)
  return out

# Cache this.
# Satisfies: if y = shadow(x) then x ~<= y, and if x' = shadow(y),
# then if x' <= x then x ~= y, or if x' > x then x > y.
#
# i.e.   x ≲ y ≲ x' ≤ x   →   x ≅ y   (same 'block')
# i.e.   x ≲ y ≲ x' > x   →   x < y

def get_shadow(AB, x, other):
  m = functools.reduce(mrca,
                       (get_shadow(c) for c in get_inferiors(x)),
                       None)
  if m == None:
    return record_match(x, other) or None
  return m

# RCC5 decision procedure, where x and y are in different sources.

def cross_ge(AB, x, y):
  return not outlier(AB, x, y)

def cross_eq(AB, x, y):
  # see if they're peers ??
  return cross_ge(x, y) and cross_ge(y, x)

def gt(AB, x, y):
  return cross_ge(x, y) and not cross_ge(y, x)

# -----------------------------------------------------------------------------
# Cache every node's level (distance to root)
#   simple recursive descent from roots

# Level is contravariant with RCC5: x < y implies level(x) > level(y)

def cache_levels(S, roots):
  def cache(x, n):
    S.set_level(x, n)
    for c in S.get_inferiors(x):
      cache(c, n+1)
  for root in roots:
    cache(root, 1)



if __name__ == '__main__':
  src = load_source(csv.reader(sys.stdin), "A")

  props = (prop.get_property(label)
           for label in ("taxonID", "canonicalName", "parentNameUsageID"))
  writer = csv.writer(sys.stdout)
  for row in emit_rows(src, props): writer.writerow(row)

