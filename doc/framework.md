# List Tools what I'm talking about

**OUT OF DATE**

## Overview

The 'tools' or 'scripts' or 'programs' in this repository manipulate
tabular data in the form of CSV files.  Many of them act as filters
and can be invoked from the Unix shell.  They can be sequenced into
pipelines using '|' if desired.  CSV files uniformly have a single
header row giving column names.

I have found it convenient to run the tools by using `make` applied to
a suitable 'makefile'.  This way, intermediate files can be cached
in the local file system, avoiding expensive regeneration steps when
inputs have not changed.

For a complete example, see [example.md](example.md).

### Polarity note

As is customary in computer science, trees grow downwards, like a
plant's root system, not upwards towards the sun.  So:

 * The root is at the top of the tree; we go up
   the tree to get to the root.
 * As sets of individual specimens/observations, taxa are therefore
   'bigger' or 'higher' toward the top/root, and 'smaller' or 'lower'
   toward the bottom/tips.
 * As an upper semilattice, the root is the 'top' node, and
   is greater than all others; the tips are smallest.
 * A < B means that A is descended from B, contained in B, or 'below'
   B.

As is customary in biology, we say 'tip' rather than computer
science's 'leaf', in case someone is tempted to take 'leaf' too
literally.

## Semantics

### Individuals

I'll use the neutral word 'individual' for the entities that we are
classifying; it's meant to subsume various terms in taxonomy and
biodiversity informatics such as 'occurrence', 'organism', 'specimen',
and 'observation'.

### Taxon concepts

Each record of each checklist file (see below) has an associated "taxon concept"
described by the fields of the record with the checklist as context.
(Full details of the taxon concept would be easily available if the
record linked to appropriate taxonomic literature, but this is often
not the case.  We just make do with what we have, which is
parent links, disjointness of sibling accepted taxa, and the ability
to search the literature by the names provided.)

Say that an individual "falls under" a taxon concept to mean that the
individual is classified under, or belongs to, that taxon concept, or
that the taxon concept applies to the individual.

I'll write '_Hyla_ sec. C' (for example) to denote the taxon concept described in
source C under the name '_Hyla_'.  (C might be a checklist.)

The individuals falling under a taxon concept are called its 'extension'.

According to the way the word 'taxon' is typically used, a taxon is
tied to a name, and its extension changes following speech acts such
as 'redescription' or 'lumping' or 'splitting'.  That is, individuals
enter and leave a taxon as a result of human activity.  By contrast, a
taxon concept is rigid: whether an individual falls under a taxon
concept does not change over time and is insensitive to what happens
in the taxonomic literature.

### RCC-5

RCC-5 (region connection calculus) is a simple system for classifying
relationships between 'regions'.  What constitutes a region depends on
how RCC-5 is being used - we could be talking sets, or geographic
regions - but in our case we're talking about taxon concepts.  
Exactly one of the five RCC-5 relationships holds between any
two regions.  These relationships are
 * A = B: A and B are the same, at least inasmuch as the same individuals fall under both
 * A < B: A is contained in B but isn't the same as B; 
   the individuals that fall under A all fall under B but not vice versa
 * A > B: same as B < A
 * A ! B: A and B are disjoint; no individual falls under both
 * A >< B: A and B overlap but neither is contained in the other.
   There are individuals falling under both A and B, and under just A,
   and under just B.
   By convention we'll refer to this situation as 'overlap'.

In cases where the precise RCC-5 relationship isn't known, we can also write
 * A <= B: A < B or A = B
 * A >= B: A > B or A = B
 * A not! B: A and B intersect; equivalently, are not disjoint;
   equivalently, it is not the case that A ! B;
   equivalently, at least one individual falls under both
 * A ? B: RCC-5 relationship is unknown

Applied to taxon concepts, we write 
 * A = B when extension(A) = extension(B)
 * A < B when extension(A) ⊂ extension(B) (A's extension is a proper subset of B's extension)
 * A > B when extension(B) ⊂ extension(A)
 * A ! B when extension(A) ∩ extension(B) = ∅
 * A >< B otherwise


### Exemplars

Call an individual an _exemplar_ for a comparison of checklists A and
B if taxon concepts in both A and B are known with the property that
the individual falls under both of the taxon concepts.  The individual
is proof that the extensions of the taxon concepts intersect
(i.e. have one or more individuals in common).

'Exemplar' is a semantic ideal rather than an operational notion.  To
make the exemplar idea practical when we don't have direct information
about individuals and the taxon concepts under which they fall, we might
use type specimens as the individuals in the analysis, and take the
type specimens to be exemplars when their classification is known in both checklists,
i.e. when their associated names are the same or otherwise "match".
This is explained below under the `exemplar.py` command.

The pragmatic representation of an exemplar, when exemplars are based
on type specimens, is therefore a set of "matched" records sharing the
same type specimen, and containing at least one record from each
checklist.

### Exemplar sets

The central problem in comparing checklists is comparing the taxon
concepts of one with the taxon concepts of the other.

We take sets of exemplars to be computable approximations to
taxon concepts, in the sense that if

&nbsp;&nbsp;&nbsp;&nbsp;S = E(C) = {e: e is an exemplar in C}

for taxon concept C, and similarly

&nbsp;&nbsp;&nbsp;&nbsp;T = E(D) = {e: e is an exemplar in D}

for taxon concept D,
then the set relationship between S and T suggests
an RCC-5 relationship between C and D and vice versa.
If S and T are different, e.g. S ⊂ T, then "suggests" means it's a
good heuristic bet that C < D (and similarly for the other relationships >, ><,
and ! and their set equivalents).  S = T is not informative on its
own, but there may be ways to
deduce the relationship between C and D based on other information in
the checklists (such as parent links and process of elimination).
Having "enough" exemplars reduces the chance that S = T when C ≠ D,
in which case C = D is a good heuristic bet (parsimonious
assumption) when S = T.


