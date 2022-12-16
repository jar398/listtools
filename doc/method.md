
# RCC-5 decision method

The central technical problem is figuring out which RCC-5 relationship
holds between a taxon in one checklist and a taxon in another checklist.

So we might ask for the RCC-5 relationship between taxa x and y, where
x is the taxon in checklist A that some record r1 stands for, and y is
the taxon in checklist B that some record r2 stands for.

We may not know much about these taxa, and we may even have to
speculate on their properties.  But we do know what the checklists say.

Of course we don't have enough information to discern the checklist
authors' intent, so the result is heuristic.  Even if two checklists
are identical in every way, it is still possible for taxon comparison
to be underdetermined, because of some unknown context bearing on the
checklists.  The approach is based on the given hierarchies in the
checklists together with judicious use of names that are common
between the two checklists.

## Overview

1. Checklist preparation including scientific name parsing
1. A 'workspace' holds two checklists
1. Record matches by name (hierarchy insensitive)
1. Heuristic exemplar choice
1. Exemplar set calculation
1. Calculation of 'reflections' (smallest >= taxon in other checklist)
1. 'Theory' = ability to calculate RCC-5 relationship between arbitrary nodes
1. Applications: alignments, merged checklists, diffs

## Checklist preparation

start.py - normalize a checklist file (change to CSV if TSV; check
ids; check consistency of acceptedNameUsageID with taxonomicStatus;
sanity checking on canonicalName and scientificName)

extract_names.py - prepare list of scientific names, to be fed to gnparser

use_gnparse.py - move columns from the gnparser output into the checklist
.csv file (to separate author/year from genus+epithet, and for stemming, in particular)

## Workspace

All subsequent work operates on a data structure holding two
checklists and values calculated from them in combination.

## Record matches

From checklists A and B, form a table of mutually unique matches
between rows of A and rows of B, capturing the most likely
correspondences between rows in the absence of any consideration of
the checklists' hierarchies (parentNameUsageID and
acceptedNameUsageID).  The important thing is the contract, not the
details.  One could easily replace the matcher with another matcher,
in the quite likely event another is better.

Ambiguities are treated the same as non-matches.

A match is sought based on a sequence of columns taken one at a time
in order: 

    scientificName, type, canonicalName, canonicalStem, managed_id

So, if exactly one record x in A has scientificName S, and exactly one
record y in B has scientificName S, then we have a record match.  If
not, then we do the same check but with the values in the `type`
column, then "canonicalName", and so on.

The purpose of the `type` column (probably a poor choice of name) is
to allow species name matches where there is a genus change, meaning
there is no exact match between the scientificName or canonicalName.
The `type` is synthesized from information extracted by gnparse: the
year (of publication), the epithet (could be specific or
infraspecific) as stemmed by gnparse, and the last name of the first
author.

In theory, two records could match uniquely based on their `type` yet be
completely unrelated, but in practice this hasn't yet happened.  And
of course this would be a risk for any column.

`managed_id` allows records from checklists using a common managed
taxonID space (such as NCBI Taxonomy Ids or GBIF taxon ids) to match
by taxonID when all else fails.

## Exemplar choice

An "exemplar" is intended to be (or intended to 'represent', if you
are uncomfortable with realism) a single specimen, such as a type specimen, whose
classification in both of the two source checklists is known.  Our
RCC-5 comparison between taxa in different checklists will rely, in
part, on set comparisons of exemplar sets.

Note that just because two taxa (from the same or different
checklists) have the same exemplar, does not mean they are the same taxon.
Sometimes they might be, sometimes not, depending on their circumscriptions.

For each exemplar there are two taxa, one in each checklist, each of which is
the smallest taxon in the checklist known to contain that exemplar.
Let's call these the exemplar's "minimal" taxa.

Exemplars are chosen as follows.  If (x, y) is a name-matched pair
(see above), and either 
no descendant of x belongs to a name-matched pair OR
no descendant of y belongs to a name-matched pair, then we hypothesize an
exemplar, with x and y as its minimal taxa.

## Exemplar set calculation

For each taxon in each checklist, we cache (in the workspace) the set
of exemplars that occur in that taxon (i.e. at or below it in the
checklist's hierarchy).  This is a simple postorder union computation.
These sets are called "blocks" because it is sometimes the case that
multiple taxa (even within a single checklist) have the same exemplar
set, so they must be distinguished in some other way.

The word "block" is short for "block of a partition" where the
equivalence relation on the partition is exemplar set equality.

A taxon whose exemplar set is empty is called "peripheral".

## Cross-MRCA calculation

For each taxon in each checklist, we cache the smallest (fewest
members, most 'tipward' in the hierarchy) taxon that whose exemplar
set is equal to or larger than that of the given taxon.  This is done
as a simple postorder MRCA computation grounding out with the
exemplars found on their minimal taxa.

## "Reflection" calculation

This is similar to cross-MRCA but attempts to respect taxon inclusion
and not just exemplar set inclusion.  That is, it distinguishes
between taxa within a block, unless a cross-MRCA which is common among
all members of a block.

The calculation of the reflection of a record x has four cases,
checked in sequence:

1. If x is peripheral (has an empty exemplar set), x's
reflection is the reflection of the smallest ancestor taxon of x that has a
nonempty exemplar set.
1. Let y = x's cross-MRCA.  If y has a larger exemplar set than that
of x, we take y to be x's reflection.
1. If x has a name match z in its block, we take z to be x's
reflection.
1. Let y' = the highest (largest) ancestor of y in the x's block.
If: (a) neither x nor y' have matches by name inside the block, and
(b) x's parent is also outside the block, then we take the 
to be the reflection to be y'.  This enables "topological" matches between
some records that have unrelated names.
1. Otherwise, we have a difficult, obscure, and uninteresting situation that
I won't to explain (exercise to reader: draw a picture), and the method
gives the answer as unknown.  The 'true' reflection is y or some ancestor in the same 
block; we just haven't bothered to figure out which one (and that is 
not always even possible).
It would be possible to improve on this, at a cost in
complexity, and with little benefit because the situation doesn't come
up much in practice.

## RCC-5 decision procedure

While RCC-5 is in principle a matrix keyed by taxa from the two
checklists, it makes no sense to ramify it as an actual stored matrix.
Checklists can be quite large, so this is not practical.  Fortunately
it isn't necessary either.  Instead we have a 'virtual matrix'
consisting of a two-argument function that can compute which RCC-5
relationship holds between a given record in one checklist and a given
record in the other.

The cases for determining the relationship are roughly as follows:
1. Sibling synonyms can relate to one another in any manner at all.
1. Each synonym is <= to its accepted taxon; we don't know whether < or =.
1. Peripheral nodes require degenerate special handling.
1. If the nodes are in different blocks, we use the result of the
exemplar set comparison into an RCC-5 relationship.
1. For comparisons within a block, the reflections are consulted.
x is compared to y's reflection, and y to x's reflection, and the relation is the conjunction of these two resulting within-checklist relationships.

## Applications

1. Euler/X alignment
1. Merge - create a new checklist combining the taxa of the two inputs
1. Diff - show what changes would need to be made to the first in order to obtain the second

# Previous work

I learned a lot of these ideas by working on Smasher, but the current approach to alignment is much more careful and principled.  Smasher was too aggressive in its collapsing of synonyms, and it was limited in that it only cared about merging taxonomies.
Rees JA, Cranston K. Automated assembly of a reference taxonomy for phylogenetic data synthesis. Biodiversity Data Joutnal 2017;(5):e12581. Published 2017 May 22. https://doi.org/10.3897/BDJ.5.e12581

# Acknowledgment

Mark Holder taught me that is was practical to use sets for these calculation, in spite of the large sizes of some of these sets.  His implementation uses bit sets, while mine uses python3 sets (whose details I do not know).
https://doi.org/10.7717/peerj.3058  (I think that would be the one?)
