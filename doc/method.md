
# RCC-5 decision method

The central technical problem is figuring out which RCC-5 relationship
holds between a node in one checklist and a node in another checklist.  For example

Of course we don't have enough information to discern the checklist
authors' intent, so the result is heuristic.  Even if two checklists
are identical in every way, it is still possible for taxon comparison
to be underdetermined.  The approach is based on the given hierarchies
in the checklists together with judicious use of names that are common
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

use_gnparse.py - move columns for gnparser output into the checklist
.csv file (stemming)

## Workspace

A data structure holding two checklists and values calculated from
them in combination.

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

    ["scientificName", "type",
     "canonicalName", "canonicalStem",
     "managed_id"]

So, if exactly one record x in A has scientificName S, and exactly one
record y in B has scientificName S, then we have a record match.  If
not, then we do the same check but with the values in the "type"
column, then "canonicalName", and so on.

The purpose of the "type" column is to allow species name matches
where there is a genus change, meaning there is no exact match between
the scientificName or canonicalName.
The "type" is synthesized from information extracted by gnparse: the
year (of publication), the epithet (could be specific or
infraspecific) as stemmed by gnparse, and the last name of the first
author.

In theory, two records could match uniquely based on their type yet be
completely unrelated, but in practice this hasn't yet happened.  And
of course this would be a risk for any column.

## Exemplar choice

An "exemplar" is intended to be (or intended to 'represent', if you
prefer) a single specimen, such as a type specimen, whose
classification in both of the two source checklists is known.  Our
RCC-5 comparison between taxa in different checklists will rely, in
part, on set comparisons of exemplar sets.

Note that just because two taxa (from the same or different
checklists) have the same type, does not mean they are the same taxon.
Sometimes they might be, sometimes not.

For each exemplar we identify a taxon, one in each checklist, that is
the smallest taxon in the checklist known to contain that exemplar.
Let's call these the exemplar's "minimal" taxa.

Exemplars are chosen as follows.  If (x, y) is a name-matched pair
(see above), and either there are no name-matched pairs tipward of x
in the hierarchy OR there are no name-matched pairs tipward of y in
the hierarchy, then we create an exemplar, with x and y as their
minimal taxa.

## Exemplar set calculation

For each taxon in each checklist, we cache the set of exemplars that
occur in that taxon (i.e. at or below it in the checklist's
hierarchy).  This is a simple postorder union computation.  These sets
are called "blocks" because it is sometimes the case that multiple
taxa (even within a single checklist) have the same exemplar set, so
they must be distinguished in some other way.

The word "block" is short for "block of a partition" where the
equivalence relation on the parition is exemplar set equality.

## Cross-MRCA calculation

For each taxon in each checklist, we cache the smallest (fewest
members) taxon that whose exemplar set is larger than or equal to that
of the given taxon.  This is done as a simple postorder MRCA
computation grounding out with the exemplars found on their minimal
taxa.

## "Reflection" calculation

This is similar to cross-MRCA but attempts to respect taxon inclusion
and not just exemplar set inclusion.

The calculation of the reflection of a record x has four cases,
checked in sequence:

1. For "peripheral" taxa (x has an empty exemplar set), x's
reflection is the reflection of the smallest ancestor taxon with a
nonempty exemplar set.
1. Let y = x's cross-MRCA.  If y has a larger exemplar set than that
of x, we take y to be the reflection.  (Otherwise the exemplar sets
are the same, i.e. x and y are in the same "block".)
1. If x has a name match z in the same block, we take z to be the 
reflection.
1. Let y' = the highest (largest) ancestor of y in the same block.
If: (a) neither x nor y' have matches by name inside the block, and
(b) both x's parent is also outside the block, we take the cross-MRCA
to be the reflection.  This enables some "topological" matches between
some records with unrelated names.
1. Otherwise, we have a difficult and obscure situation that
I won't bother to explain (exercise to reader: draw a picture), and we
give the answer as unknown.  It is just y or some ancestor in the same 
block.

## RCC-5 decision procedure

While RCC-5 is in principle a matrix keyed by taxa from the two
checklists, it makes no sense to ramify it as an actual stored matrix.
Checklists can be quite large, so this is not practical.  Fortunately
it isn't necessary either.  Instead we have a 'virtual matrix'
consisting of a two-argument function that can compute what the RCC-5
relationship is between a record in one checklist and a record in the
other.

The cases are roughly as follows:
1. Sibling synonyms can relate to one another in any manner at all.
1. Each synonym is <= to its accepted taxon; we don't know which.
1. Peripheral nodes require degenerate special handling.
1. If the nodes are in different blocks, we turn the set comparison into an RCC-5 relationship.
1. For comparisons within a block, the reflections are consulted.
x is compared to y's reflection, and y to x's reflection, and the relation is the intersection of the two resulting within-checklist relationships.

# Previous work

Rees JA, Cranston K. Automated assembly of a reference taxonomy for phylogenetic data synthesis. Biodivers Data J. 2017;(5):e12581. Published 2017 May 22. https://doi.org/10.3897/BDJ.5.e12581

# Acknowledgment

Mark Holder taught me that is was both convenient and practical to use sets for these calculation.
[need to cite him - might be https://doi.org/10.7717/peerj.3058 ]
