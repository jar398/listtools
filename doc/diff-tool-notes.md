
Notes on what a 'diff tool' might be

# Processing steps

Some of these steps correspond to commands that a part of the
`listtools` suite, while others are just code belonging to particular
tools (currently `merge.py` but potentially elsewhere in addition or
instead)

## Cleaning

`start.py` performs a miscellaneous set of functions to ensure
well-formed and canonical inputs:
  1. converts TSV to CSV
  2. checks for presence and uniqueness of primary keys (`taxonID`)
  3. ensures that the `canonicalName` and `scientificName` names contain canonical and
     scientific names (without authority and with authority)
  4. checks for consistency between `taxonomicStatus` and `acceptedNameUsageID`
  5. adds `species` as value of `taxonRank` rank if `taxonRank` is missing and
     canonical name is a binomial
  6. adds a `managed_id` if desired, based on a specified prefix and
     the `taxonID` value, which is useful when matching between
     different versions of a sources such as NCBI or GBIF that has
     managed identifiers.

## Name parsing and stemming

We run `gnparse` on each input to augment an input file with
additional columns based on parsing the `scientificName`, and to set
the `canonicalName` from the `scientificName` when the `canonicalName`
is missing.  One of the additional columns is a stemmed version of the
`canonicalName`, which permits matching of -a names to -us names (for
example).

## Record matches

`match_records.py` takes two prepared CSV files and generates a mapping
between their records, based on the key-like fields of the records
including `canonicalName`, `scientificName`, stemmed name, and
`managed_id`, all of which are optional.  Not all records are matched
necessarily, but when they are a single record is matched to a single
record.  Unambiguous matches can be based on any key-like field, but
if no field allows a unique match, no match is made at all.

## Tipward record matches

Record matches include both those for tips (usually species or
subspecies) and internal nodes (usually higher taxa).  Matching must
begin with just the tips since species can get moved from one genus or
family to another, making higher taxon names unreliable.

So, there is a scan over each of the two inputs checklists to find the
"most tipward" record matches, which form a subset of the complete
list of record matches.  These tipward matches are important, and such
a match will be referred to as a "TRM".  These are not necessarily
tips, since their descendents might be unmatched, but they have the
property that no TRM is a descendant of any other TRM.  Considered as
taxon extensions, they are mutually disjoint.


## The TRMs subtended by a higher node

For each every higher taxon record (node) we can talk about all the
tips that are subtended by that node, or (to say the same thing)
descended from it.  If we take two higher taxa, then a difference in
their tip sets is not very informative, because a tip node might be
present in one checklist and not the other, or might be a tip in one
but not in the other.

However, if we consider instead all the TRMs (which are not
necessarily tips), we know that all such matches likely represent taxa
that are in both checklists, so that a match is good evidence that the
higher taxa are the same concept (or very similar), and a mismatch is
good evidence that the higher taxa are different concepts (because we
can point to a matched node pair where a node is subtended by one of the
higher nodes, and not by the other).


## Matching higher nodes

Matches between higher nodes based on their TRM sets are not
necessarily unique.  In general we want to match a chain of higher
nodes in one checklist, all having the same TRM set, with a similar
chain in the other checklist.  Trivial, unambiguous chains (length 1)
can be matched with high likelihood that the same taxon concept is
intended by both nodes/records, even when there is no record match.
If all the nodes in a chain have record matches in the other checklist
having the same TRM set (i.e. in the corresponding chain), in the same
order, we have a well-evidenced match between the two chains.  In all
other cases the ambiguity has no good automatic resolution and should
be flagged as requiring user input, or resolved arbitrarily with a
warning.


# RCC-5 articulations

Having identified TRMs and higher record matches based on TRM
subtension, we are in a position to determine the RCC-5 relationship
between any node in one checklist and any node in the other.

(For within-checklist questions, we can answer RCC-5 questions by
observing the hierarchy in the usual way: a child may be assumed `<`
its parent, accepted children of the same parent are disjoint, and we
do not know the relationships between any unaccepted name (synonym) and
its siblings.)

One obvious output format would be a set of articulations between
nodes of the two checklists, but which articulations should be
included?  If the checklists have length about n, then there are about
n^2 articulations, which is neither practical nor desirable.  We need
a subset from which, together with the checklists themselves, all
other articulations logically follow; and beyond that there might be
articulations that are especially surprising or interesting.

Here are some articulations that I have identified as interesting.
Assume the checklists are A and B, and that x is in A and y is in B,
and consider also the reverse comparison (B compared to A).

  1. x = y  where x is a tipward record match to y
  2. x = y  where x and y unambiguously subtend the same TRMs (see above)
  3. x < y  where x has no equivalent in B ('addition')
  4. x < y  where y has no equivalent in A ('refinement')
  5. x >< y  where x < parent(y)  ('inconsistency')
  6. x R y  where x and y are record matches and R is not =


# 'Diff' file format

I don't know what the most useful 'diff' format would be.  I would
like to see reports generated in two steps:

  1. Generate a set of articulations, as above, that captures the
     result of alignment process such as the one outlined here; one
     might call this a 'basic report'.
  2. A set of tools that take a 'basic report' as input and
     generates whatever kind of report is needed: summaries,
     diffs (similar to MDD diffs?), merges, diagrams, etc.
