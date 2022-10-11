
Notes on what a 'diff tool' might be

# Processing steps

Some of these steps correspond to commands that a part of the
`listtools` suite, while others are just code belonging to particular
tools (currently `merge.py` but potentially elsewhere in addition or
instead)

## Input

Two checklists in Darwin Core format.  Rows are for either accepted
names or synonyms.  Each synonym links to its corresponding accepted
row via `acceptedNameUsageID` per Darwin Core.  The `canonicalName`
extension is also processed; this column is used in GBIF and gives the
`scientificName` minus the authority and year information (so, usually
a binomial or uninomial or trinomial, but also _var._, temporary taxon
names, and so on).

Within each checklist, it is assumed that every record denotes a
distinct extension, i.e. for any pair of records there is some
specimen that falls under one record but not the other.  This may not
be always be true (e.g. when considering objective synonyms in the
absence of lumping or splitting), but I don't think this will matter
much in practice; it is a simple assumption to explain and reason about.

Between checklists, extensions may or may not be shared, and in
particular a single name might be found through analysis of hierarchy
or synonyms to have different extensions in the two inputs.  The tool
is dealing centrally with taxon concepts, only using names
heuristically to try to match concepts between the two checklists.

Assume throughout the following that 'A' and 'B' are the two
checklists that are the inputs to the tool.


## Cleaning

On each checklist input, `start.py` performs a miscellaneous set of functions to ensure
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

Things to note:
  - a synonym can match either a synonym or an accepted
  - a subspecies might match a species or synonym, if no better option is available

For example, a record in checklist A with `canonicalName` _Saguinus
labiatus thomasi_ (a subspecies) could record-match to a record with
`canonicalName` _Saguinus thomasi_ in B, but only if _Saguinus
thomasi_ does not occur in A (which it probably does not).

## Tipward record matches

Record matches include both those for tips (usually species,
subspecies, or synonyms) and internal nodes (usually higher taxa).  Matching must
begin with just the tips since species can get moved from one genus or
family to another, making higher taxon names unreliable.

For purposes of this analysis, synonyms are treated the same as other
descendants.  For example, if a species S has a synonym T that is a
record match to T', then T to T' will be a tipward record match, while
S isn't.

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

Note: If a species has synonyms or subspecies that are TRMs, then the
species will not itself be a TRM.  This may seem odd - it implies that
we might be matching species records to one another entirely on the
basis of their synonyms, with the species names not mattering very
much.


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

For within-checklist questions, we can answer RCC-5 questions by
observing the hierarchy in the usual way: 
 * an accepted child is `<` its parent
 * accepted children of the same parent are disjoint (`!`)
 * a synonym (unaccepted name) is either `<` or `=` its parent, i.e. `<=`
 * we do not know the relationships between a synonym and any or its siblings;
   it could be any of the five RCC-5 relationships.

These assumed within-checklist relations, together with the
equivalences between tipward record matches, are the basis on which
articulations between the two checklists are inferred.

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
  3. x < y  where x and its descendants have no equivalent in B ('addition' of x to y)
  4. x < y  where y has no equivalent in A (y is a 'refinement' of x's parent)
  5. x >< y  where x < parent(y)  ('inconsistency', maximal)
  6. x R y  where x and y are record matches, and R is one of the RCC-5
     relations other than =

[This documentation is not current]


# 'Diff' file format

I don't know what the most useful 'diff' format would be.  I would
like to see reports generated in two steps:

  1. Generate some set of articulations, as above, that captures the
     result of alignment process such as the one outlined here; one
     might call this a 'basic report'.
  2. A set of tools, each of which takes a 'basic report' as input and
     generates whatever kind of report is needed: summaries,
     diffs (similar to MDD diffs?), Euler/X syntax, merges, diagrams, etc.


# How lumping and splitting are detected

If we suppose that checklist B is a later 'version' of a curated
checklist, and checklist A is an earlier 'version', we can talk about
differences between A and B as resulting from 'changes' or 'events' to
the curated checklist.  The most common and important of these events
are 'lumping' and 'splitting' events.  A simple lumping or a splitting
event would be reflected in the analysis as follows:

* Suppose species S in checklist A has synonym or subspecies T
* Suppose that T is a tipward record match to T'
* Suppose that S' and T' are siblings
* Suppose that S and S' are record matches (not tipward)
* We can infer that 
  the name (of S in A, S' in B) has been split into S' and T' in B.  
  That is, the name refers to different 'concepts' in the two checklists.

Similarly, if the direction of change is reversed (B changed to become
A), we would say the S' and T' are lumped to form S.

An obvious limitation is that even if S' + T' is a split of S, if
there is no record T in A matching T', then we do not know whether the
new record T' is split off from S, or split from one of S's sibling
species.  And of course we also can't distinguish the case where T' is
split from S or a sibling from the case where T' is "wholly new"
(newly discovered).

The presence of T is a matter of luck and would typically result from
the split undoing a previous lump, i.e. the restoration of a previous
classification in which T had been accepted, but later lumped in with
S.  Or, in a very similar manner, it could also result from a
checklist being aggregated from multiple classications with some
classifications taking T' as a synonym and others taking it as
accepted, with the synonym classification given priority, initially,
and the accepted classification gaining priority later on.

In the aggregators (NCBI, GBIF, and so on), this problem generally applies to splitting events
but not to lumping events, because when records are lumped, the
non-type source records are demoted from accepted to synonym, and the
synonym records are retained in the newer checklist version.  If our
checklists do not have curated synonyms or subspecies, then lumps and
splits will only be detected at ranks higher than species (genus, family, etc.).
