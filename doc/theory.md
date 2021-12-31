
# Terminology, philosophy, method

## Terminology - rows, taxa, extensions, interpretation

Many of the tools are completely generic over tabular data, but a few
are specific to biodiversity information in the form of "Darwin Core"
files.

When I speak of a Darwin core (DwC) file I take this to mean (for
purposes of these tools) a CSV file where each record (row) has
information connected to a taxon.  More precisely, I take a record to
refer to what I'd call an _extension_ (of a taxon description).  An
extension is simply the set of organisms/specimens/observations as
described or circumscribed or referenced somewhere, perhaps in a
database.  Each record is itself a little taxon description, and may
contain references to outside sources to help constrain an extension.

I tend to use 'taxon' and 'extension' interchangeably but they are
slightly different; a taxon is something subject to taxonomy
(classification) while an extension is specifically a set of organisms,
and may or may not be a taxon.

If the description or reference is vague or incomplete that can be
unfortunate but it is not necessarily a problem.  The names can still
be used for searching outside sources.  Different people might
comprehend the circumscription differently, but they should try to
remain open minded pending further information about what was meant,
and the differences can be ironed out with research if necessary.
([Model theory](https://en.wikipedia.org/wiki/Model_theory) is how I
understand this approach: We do not say exactly which extension we
mean, but we talk about an extension in terms that constrain meaning
adequately for the task at hand.)

Call the extension associated with a record the _interpretation_ of
the record.  The record is to be interpreted in the context of all the
other records in the file in which it occurs, and what they say about
one another, and in the context of whatever we know about the origin
of the file itself.  In the same way we can speak of an interpretation
of a table (set of Darwin Core records) as the harmonious simultaneous
interpretation of all of its records.

Each record contains one or more names, 'identifiers', or name- or
identifier-like strings.  In some files a single name might be used
for multiple distinct extensions (homonyms), but if so each "way" will
have its own record.  Some columns such as `taxonID` may contain
record identifiers unique within the file, but not necessarily
univocal (with a common interpretation) from one file to the next.  A
`taxonID` identifies both a record and its associated taxon/extension
(according to interpretation of the file/database in which it occurs).

Other Darwin Core columns containing record identifiers have
column names that contain 'usage' as a morpheme,
e.g. `parentNameUsageID`.  This shows a confusion about how meaning
and reference work.  Yes, corresponding to each record there are
tokens (names and 'identifiers') whose usage (pattern of use) we want
to track so that we can figure out how to interpret them.  But the
name or identifier most naturally and usefully identifies a taxon, the
biological subjects of our taxonomic efforts.

As with any linguistic entity, meaning can change over time.  Again,
this is both inevitable and not as problematic as many seem to
believe.  The best and only corrective to confusion is citing one's
sources.

## Darwin Core

For algorithmic purposes, Darwin Core is a representation of trees (or
rather forests) with data at each node.  Each node has a locally
unique identifier (`taxonID`) and at most one parent
(`parentNameUsageID`).  Local data includes a shorter taxonomic name
(`canonicalName`) that is often a Linnaean binomial (e.g. "Mustela
erminea"), and a longer 'scientific name' that is the canonical name
followed by a reference (given as author and year, which is usually
enough, for a particular name, to track down the intended publication)
(e.g. "Mustela erminea Linnaeus, 1758").

This is a simplification.  One elaboration that one should know about
is _accepted_ versus unaccepted names or records.  Accepted
names/records are those that can participate in the strict hierarchy,
including parent/child relationships and the understanding that the
children of a single parent have mutually disjoint extensions.

Other names/records, the unaccepted ones, are sometimes called
'synonyms'.  Although the taxa they designate are not different from
any other taxa, overlap between sibling synonyms is permitted.  The
distinction is encoded in Darwin Core in two ways.

   1. The `acceptedNameUsageID` of an accepted record is either empty
      or the `taxonID` of the record itself.  Synonym records are
      those that are not like this.
   1. The `taxonomicStatus` field indicates acceptance via vocabulary
      such as `accepted`, `valid`, `accepted name`, and perhaps others.
      A variety of field values (not these) indicate synonyms, e.g.
      `synonym` or `subjective synonym`.

A well-formed Darwin Core table uses these two methods consistently
with one another.

## Reasoning about extensions

To simplify in unimportant ways, when an author wants to coin a new
name for a taxon, by universal convention they must designate a single
specimen to be the name's "type specimen".  By another universal
convention ('priority rule') the name of a taxon must be that of the
oldest designated type specimen that the group contains.  Different
authors or articles can use a taxon name for various different taxa,
but by the priority all these taxa contain the name's type specimen.
These two conventions are helpful computationally.

In taxonomy it is very helpful to be able to relate taxa to one
another, given their names or records.  Important relations include
the RCC-5 relationships: two records x and y designating taxa u and v
such that either u = v (equivalence), u < v (proper inclusion, also
written ⊂), u > v (⊃), u >< v (conflict: u and v intersect but are not
equal and neither contains the other), or u ! v (u and v are
disjoint).  Once taxa are arranged by inclusion in a tree, deciding
these relationships is very easy just by reading them off from the
tree, because a child/parent edge is interpreted as set inclusion of
the child's extension in the parent's extension.  In particular,
conflict cannot occur inside a single tree.

The merge operation is more complicated.  Here we start with _two_
trees and want to be able to decide statements of relationships
between all records in both trees.  The extensions may each contain
many specimens, many of which are not even collected much less known,
so we cannot enumerate the extensions.  But we can use the type
specimens as a proxy for the 'small' extensions and estimate
relationships based on those.  Each taxon is then interpreted as a set
of type specimens.

In fact, without loss of generality, we needn't consider all the
types, but only those that are in common between the two trees.  These
are called "mutual tipward record matches" or MTRMs because each is
near the tips (leaves) of both trees.  (There are rare situations in
real Darwin Core files where types are associated with taxa that
contain taxa with their own types.  Currently those types are
disregarded in the calculations.

Now we can 'estimate' an extension via its MTRM: it's the set of MTRMs
that are reachable via 'tipward' paths from the taxon.  (If the MTRM
set is empty it treated separately; it is in effect pruned off and
then reattached when needed.  These are called 'trivial' taxa.)  The
real extension is one that contains its designated MTRMs and no others
(among the set of MTRMs for a given pair of trees).

By relying on type specimens in this way, our commitment to names is
kept to a bare minimum.  Names of higher taxa play almost no role; it
is containment of types that determines how higher taxa relate to one another.

When there are two trees, conflict between their taxa (><) can
certainly occur.  In Newick notation this looks like

  - ((a,b)e,c)f vs.
  - (a,(b,c)g)f

Heuristics are applied when two taxa have the same MTRM set, as occurs
in 'monotypic' chains a - b - c with no nontrivial branching to the
side.  In Newick this would look like a merge of, say,

  - (((a)b)d)e  with
  - (((a)c)b)e

Heuristics are required because there may be multiple ways to order
and equate records have the same MTRM set; the ordering is not forced.
It may be facilitated by record matches (as when `b` in the first tree
is determined to be equivalent to `b` in the second tree).  This is
the only situation (I think) other than those involving synonyms where
the RCC-5 theory of two trees can be incomplete.

# Record matching

Hmm... long story, probably uninteresting.  The multiplicity of
name-like information in the record (canonical, scientific, type,
managed identifier) is exploited to get matches in a variety of ways,
and to disambiguate when necessary.  This is critical because MTRMs
are selected based on record matches.

TBD

# Synthesizing a single tree from two

The 'merge' operation.  To select a tree that is a best possible tree
under some measure, conflicts have to be resolved by breaking select
parent/child links, and all incompleteness (see 'monotypic chains'
above) has to be eliminated by making decisions that force a linear
ordering.

TBD
