# List tools user guide

**OUT OF DATE**

## Overview

The 'tools' or 'scripts' or 'programs' in this repository manipulate
tabular data in the form of CSV files.  Many of them act as filters
and can be invoked from the Unix shell.  They can be sequenced into
pipelines using '|' if desired.  CSV files uniformly have a single
header row giving column names.

I have found it convenient to run the tools by using `make` applied to
a suitable 'makefile'.  This way, intermediate objects can be cached
in the local file system, avoiding expensive regeneration steps when
inputs have not changed.

For a complete example, see [example.md](example.md).

## Installation
<a name="installation"></a>

The tools require no particular installation step.  They can be
invoked from the shell out of a clone of this `git` repository.

Python 3 is assumed.

Install the python3 `regex` module using `pip3`.

Also required is `gnparser`, which has download instructions
[here](https://github.com/gnames/gnparser#installation).


## Semantics

### Individuals

I'll use the neutral word 'individual' for the entities that we are
classifying; it's meant to subsume various terms in taxonomy and
biodiversity informatics such as 'organism', 'specimen', and
'observation'.

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


## File formats

### Checklist file format

A checklist is given as a tabular text file.  If the extension is
`.csv` it is assumed to be CSV (comma-separated); otherwise it's
assumed to be TSV (tab-separated).  Each row other than the header is a 'record'.

A checklist should use [Darwin Core
Taxon](https://dwc.tdwg.org/terms/#taxon) column headings for
information of special significance to these commands.  Other columns
can be included and will be ignored or passed through as appropriate.

A prefix `dwc:` on column headings is optional.

`canonicalName` is an additional column that is not from Darwin Core
but is important.

`managed_id` is a column used internally to the tools but is not
relevant to most users.  (It is used for distinguishing the case of
aligning checklists with a 'managed' identifier space (e.g. versions
of GBIF or NCBI Taxonomy) from checklists that could have accidental
identifier collisions.)


Here are Darwin Core headings used by one or more of the tools.

 - `taxonID`: the record's primary key, uniquely specifying a record within a checklist
 - `scientificName`: the full taxonomic name, with or without authorship information
 - `canonicalName`: the taxonomic name without authorship
 - `scientificNameAuthorship`: the authorship (e.g. `Smith, 1825`).
   Optional but the matcher is more accurate if authorship is present
   either here or in `scientificName`
 - `namePublishedInYear`: year of publication (e.g. `1825`) (not
     currently used; optional)
 - `taxonRank`: `species`, `subspecies`, `genus`, and so on.
   Important but not required
 - `taxonomicStatus`: if `accepted`, `valid`, or `doubtful`, the name is
   to be considered not a synonym in this checklist.  Otherwise
   (e.g. `synonym`) it is treated as a synonym.
   Case matters.
 - `parentNameUsageID` - `taxonID` of the parent record, or empty if a root
 - `acceptedNameUsageID` - `taxonID` of non-synonym record of which
     this is a synonym, empty if
     not a synonym

One of `scientificName` or `canonicalName` must be given.

Putting the `taxonID` = the `acceptedNameUsageID` is another
way to indicate something is not a synonym.


### Exemplar file format

Some commands read or write files that specify exemplars.  An exemplar
file can be either computed by the `exemplar.py` command or provided
independently if there is some other way.  This section
should be of interest if you don't want or don't need to use
`exemplar.py`, for example if you have your own name matcher.

An exemplar file has one output row for each statement
that an exemplar falls under some taxon concept.

 - `checklist`: 0 for the A checklist, 1 for B
 - `taxonID`: the checklist record for a taxon concept
 - `exemplar id`: identifies an exemplar, locally to this file
   (not globally)
 - `canonicalName`: given only for purposes of human inspection or debugging

By construction, each exemplar id will have at least one exemplar file
row giving an A record for a taxon concept the exemplar falls under, and one
giving a B record for a taxon concept it falls under.

If preparing your own exemplars rather than using `exemplar.py` to do so:

First, note that if an exemplar belongs to a taxon concept, it
also belongs to that taxon concept's ancestors.  It is not necessary
to list all the ancestors as taxon concepts containing the exemplar.
However, it is not harmful to include a few extras (e.g. the species when
the type subspecies of the species is listed).

Second, it is desirable to list all of the most specific taxon
concepts containing an exemplar and not just some of their ancestors.
E.g. if an exemplar belongs to a species then it may also belong to
the species's type subspecies, and the subspecies should be listed as
a taxon concept containing the exemplar.

In a checklist of genera with no species, the exemplars would be
chosen one per matched genus name pair.  In general, one should match
matchable names whenever no descendant matches, regardless of rank.


## The tools

To learn about a given tool, please supplement what follows with
whatever you learn from using the `--help` option, e.g.

    src/clean.py --help

Some tools operate on arbitrary CSV files, while some assume they're
working with Darwin Core files.

### clean

A Darwin Core data source should be run through `src/clean.py` first
for some modest validation and DwC-specific cleaning.

    src/clean.py --input A.dump/Taxon.tsv >A-clean.csv

`--input` followed by a file name specifies the location of the
uncleaned Darwin Core file.  This is assumed to be in CSV
(comma-separated values) format if the file name ends in `.csv` and
TSV format otherwise.  The output is always CSV.

`--pk columname` specifies the column containing primary keys; it
defaults to `taxonID`.  If there is no such column, one is added, and
if a primary key is missing, a new one is generated by hashing the
contents of the row.

Values in the primary key column that occur there more than once
are detected and flagged.

Data cleaning is performed for the following columns:
 - `canonicalName`, `scientificName` - if a scientific name (one that has an
   authority) is found in the `canonicalName` column, and the
   `scientificName` column is empty, then the scientific name is moved
   to the `scientificName` column.  Similarly, if a non-scientific name
   is found in the `scientificName` column, it's moved to the `canonicalName` column.
 - `acceptedNameUsageID` - if it just replicates the `taxonID`, clear it
 - `taxonomicStatus` - flag if a record with taxonomic status
   `synonym` does not have an accepted taxon id, or if one with status
   `accepted` does
 - `source` - cleanup specific to EOL DH 0.9 and smasher - remove
   source record ids other than the first
 - `Landmark` - recode values, change to `landmark_status` - EOL specific cleanup

`--managed prefix:column` is for designating use of managed identifier
spaces.  If one column contains, say, NCBI taxids, or GBIF taxon
identifiers, or anything similar, that are stable across versions of
the source (as opposed to being idiosyncratic to one version), then
the column should be copied to the `managed_id` column.  This
operation is what this feature is for.  For examples, `--managed
gbif:taxonID` means that the taxonID column contains managed GBIF
taxon 'identifiers' and the `managed_id` column will contain 'managed
identifiers' (an idea I made up).  E.g. if a row's taxonID contains
`359` then the string `gbid:359` will be placed in the `managed_id`
column.  This will then be used for matching operations (well... not
currently... but it has done so in the past).

### find_taxa

Locates the Darwin Core Taxon file within a .zip file.  Let `A.dump/`
be a directory containing the files resulting from unzipping the .zip
file.  Then:

    src/find_taxa.py A.dump

writes the name of the taxon file within the .zip file, e.g. 

    A.dump/Taxon.tsv

It's usually clear by inspection which file is the taxon file, but
spelling details vary.  This command is intended for use in scripts.

The taxon file is suitable as input to `clean.py`:

    src/clean.py `src/find_taxa.py A.dump`

TBD: This currently operates by examining file names heuristically.
But it really ought to look in the meta.xml file, which provides the
taxon file name explicitly.


### extract_names

This extracts `scientificName`s in a form suitable for consumption by `gnparser`.
If there is no `scientificName` then the `canonicalName` is extracted.

    src/extract_names.py < A-clean.csv > A-names.txt
    gnparser -s < A-names.txt > A-gnparsed.csv

### use_gnparse

This consumes the output of `gnparser` and combines it with the table
that was the input to `extract_names`, enriching the table with the addition of 
new columns copied from the `gnparser` output.

    src/use_gnparse.py < A-gnparsed.csv > A.csv

### exemplar
<a name="exemplar"></a>

Writes a file giving memberships of exemplars in taxon concepts to
standard output.

 * `--A` filename  - the A checklist.
 * `--B` filename  - the B checklist.

For example,

    src/exemplar.py --A A.csv --B B.csv > AB-exemplars.csv

To do this, noting that each record in each checklist has a taxonomic
name, we posit that each record has a type specimen (the type specimen
for its name).  The type specimen is then an individual that, if found
in taxon concepts in both checklists, can be taken to be an exemplar.

To determine that a type specimen in A is an exemplar we must find a record
in B that has the same type specimen.  That is, the
records must be 'matched'.  (Note that it is records, not their taxon
concepts, that are matched.)  The matcher understands changes in
genus, changes in gender, and other vagaries of the nomenclatural
system leading to name changes, as well as the possibility of collisions (homonyms).
But most of the time name matches are exact.

The meaning of an output row is that the individual (exemplar)
identified by `exemplar id` falls under the taxon concept / record identified
by `taxonID` in the indicated checklist.  (Of course the same exemplar
can also fall under other taxon concepts, in particular taxon concepts
in the other checklist.)

Sample output: [col-19-23-exemplars.csv](col-19-23-exemplars.csv)

When the exemplars are type specimens, they play a role similar
to that played by protonyms in the Pyle/Remsen formulation.


### align
<a name="align"></a>

    src/align.py --A A.csv --B B.csv --exemplars AB-exemplars.csv

This writes an analysis report to standard output.  Concepts are
imputed for all of the A and B records, and RCC-5 relationships
between those concepts are estimated.

If you're doing regression analysis, think of A as the 'baseline' 
checklist and B as 'proposed successor'.  (This is not the only use case.)

If `--exemplars` is not given, the exemplars are computed just as the
`exemplar.py` command would.

The checklist inputs should be derived through a pipeline that begins with `clean.py`.

When there is a B concept of the same taxonomic
name (allowing for non-semantic changes like genus moves and
spelling corrections) the B concept is called the A concept's "buddy".

The output (to standard output) of `align.py` has these columns (subject to change):
 - `A taxon id` - The taxon id of an A row
 - `A taxon name` - The canonicalName of that A record (for human
   consumption).  A suffixed `*` indicates a synonym.
 - `operation` - short description of what "happens" to the name
   and/or concept as one "changes" the A checklist to its successor, the B checklist.
   (This is not to say that a temporal order between the
   checklists is required in reality.)
   Work in progress.
 - `B species that intersect` - 
   If the A record indicates rank 'species', this is a semicolon-separated
   list of relationship/id/name for
   B taxon concepts with rank 'species' that are inferred to intersect the A 
   taxon concept.
   The relationship is the RCC-5
   relationship of the A concept to the B concept, the id
   is the taxon id of the B concept's record, and the name is
   canonical name from the B record.
   An initial `.` prevents Excel from treating values as formulas.
   A value of `-` means there may be intersecting concepts but the list was not computed 
   because the A row was not for a species.
 - `LUB in B` - 
   The relationship/id/name of the A concept to its least upper bound (LUB)
   in the B checklist.  The least upper bound is the smallest B concept
   that contains the all of the A concept.
   An initial `.` prevents Excel from treating values as formulas.
   If the A
   and B names are accepted, the A concept is either the same (RCC-5 =)
   as the B concept or smaller (RCC-5 <) than it.  For synonyms it may
   be hard to tell what the precise relationship is so it <= or ? will show.
 - `A and B` - list of ids of exemplars that occur in both the A
   concept and the name's B concept (the "buddy" concept)
 - `A not B` - list of ids of exemplars that occur in the A
   concept bot not the name's B concept
 - `B not A` - list of ids of exemplars that occur in the 
   concept bot not the name's A concept


A name written with an asterisk (e.g. `Rana pipiens*`) indicates a synonym.

The canonical names in the output are there for human readability.
For a more compact ('normalized') report they might be omitted, and obtained as
needed from the checklists using the taxon id.

Report rows for synonyms that match synonyms or that have no match
nothing are suppressed.  (this may need to be revised.)

How to read the report:

If you're mainly concerned with the impact on a data set of advancing
from one version of a checklist to the next (from A to B), then focus
on lines with semicolons, i.e. intersections with more than
one B concept.  These rows indicate splits, meaning that data using the
taxon name in A would have to be re-curated to use the correct
intersecting species (concept) in B, if one wanted to make the data consistent
with checklist B.

Probably of most interest for understanding the impact of changes in
taxonomy going from A to B are the rows with multiple concepts given in
the intersecting concepts column.  These are situations where an A
species concept corresponds to multiple B concepts,
potentially creating mislabeled data or requiring re-curation to
replace each use of an A species name with the appropriate B species
name.

It may be that a record in B having a name not occurring in A is
intended to indicate a new concept, call it Y, split off from a
concept in A, call that one X.  Sometimes Y can be connected with X
via a synonym or subspecies in A, but if not there simply isn't enough
information in the checklists to permit the inference that 
that the B name is the result of a 'split', i.e. X < Y.

(In the future maybe we can come up with a good way to add this
information so that the checklist comparison can be complete.  For
example, the A checklist could be amended with the new B name
as a heterotypic synonym for the A name.)


[TBD: the report should give some indication of name changes; otherwise
they are hard to detect on a quick scan.  Maybe a separate
column with 'rename', 'lump', 'split' information.]

Example (excerpt of a larger comparison):

| | A taxon id | A taxon name | B species that intersect | LUB in B | exemplar ids |
|---|---|---|---|---|---|
| (1) | 35492802 | Platyrrhinus lineatus |  > 4JY2M Platyrrhinus lineatus; >< 4JY2R Platyrrhinus umbratus |  < 6S3Q Platyrrhinus | 3935;3936 |
| (2) | 35504725 | Platyrrhinus lineatus nigellus |  <= 4JY2R Platyrrhinus umbratus |  = 855Z5 Platyrrhinus lineatus nigellus* | 3935 |
| (3) | 35504048 | Platyrrhinus lineatus lineatus |  = 4JY2M Platyrrhinus lineatus |  = 4JY2M Platyrrhinus lineatus | 3936 |
| (4) | 35492801 | Platyrrhinus infuscus |  = 4JY2L Platyrrhinus infuscus |  = 4JY2L Platyrrhinus infuscus | 5577 |
| (5) | 35492805 | Platyrrhinus vittatus |  = 4JY2S Platyrrhinus vittatus |  = 4JY2S Platyrrhinus vittatus | 3937 |
| (6) | 35492804 | Platyrrhinus umbratus |  > 874KL Platyrrhinus aquilus; >< 4JY2R Platyrrhinus umbratus |  < 6S3Q Platyrrhinus | 3938;3939;3940 |
| (7) | 35505347 | Platyrrhinus umbratus aquilus |  = 874KL Platyrrhinus aquilus |  = 874KL Platyrrhinus aquilus | 3938 |
| (8) | 35504727 | Platyrrhinus umbratus oratus |  <= 4JY2R Platyrrhinus umbratus |  = 855Z8 Platyrrhinus umbratus oratus* | 3939 |
| (9) | 35504049 | Platyrrhinus umbratus umbratus |  < 4JY2R Platyrrhinus umbratus |  < 4JY2R Platyrrhinus umbratus | 3940 |

(1) _P. lineatus_ sec. A is not in the B checklist, but
it fully contains _P. lineatus_ sec. B, and
it contains some of _P. umbratus_ sec. B (it overlaps (><) but does not contain it).  The nearest (smallest) B 
concept covering all of _P. lineatus_ sec. A is the genus _Platyrrhinus_.

(2) _P. lineatus nigellus_ sec. A is strictly contained in
_P. umbratus_ sec. B, i.e. it has been moved out of _P. lineatus_.

(3) _P. lineatus lineatus_ sec. A promoted to species, i.e. its name in B is _P. lineatus lineatus_.

(4), (5) Species carried over, same concept in both checklists

(6) _P. umbratus_ sec. A (the concept) is not in the B checklist, but is represented by
_P. aquilus_ sec. B (which it contains) and by part of _P. umbratus_ sec. B 
(that is, it overlaps _P. umbratus_ sec. B without containing it completely).

(7) Subspecies _aquilas_ promoted to species; a clerical change in how the concept is named.

(8) _P. umbratus oratus_ sec. A doesn't have its own record in B except as a synonym.
It has been lumped into _P. umbratus_ sec. B (which, remember, differs from
_P. umbratus_ sec. A).  The B checklist has a synonym record for it.k

(9) _P. umbratus umbratus_ sec. A has been lumped into
_P. umbratus_ sec. B, and there is no synonym record in B.  
(Personally, I'm of the opinion that when a checklist is revised. the revision
should always have synonym records deprecated names in previous checklists.  But
this is not always the case.)

Sample output: [col-19-23-report.csv](col-19-23-report.csv)


### ncbi_to_dwc

Converts an NCBI taxdmp .zip file to DwC.  For example, suppose we fetch
a zip file and unzip it, e.g.:

    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2015-05-01.zip
    unzip -d dump taxdmp_2015-05-01.zip

Then we can convert the dump directory to a single DwC .csv file:

    src/ncbi_to_dwc.py dump >taxdmp_2015-05-01.csv

Columns in the DwC output:
`taxonID`, `NCBI Taxonomy ID`, `parentNameUsageID`, `taxonRank`,
`acceptedNameUsageID`, `scientificName`, `canonicalName`,
`taxonomicStatus`, `nomenclaturalStatus`

### project

`project` (emphasis on second syllable) drops columns from the input.
With `--keep`, it drops all columns other than those specified:

    src/project.py --keep a,b,c < foo.csv > less-foo.csv

With `--drop`, it drops particular columns, keeping all the rest:

    src/project.py --drop a,b,c < foo.csv > less-foo.csv


### subset

`subset` generates a subset of the rows of a given file, starting from
a specified root.

    src/subset.py --root 40674 < all.csv > some.csv

It assumes the usual Darwin Core hierarchical structure, given by
these columns: 

 - `taxonID` identifies a row (and a taxon, which is probably an extension)
 - `parentNameUsageID` indicates the record for the parent in the hierarchy
 - `acceptedNameUsageID` indicates an accepted name that might replace an
   non-accepted one (although in fact there is more meaning to this attribute).
 - `taxonomicStatus` is `accepted` for an accepted (non-synonym) row, 
   `synonym` for a synonym row

You can specify the root using its `canonicalName` it that's unique:

    src/subset.py --root Mammalia < all.csv > some.csv

### sortcsv

`sortcsv` sorts standard input by the contents of the given `--pk`
(pk = primary key) column, and writes the result to standard output.

    src/sortscv.py --pk taxonID < foo.csv > foo-sorted.csv


### newick

An extremely rudimentary Newick notation parser.  

The following accepts a Newick string on the command line and emits a
CSV table with basic Darwin Core columns:

    $ src/newick.py --newick "(a,b)c"
    taxonID,canonicalName,parentNameUsageID,acceptedNameUsageID,taxonomicStatus
    2,a,1,,accepted
    3,b,1,,accepted
    1,c,,,accepted
    $ 

The following reads CSV (with header row) from standard input and
writes a Newick string to standard output:

    $ cat <<EOF | src/newick.py
    taxonID,canonicalName,parentNameUsageID,acceptedNameUsageID,taxonomicStatus
    2,a,1,,accepted
    3,b,1,,accepted
    1,c,,,accepted
    EOF
    (a,b)c
    $ 

Branch length syntax is not handled.  Newick escape sequences are not handled.

This program supports an idiosyncratic syntax for synonyms (useful
since the main use for this feature is testing): an asterisk `*`
suffixed to a name says that the name is to be considered a synonym
(i.e. not accepted).


## Typical pipeline

A typical processing pipeline to generate a "align" style report would be:

 1. Obtain two .zip files, DwCAs for the two checklists
 1. Find the Taxon file names in the DwCAs using `find_taxa.py`
 1. (Optional: use `subset.py` to extract subtree of interest)
 1. Prepare checklists for further processing using `clean.py`
 1. Use `gnparse` to obtain stemming and other information useful to `exemplar.py`.
    This requires multiple steps.  For each checklist:
     1. Prepare list of scientific names by applying `extract_names.py` to checklist
     1. Run `gnparse` on that list
     1. Fold the `gnparse` output into checklist with `use_gnparse.py`
 1. Run `exemplar.py` on the checklists to obtain file `exemplars.csv`
 1. Apply `align.py` to the checklists and to `exemplars.csv` to obtain report

