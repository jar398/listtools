# List tools user guide

## Overview

The 'tools' or 'scripts' or 'programs' in this repository manipulate
tabular data in the form of CSV files.  Most of them act as filters
and can be invoked from the Unix shell.  They can be sequenced into
pipelines using '|' if desired.  CSV files uniformly have a single
header row giving column names.

I have found it convenient to run the tools by using `make` applied to
a suitable 'makefile'.  This way, intermediate objects can be cached
in the local file system, avoiding expensive regeneration steps when
inputs have not changed.

For my testing I use python 3, `bash`, and GNU `make` exclusively, all
running on Debian Linux.

## Terminology - rows, taxa, extensions, interpretation

Many of the tools are completely generic over tabular data, but a few
are specific to biodiversity information in the form of "Darwin Core"
files.

When I speak of a Darwin core (DwC) file I take this to mean (for
purposes of these tools) a CSV file where each record (row) has
information connected to a taxon.  More precisely, I take a row to
refer to what I'd call an _extension_ (of a taxon description).  An
extension is simply the set of organisms/specimens/observations as
described or circumscribed or referenced somewhere, perhaps in a
database.  Each row is itself a little taxon description, and may
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
of the file itself.

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


## Installation
<a name="installation"/>

Python 3, `bash`, and Gnu `make` are all assumed.

Install the python3 `regex` module (although it's easy to revert to
the built-in `re` if you prefer not to).

One of the `Makefile` rules uses `gnparser`, which has download
instructions [here](https://github.com/gnames/gnparser#installation).
`Makefile` assumes that `gnparser` can be located via your shell's `PATH`.
It is easy to remove this dependency if desired, with only a mild
degradation in effectiveness.


## The Makefile

There is no need to use the Makefile since there is no build step for
this project.  However, it contains a variety of canned pipelines for
creating artifacts that are interesting and illustrative.  For example
`make ncbi-report` downloads two versions of the NCBI (Genbank)
taxonomy and produces comparison reports.  It may be helpful to consult
the Makefile for examples of use of the various tools and ideas on how
to string them together.

Read the comments in the [`Makefile`](../src/Makefile) for other things to try.


## The tools

### start

Any Darwin Core data source should be run through src/start.py first,
for some modest validation and DwC-specific cleaning.

    src/start.py --pk taxonID <foo-raw.tsv >foo.csv

`--pk` specifies the column containing the table's primary keys.  If
there is no such column, one is added, and if a primary key is
missing, a new one is generated by hashing the contents of the row.

Reused putatively-primary keys are detected and flagged.

Data cleaning is performed for the following columns:

 - `canonicalName`, `scientificName` - if a scientific name (with
   authority) is found in the canonicalName column, and the
   `scientificName` column is empty, then the scientific name is moved
   to the `scientificName` column (and so on)
 - `acceptedNameUsageID` - if empty, copy the `taxonID`
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
column.  This will then be used for matching operations.


### sortcsv

`sortcsv` sorts standard input by the contents of the given `--pk`
column, and writes the result to standard output.

    src/sortscv.py --pk taxonID <foo.csv >foo-sorted.csv


### project

`project` (emphasis on second syllable) drops columns from the input.
With `--keep`, it drops all columns other than those specified:

    src/project.py --keep a,b,c <foo.csv >less-foo.csv

With `--drop`, it drops particular columns, keeping all the rest:

    src/project.py --drop a,b,c <foo.csv >less-foo.csv


### extract_names

This extracts scientific names in a form suitable for consumption by `gnparse`.


### use_gnparse

This consumes the output of `gnparse` and combines it with the table
that it just processed, enriching the table with the addition of extra
columns that can be used in further analysis.

 * `canonicalName` is added if it is not already present
 * `canonicalStem` is added (see the `gnparse` documentation)
 * `year` is added
 * `type` is added - this is a special string used for matching, and has 
   the form `TS|yyyy|eeeee|aaa` with components generated from
   parsing the scientific name.  `eeeee` is the _last_ epithet,
   either subspecific, specific, or generic.  `aaa` is supposed to be
   the last name of the first author (but I bet it's sometimes
   wrong). `TS` is supposed to remind us that the string is supposed
   to be associated with the type specimen or type series, and not any
   particular taxon containing them.


### match_records

Given input checklists A and B, finds the best unique mutual matches
between the A records and the B records.

 * `--A` filename  - the A input.
 * `--B` filename  - the B input.
 * `--pk` K - specify primary key; default `taxonID`
 * `--index` columns names - lists columns to be used for comparison,
      in priority order.

Standard input is the A input, but it can be specified as the B input
by giving `-` as the filename.

The output (to standard output) has these columns:

 - `match_id` - the id of an A record, or empty if the row corresponds
   to an unmatched B record.
 - `relation` - an RCC5 relationship, but really either `=`, `?` (no
   information), or empty (nothing to compare to).
   (Mistake, that ought to be `relationship`.)
 - `taxonID` - the id of a B record, or empty if the row corresponds
   to an unmatched A record.
 - `basis_of_match` - a textual report of what happened, with
   the deciding column name is a match was made, and a diagnostic
   if no match was made and there was a near miss.

If the primary key is specified as something different, the first
three columns reflect what's given.

Records are matched by the values in their 'index' columns, which
should be given in priority order.  That is, if the index columns are
x, y, and z, then an attempt to match records is made first on the
values in the x column.  If either value is missing, or if the
match is ambiguous, then the y column is consulted, and so on.

In the case of an ambiguity, i.e. multiple A records matching a B
record with the same match score or vice versa, none of the records is
matched.  This fact is recorded in the `remark` column of the delta.
A situation like this should provoke manual review; it could mean
that the inputs are ill-formed, or it could mean the scoring algorithm
simply doesn't have enough information to choose between the matches.

The B input should not contain a `mode` column, and if it does, the
column is ignored.

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

### subset

`subset` generates a subset of the rows of a given file, starting from
a specified root.

    src/subset.py --root 40674 <all.csv >some.csv

It assumes the usual Darwin Core hierarchical structure, given by
these columns: 

 - `taxonID` identifies a row (and a taxon, which is probably an extension)
 - `parentNameUsageID` indicates the record for the parent in the hierarchy
 - `acceptedNameUsageID` indicates an accepted name that might replace an
   non-accepted one (although in fact there is more meaning to this attribute).
 - `taxonomicStatus` is `accepted` for an accepted (non-synonym) row, 
   `synonym` for a synonym row

You can specify the root using its `canonicalName` it that's unique:

    src/subset.py --root Mammalia <all.csv >some.csv

### merge

Replaces `align` TBD

### align

Hierarchy alignment, seeded with record matches (`match_records`,
above).  Work in progress.

 - Input files: A; B; A/B record matches (see `match_records`)
 - Standard output: the 'sum' of A and B.
   Columns:
     - `taxonID` - a unique primary key within the output.  Usually but not 
       always either the A taxonID or the B taxonID.
     - `taxonID_A` - a taxonID from the A input file, or empty
     - `taxonID_B` - a taxonID from the B input file, or empty
     - `parentNameUsageID` - a taxonID in the output, or empty at the root(s)
     - `acceptedNameUsageID` - a taxonID in the output, when the row is a synonym.
     - `canonicalName` - this is mainly for debugging and ease of inspection.
       A `canonicalName` copied from the B input, if taxonID_B is present,
       or from the A input otherwise.
     - `remark` - annotation concerning how the alignment was made, or why it wasn't

Invocation:

    src/align.py < A.csv --target B.csv --matches AB_RM.csv > AB.csv

where `AB_RM.csv` was generated by `match_records`.

### newick

An extremely rudimentary Newick notation parser.  

The following
accepts a Newick string on the command line and emits a CSV table with
columns `taxonID`, `parentNameUsageID`, and `canonicalName`:

    src/newick.py "(a,b)c"

The following reads CSV (with header row) from standard input and
writes a Newick string to standard output:

    src/newick.py

Branch length syntax is not handled.  Newick escape sequences are not handled.

TBD: be able to take a Newick string from standard input.

### util

This module is not invoked on its own.  It contains a few internal
utilities shared among the various tools.


## Record-based delta/diff/patch feature

(These programs are inherited from a previous round of development and
may not work any more!  Note that they operate only on records,
insensitive to hierarchy.)

### delta

`delta` compares two files (call them A and B), matching rows of one
to rows of the other, and generating a "delta" (call it B-A) which
describes the differences between A and B.  A delta is an annotated
report listing all unmatched and changed rows, not including exact row
matches.

    src/delta.py <a.csv --target b.csv \
                        --matches m.csv \
                        --pk y \
                        --index x,y,z \
                        --managed z,y,d,e,f

`--matches` names a file that was generated by `match_records` (or
another tool generating the same kind of file).  If it's omitted then
`match_records` is invoked implicitly.

`--pk` specifies the primary key column for both inputs.

Each row of the delta comes from A only, from B only, or from matched
rows in A and B.  Because of the use of deltas in patching, these
three types of row are given labels `remove`, `add`, and `update` in
the `mode` column of the delta.  The primary key in the delta is taken
from A in the case of `remove` and `update` records.  The primary key
from the B file is given in `new_pk` of the delta.

The output contains only the managed columns (`--managed`), and matched
rows are considered updated only if one of the managed columns
differs.

The matches are done on only a row-by-row basis and are not sensitive
to hierarchy or other sources of meaning.  Hierarchy and
synonym-to-accepted links are treated the same as any other field and
do not require identical contents, meaning that the overall comparison
is not truly sensitive to hierarchy.  Rows may be matched even
when consideration of hierarchy would require them to be interpreted
as distinct taxa.  For hierarchy sensitive comparison see `align`, below.

### apply

Applies a sorted delta, B-A, to a sorted file A (which typically would be
the A file from which the delta was generated), generating a file B′.
B′ will be projection of B to the given 'managed' columns, with
perhaps the rows in a different order.

    src/apply.py --delta ba.delta --pk taxonID \
                 < a.csv > b2.csv

The delta would probably be produced by the delta tool.  It should
have columns `mode` and `new_pk` (see above).

### scatter

Given a delta, generate a directory containing one file for each of
the three kinds of record in the delta (add, udpate, remove).

    src/scatter.py <delta.csv --dest delta

writes `delta/add.csv`, `delta/update.csv`, `delta/remove.csv`.  These
can be fed to appropriate database commands to incrementally update a
database to a new version of a table (intended for EOL mainly).

### hierarchy

This is specific to EOL.  It applies a taxon id (`taxonID`) to 
'page id' mapping to a file full of records to generate a file
with one record per page, giving the parent of each taxon.

The resulting taxon list can be subjected to `delta` and `scatter` to
incrementally update an in-database hierarchy, etc.

Columns expected: `taxonID`, `parentNameUsageID`,
`acceptedNameUsageID`, `taxonomicStatus`


## Reading the output

I use a few conventions for the progress output generated by the tools:

 - Lines beginning with a run of hyphens `--` are intended to be read by the
   general user, and describe what is going on in terms that ought to
   be easily understood
 - Lines beginning with asterisk `*` have to do with the detection of bugs or
   peculiar or confusing situations that might demand some attention
 - Lines beginning with a hash `#` are directed at me, and I don't
   expect anyone else to understand them
 - Other lines come from code that I didn't write, such as `make` or `wc`

