# List tools user guide

## Overview

The 'tools' or 'scripts' or 'programs' in this repository manipulate
tabular data in the form of CSV files.  Most of them act as filters
and can be invoked from the Unix shell.  They can be sequenced into
pipelines using '|' if desired.  CSV files uniformly have a single
header row giving column names.

I have found it convenient to run the tools by using `make` applied to
a suitable `makefile'.  This way, intermediate objects can be cached
in the local file system, avoiding expensive regeneration steps when
inputs have not changed.

For a complete example, see [example.md](example.md).

## Installation
<a name="installation"></a>

The tools require no particular installation step.  They can be
invoked from the shell out of a clone of this `git` repository.

Python 3 is assumed.

Install the python3 `regex` module.

Also required is `gnparser`, which has download instructions
[here](https://github.com/gnames/gnparser#installation).


## Semantics

### Individuals

I'll use the neutral word 'individual' for the entities that we are
classifying; it's meant to subsume various terms in taxonomy and
biodiversity informatics such as 'organism', 'specimen', and
'observation'.

### Taxon concepts

Each record of each checklist has an associated "taxon concept"
determined by the fields of the record with the checklist as context.

Say that an individual "falls under" a taxon concept meaning that the
individual is classified under, or belongs to, that taxon concept.

We may not know much about the taxon concept, and it would be
challenging to find out (it would usually require a literature
search).  However, that doesn't prevent us from reasoning about it
using information in the checklist.

The individuals falling under a taxon concept are its 'extension'.

According to the way the word 'taxon' is typically used, a taxon is
tied to a name and its extension changes following speech acts such as
'redescription' or 'lumping' or 'splitting'.  By contrast, a taxon
concept is rigid: whether an individual falls under a taxon concept
does not change over time and is insensitive to what happens in the
taxonomic literature.

### Type specimens

Each record of each checklist has a name ('scientific' if it includes authority 
information; 'canonical' if not),
and the name has a designated individual (its type specimen,
or just 'type') according to
the rules of biological nomenclature.

Therefore, each record has an associated type, via its name, in context.

The type associated with a record falls under
the taxon concept associated with that record.

We generally don't know much about a type - we don't know where it is
housed, or who collected it, and so on.  But that doesn't matter since
algorithmically we'll only be concerned with membership of types in
taxon concepts: 
 - whether or not any particular type belongs to any
   given taxon concept (as represented by a record from a checklist), 
 - for each type which taxon concepts contain it,
 - and for each taxon concept which types it contains.


### Exemplars

For each type specimen, we can form the group of records, some in A and
some in B, with which the type specimen is associated.  If the group contains at least
one record from each checklist, call the type specimen an "exemplar"
with respect to checklists A and B.

We take sets of exemplars to be computable approximations to (unknown)
taxon concepts, in the sense that if exemplar set S = E(C) = {e: e is
an exemplar in C} for concept C, and similarly T = E(D), then an RCC-5
relationship (see below) between concepts induces the same RCC-5 relationship
between the exemplar sets.  The opposite also holds, except in some
cases where S = T, in which case contextual information can be used to
infer a distinction between C and D [work in progress].

Exemplars play a role similar to that played by protonyms in the
Pyle/Remsen formulation.

### RCC-5

RCC-5 (region connection calculus) is a simple system for classifying
relationships between 'regions'.  What constitutes a region depends on
how RCC-5 is being used - we could be talking sets, or geographic
regions - but in our case we're talking about taxon concepts.  The
RCC-5 relationships, exactly one of which holds for any
two regions, are
 * A = B: A and B are the same; the same individuals fall under both
 * A < B: A is contained in B but isn't the same as B; 
   the individuals in A all fall under B but not vice versa
 * A > B: same as B < A
 * A ! B: A and B are disjoint; no individual falls under both
 * A >< B: A and B overlap but neither is contained in the other.
   There are individuals falling under A and B, and under just A,
   and under just B.
   By convention we'll refer to this situation as 'overlap'.

In cases where the precise RCC-5 relationship isn't known, we can also write
 * A <= B: A < B or A = B
 * A >= B: A > B or A = B
 * A not! B: it is not the case that A ! B (i.e. at least one individual falls under both)
 * A ? B: RCC-5 relationship is unknown


## The tools

To learn about a given tool, please supplement what follows with
whatever you learn from using the `--help` option, e.g.

    src/clean.py --help

Some tools operate on arbitrary CSV files, while some assume they're
working with Darwin Core files.

### find_taxa

Locates the Darwin Core taxon file within a .zip file.  Let `A.dump/`
be a directory containing the files resulting from unzipping the .zip
file.  Then:

    src/find_taxa.py A.dump

writes the name of the taxon file within the .zip file, e.g. 

    A.dump/Taxon.tsv

It's usually clear by inspection which file is the taxon file, but
spelling details vary.  This command is intended for use in scripts.

The taxon file is suitable as input to `clean.py`, see below.

TBD: This currently operates by examining file names heuristically.
But it really ought to look in the meta.xml file, which provides the
taxon file name explicitly.

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


### extract_names

This extracts `scientificName`s in a form suitable for consumption by `gnparser`.
If there is no `scientificName` then the `canonicalName` is extracted.

    src/extract_names < work/A-clean.csv > work/A-names.txt
    gnparser -s < work/A-names.txt > work/A-gnparsed.csv

### use_gnparse

This consumes the output of `gnparser` and combines it with the table
that was the input to `extract_names`, enriching the table with the addition of 
columns from the `gnparser` output.

    src/use_gnparse < work/A-gnparsed.csv > work/A.csv

### exemplar
<a name="exemplar"></a>

This is a heuristic name matcher.  Its purpose is to find groups of
taxon concepts that share a type specimen.
It understands changes in genus, changes
in gender, and other vagaries of the nomenclatural system.

Given input checklists A and B in Darwin Core form, finds groups of A
and B records (at least one of each) such that all members of each
group share a type specimen.

 * `--A` filename  - the A checklist.
 * `--B` filename  - the B checklist.

Writes a file giving exemplar group membership to standard output.

    src/exemplar --A work/A.csv --B work/B.csv > work/AB-exemplars.csv

There is one output row for each A or B checklist record whose associated
type specimen is an exemplar.
 - `checklist`: 0 for the A checklist, 1 for B
 - `taxonID`: the checklist row for a taxon concept
 - `exemplar id`: identifies an exemplar, locally to this analysis (not global).

The meaning of an output row is that the type specimen of the
indicated taxon concept is the exemplar identified by `exemplar id`.
(Of course the same exemplar can also be the type specimen of other
taxon concepts.)

## plugin
<a name="plugin"></a>

    src/plugin.py --A work/A.csv --B work/B.csv --exemplars work/AB-exemplars.csv

This writes an analysis report to standard output.

If you're doing regression analysis, think of A as the 'baseline' 
checklist and B as 'proposed successor'.  (This is not the only use case.)

If `--exemplars` is not given, the exemplars are computed just as the
`exemplar.py` command would.

The output (to standard output) has these columns (subject to change):
 - `A taxon id` - The taxon id of an A row
 - `A taxon name` - The canonicalName of that A record (for human consumption)
 - `B species that intersect` - 
   If the A record indicates rank 'species', this is a semicolon-separated
   list of relationship/id/name for
   B taxon concepts with rank 'species' that are inferred to intersect the A 
   taxon concept.
   The relationship is the RCC-5
   relationship of the A concept to the B concept, the id
   is the taxon id of the B concept's record, and the name is
   canonical name from the B record.
   A value of `-` means there may be intersecting concepts but the list was not computed 
   because the A row was not for a species.
 - `LUB in B` - 
   The relationship/id/name of the A concept to its least upper bound (LUB)
   in the B checklist.  The least upper bound is the smallest B concept
   that contains the A concept.  
   If the A
   and B names are accepted, the A concept is either the same (RCC-5 =)
   as the B concept or smaller (RCC-5 <) than it.  For synonyms it may
   be hard to tell what the precise relationship is so it <= or ? will show.
 - `exemplar ids` - 
   If the A record is indicates a species, a list of exemplar ids for the exemplars
   belonging to that species, otherwise `-`

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

A taxon id | A taxon name | B species that intersect | LUB in B | exemplar ids |
---|---|---|---|---|
35492802 | Platyrrhinus lineatus |  > 4JY2M Platyrrhinus lineatus; >< 4JY2R Platyrrhinus umbratus |  < 6S3Q Platyrrhinus | 15672;15673;15674;15675 |
35504725 | Platyrrhinus lineatus nigellus | - |  < 4JY2R Platyrrhinus umbratus | -
35504048 | Platyrrhinus lineatus lineatus | - |  = 4JY2M Platyrrhinus lineatus | -
35492801 | Platyrrhinus infuscus |  = 4JY2L Platyrrhinus infuscus |  = 4JY2L Platyrrhinus infuscus | 15676
35492805 | Platyrrhinus vittatus |  = 4JY2S Platyrrhinus vittatus |  = 4JY2S Platyrrhinus vittatus | 15677
35492804 | Platyrrhinus umbratus |  > 874KL Platyrrhinus aquilus; >< 4JY2R Platyrrhinus umbratus |  < 6S3Q Platyrrhinus | 15678;15679;15680;15681;15682;15683
35505347 | Platyrrhinus umbratus aquilus | - |  = 874KL Platyrrhinus aquilus | -
35504727 | Platyrrhinus umbratus oratus | - |  < 4JY2R Platyrrhinus umbratus | -
35504049 | Platyrrhinus umbratus umbratus | - |  < 4JY2R Platyrrhinus umbratus | 15678




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

    src/subset.py --root 40674 < work/all.csv > work/some.csv

It assumes the usual Darwin Core hierarchical structure, given by
these columns: 

 - `taxonID` identifies a row (and a taxon, which is probably an extension)
 - `parentNameUsageID` indicates the record for the parent in the hierarchy
 - `acceptedNameUsageID` indicates an accepted name that might replace an
   non-accepted one (although in fact there is more meaning to this attribute).
 - `taxonomicStatus` is `accepted` for an accepted (non-synonym) row, 
   `synonym` for a synonym row

You can specify the root using its `canonicalName` it that's unique:

    src/subset.py --root Mammalia < work/all.csv > work/some.csv

### sortcsv

`sortcsv` sorts standard input by the contents of the given `--pk`
(pk = primary key) column, and writes the result to standard output.

    src/sortscv.py --pk taxonID <foo.csv >foo-sorted.csv


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
