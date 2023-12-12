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

Say that an individual "falls under" a taxon concept to mean that the
individual is classified under, or belongs to, that taxon concept, or
the taxon concept applies to the individual.

I'll write '_Hyla_ sec. C' to denote the taxon concept described in
source C under the name '_Hyla_'.  (C might be a checklist.)

We may not know much about the taxon concept, and it could be
challenging to find out (it would usually require a literature
search).  However, that doesn't prevent us from reasoning about it
using information in the checklist.

The individuals falling under a taxon concept are its 'extension'.

According to the way the word 'taxon' is typically used, a taxon is
tied to a name, and its extension changes following speech acts such as
'redescription' or 'lumping' or 'splitting'.  By contrast, a taxon
concept is rigid: whether an individual falls under a taxon concept
does not change over time and is insensitive to what happens in the
taxonomic literature.

### Exemplars

An individual is an _exemplar_ for a comparison of checklists A and B
if taxon concepts in both A and B are known such that the individual
falls under both of the concepts.  The individual is proof that the
taxon concepts intersect.

It is most useful if each taxon concept is as small or as 'fine' as it
can be in its checklist; that is, it contains as few concepts as
possible from that checklist.  This is not a requirement of the
analysis, but it is a way to improve its performance.

Exemplar is a semantic ideal rather than an operational notion.  To
make the exemplar idea practical when we don't have direct information
about individuals and the taxon concepts under which they fall, we can
use type specimens as the individuals in the analysis, and take those
individuals to be exemplars when they are known in both checklists.
This is explained below under the `exemplar.py` command.

The pragmatic representation of an exemplar, when exemplars are based
on type specimens, is therefore a set of 'matched' records containing
at least one record from each checklist.

### Exemplar sets

We take sets of exemplars to be computable approximations to (unknown)
taxon concepts, in the sense that if exemplar set 

&nbsp;&nbsp;&nbsp;&nbsp;S = E(C) = {e: e is an exemplar in C}

for concept C, and similarly

&nbsp;&nbsp;&nbsp;&nbsp;T = E(D),

then an RCC-5 relationship (see below) between concepts induces the
sames RCC-5 relationship between the exemplar sets.  The opposite also
holds, except in some cases where S = T, in which case contextual
information can be used to infer a distinction between C and D [work
in progress].

### RCC-5

RCC-5 (region connection calculus) is a simple system for classifying
relationships between 'regions'.  What constitutes a region depends on
how RCC-5 is being used - we could be talking sets, or geographic
regions - but in our case we're talking about taxon concepts.  The
RCC-5 relationships, exactly one of which holds for any
two regions, are
 * A = B: A and B are the same, at least inasmuch as the same individuals fall under both
 * A < B: A is contained in B but isn't the same as B; 
   the individuals that fall under A all fall under B but not vice versa
 * A > B: same as B < A
 * A ! B: A and B are disjoint; no individual falls under both
 * A >< B: A and B overlap but neither is contained in the other.
   There are individuals falling under A and B, and under just A,
   and under just B.
   By convention we'll refer to this situation as 'overlap'.

In cases where the precise RCC-5 relationship isn't known, we can also write
 * A <= B: A < B or A = B
 * A >= B: A > B or A = B
 * A not! B: A and B intersect; are not disjoint;
   it is not the case that A ! B;
   at least one individual falls under both
 * A ? B: RCC-5 relationship is unknown


## File formats

### Checklist file format

A checklist is given as a tabular text file.  If the extension is
`.csv` it is assumed to be CSV (comma-separated); otherwise it's
assumed to be TSV (tab-separated).  Each row is a 'record'.

A checklist should use [Darwin Core
Taxon](https://dwc.tdwg.org/terms/#taxon) column headings for
information of special significance to these commands.  Other columns
can be included and will be ignored or passed through as appropriate.

A prefix `dwc:` on column headings is optional.

`canonicalName` is an additional column that is not from Darwin Core
but is important.

`managed_id` is a column used internally to the tools but it not
relevant to most users.  (It is used for distinguishing the case of
aligning checklists with a 'managed' identifier space (e.g. versions
of GBIF or NCBI Taxonomy) from checklists that could have accidental
identifier collisions.)


Here are Darwin Core headings used by one or more of the tools.

 - `taxonID`: the record's primary key, uniquely specifying a record within a checklist
 - `scientificName`: the full taxonomic name with authorship information
 - `canonicalName`: the taxonomic name without authorship
 - `scientificNameAuthorship`: the authorship (e.g. Smith, 1825)
 - `namePublishedInYear`: year of publication; part of authorship (e.g. 1825)
 - `taxonRank`: `species`, `subspecies`, `genus`, and so on
 - `taxonomicStatus`: `accepted`, `synonym`, or any of a few variants
   of this [document]
 - `nomenclaturalStatus`: 
 - `taxonRemarksCitation`: 

One of `scientificName` or `canonicalName` must be given.


### Exemplar file format

Some commands read or write files that specify exemplars.  An exemplar
file can be computed by the `exemplar.py` command or provided
independently if it can be made in some other way.  This section
should be of interest if you don't want or don't need to use
`exemplar.py`.

An exemplar file has one output row for each inference (or any kind of
statement) that an exemplar falls under some taxon concept.

 - `checklist`: 0 for the A checklist, 1 for B
 - `taxonID`: the checklist record for a taxon concept
 - `exemplar id`: identifies an exemplar, locally to this file
   (not globally)

By construction, each exemplar id will have at least one exemplar file
row giving an A record for a taxon concept it falls under, and one
giving a B record for a taxon concept it falls under.

If an individual E belongs to taxon concept T, then E is trivially
seen to fall under all ancestors of T.  Therefore the exemplar file
doesn't need to include these redundant statements about ancestors.
By the same token, it is desirable to find the most specific taxon
concepts containing E (as described above).



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

    src/extract_names.py < work/A-clean.csv > work/A-names.txt
    gnparser -s < work/A-names.txt > work/A-gnparsed.csv

### use_gnparse

This consumes the output of `gnparser` and combines it with the table
that was the input to `extract_names`, enriching the table with the addition of 
new columns copied from the `gnparser` output.

    src/use_gnparse.py < work/A-gnparsed.csv > work/A.csv

### exemplar
<a name="exemplar"></a>

Writes a file giving memberships of exemplars in taxon concepts to
standard output.

 * `--A` filename  - the A checklist.
 * `--B` filename  - the B checklist.

For example,

    src/exemplar.py --A work/A.csv --B work/B.csv > work/AB-exemplars.csv

To do this, noting that each record in each checklist has a taxonomic
name, we posit that each record has a type specimen (the type specimen
for its name).  The type specimen is then an individual that, if found
in taxon concepts in both checklists, can be taken to be an exemplar.

To determine that a type specimen is an exemplar we must find a record
in the other checklist that has the same type specimen.  That is, the
records must be 'matched'.  (Note that it is records, not their taxon
concepts, that are matched.)  The matcher understands changes in
genus, changes in gender, and other vagaries of the nomenclatural
system, such as the possibility of collisions (homonyms).

The meaning of an output row is that the individual (exemplar)
identified by `exemplar id` falls under the taxon concept / record identified
by `taxonID` in the indicated checklist.  (Of course the same exemplar
can also fall under other taxon concepts, in particular taxon concepts
in the other checklist.)

Sample output: [col-19-23-exemplars.csv](col-19-23-exemplars.csv)

When the exemplars are type specimens, they play a role similar
to that played by protonyms in the Pyle/Remsen formulation.


### plugin
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
35492802 | Platyrrhinus lineatus |  > 4JY2M Platyrrhinus lineatus; >< 4JY2R Platyrrhinus umbratus |  < 6S3Q Platyrrhinus | 3935;3936
35504725 | Platyrrhinus lineatus nigellus |  <= 4JY2R Platyrrhinus umbratus |  = 855Z5 Platyrrhinus lineatus nigellus* | 3935
35504048 | Platyrrhinus lineatus lineatus |  = 4JY2M Platyrrhinus lineatus |  = 4JY2M Platyrrhinus lineatus | 3936
35492801 | Platyrrhinus infuscus |  = 4JY2L Platyrrhinus infuscus |  = 4JY2L Platyrrhinus infuscus | 5577
35492805 | Platyrrhinus vittatus |  = 4JY2S Platyrrhinus vittatus |  = 4JY2S Platyrrhinus vittatus | 3937
35492804 | Platyrrhinus umbratus |  > 874KL Platyrrhinus aquilus; >< 4JY2R Platyrrhinus umbratus |  < 6S3Q Platyrrhinus | 3938;3939;3940
35505347 | Platyrrhinus umbratus aquilus |  = 874KL Platyrrhinus aquilus |  = 874KL Platyrrhinus aquilus | 3938
35504727 | Platyrrhinus umbratus oratus |  <= 4JY2R Platyrrhinus umbratus |  = 855Z8 Platyrrhinus umbratus oratus* | 3939
35504049 | Platyrrhinus umbratus umbratus |  < 4JY2R Platyrrhinus umbratus |  < 4JY2R Platyrrhinus umbratus | 3940

Row 1: _P. lineatus_ sec. A is not in the B checklist, but
it fully contains _P. lineatus_ sec. B, and
it contains some of _P. umbratus_ sec. B (it overlaps (><) but does not contain it).  The nearest (smallest) B 
concept covering all of _P. lineatus_ sec. A is the genus _Platyrrhinus_.

Row 2: _P. lineatus nigellus_ sec. A is strictly contained in
_P. umbratus_ sec. B, i.e. it has been moved out of _P. lineatus_.

Row 3: _P. lineatus lineatus_ sec. A promoted to species, i.e. its name in B is _P. lineatus lineatus_.

Rows 4, 5: Species carried over unchanged.

Row 6: _P. umbratus_ sec. A (the concept) is not in the B checklist, but is represented by
_P. aquilus_ sec. B (which it contains) and by part of _P. umbratus_ sec. B 
(that is, it overlaps _P. umbratus_ sec. B without containing it completely).

Row 7: Subspecies _aquilas_ promoted to species; a clerical change in how the concept is named.

Row 8: _P. umbratus oratus_ sec. A doesn't have its own record in B except as a synonym.
It has been lumped into _P. umbratus_ sec. B (which, remember, differs from
_P. umbratus_ sec. A).  The B checklist has a synonym record for it.k

Row 9: _P. umbratus umbratus_ sec. A has been lumped into
_P. umbratus_ sec. B, and there is no synonym record in B.  
(Personally, I'm of the opinion that when a checklist is revised. the revision
should always have synonym records for names previous checklists.)

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
