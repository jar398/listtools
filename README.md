# listtools
Tools for manipulating lists of things, and taxonomic checklists in particular.

## Overview

The 'tools' or 'scripts' in this repository manipulate tabular data
in the form of CSV files.  Most of them act as filters and can be
invoked from the Unix shell.  They can be sequenced into pipelines
using '|' if desired.

I have found it convenient to run the tools by using `make` applied to
a suitable 'makefile'.  This way, intermediate objects can be cached
in the local file system, and regenerated only when an input required
for the production changes.

For my testing I use python 3, `bash`, and GNU `make` exclusively.

### Terminology

When I speak of a Darwin core (DwC) file I take this to mean (for
purposes of these tools) a CSV file where each row describes a name
(such as a Linnaean binomial) or a use of a name, in the context of
all the other rows in the same file, and in the context of whatever we
know about the origin of the file itself.  I am calling such things
(names or uses of names) "usages" roughly following "taxon name usage"
or TNU, although I don't pretend that "usage" here is quite as
rigorous a notion as TNU.

The primary key in such files is `taxonID` which is a bit confusing
because the identifier for a usage row doesn't always identify a taxon
- it identifies the way a name is used, and that use may imply some
particular taxon, or not.  Multiple "usages" may correspond to the
same taxon, so if we did have a taxon identifier, it would not
identify a usage.  On the other hand, while a usage identifier (a
`taxonID`) always identifies a usage, it may not be specific enough to
identify a particular taxon.

Other Darwin Core columns containing usage identifiers have column
names that contain 'usage' a morpheme, e.g. `parentNameUsageID`.

Sometimes a given name is used in multiple ways in the same file
(homonyms, hemihomonyms, etc.), in which case the usages are
distinguished by their usage identifiers.

In the EOL internals, usage rows (records) are called "nodes", but I
avoid this word due to its various conflicting and restricting
meanings, such as graph database nodes in Neo4J.

## The tools

### start

Any Darwin Core data source should be run through src/start.py first,
for some modest validation and DwC-specific cleaning.

    src/start.py --pk taxonID <foo-raw.tsv >foo.csv

`--pk` specifies the primary key used for this table.

### sortcsv

    src/sortscv.py --pk taxonID <foo.csv >foo-sorted.csv

This should be clear...

### project

The following drops all columns other than those specified:

    src/project.py --keep a,b,c <foo.csv >less-foo.csv

The following drops particular columns, keeping all the rest:

    src/project.py --drop a,b,c <foo.csv >less-foo.csv

### diff

`diff` compares two files, which we can call A and B, matching rows of
one to rows of the other, and generating a "delta" which we might call
B-A.  A delta is an annotated report listing all unmatched and changed
rows, and not listing exact row matches.

    src/diff.py <a.csv --target b.csv \
                       --pk a
                       --index x,y,z
                       --managed y,z,d,e,f

Records are matched by the values in their 'index' columns, which
should be given in priority order.  That is, if the index columns are
x, y, and z, then at attempt to match records is made first on the
values in the x column.  If one or both values is missing, or if the
match is ambiguous, then the y column is consulted, and so on.

The output contains only the managed columns, and matched rows are
considered changed only if one of the managed columns differs.

`--pk` specifies the primary key column for both inputs.

### apply

Applies a sorted delta B-A to sorted file A (which typically would be
the A file from which the delta was generated), generating a file B′.
B′ will be projection of B to the given 'managed' columns, with
perhaps the rows in a different order.

    $P/apply.py --delta ba.delta --pk taxonID \
                < a.csv > b2.csv

### hierarchy

This is specific to EOL.  It applies a usage id to page id mapping
(provided as an argument) to a file full of usages to generate a file
with one row per page, giving the parent of each page.

The resulting page list can be subjected to a

### ncbi_to_dwc

Converts an NCBI taxdmp .zip file to DwC.  For example, suppose we fetch
a zip file and unzip it, e.g.:

    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2015-05-01.zip
    unzip -d dump taxdmp_2015-05-01.zip

Then we can convert the dump directory to a single DwC .csv file:

    src/ncbi_to_dwc.py dump >taxdmp_2015-05-01.csv

### scatter

Given a delta, generate a directory containing one file for each of
the three kinds of record in the delta (add, udpate, remove).

    src/scatter.py <delta.csv --dest delta

writes `delta/add.csv`, `delta/update.csv`, `delta/remove.csv`.  These
can be fed to appropriate database commands to modify a table in the
database.

### subset

`subset` generates a subset of the rows of a given file assuming a
usual Darwin Core hierarchical structure: `taxonID` to give a way to
name a usage row, `parentNameUsageID` to indicate the usage row for
the parent in the hierarchy, and `acceptedNameUsageID` to indicate an
accepted usage that can replace an non-accepted one.

    src/subset.py --root 40674 <all.csv >some.csv

### util

This is not invoked on its own, but just contains some useful code
shared among the various tools.


## Configuration for EOL

For EOL, the list tools are intended to be used in conjunction with
`plotter`, which provides access to the EOL content repositories.

Configure `plotter` according to its documentation.  Then, assuming
the `plotter` repository is cloned in directory `../plotter`, do the
following:

    ln -sfT ../plotter/config config

If it's elsewhere, adjust the above command, and also update
`Makefile` where it mentions `plotter`.
