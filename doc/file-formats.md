# File formats

## CSV format

Doc TBD

## Darwin Core format
<a name="DwCA"></a>

A Darwin Core Archive, for present purposes, is a `.zip` file that
contains a Darwin Core Taxon table in TSV or CSV format.  [doc
location TBD]

## Checklist file format

A checklist is given as a CSV file.  TSV files must be converted to
CSV for use most of the tools.  See [`clean`]<a href="guide#clean">.

A checklist should use [Darwin Core
Taxon](https://dwc.tdwg.org/terms/#taxon) column headings for
information of special significance to these commands.  Other columns
can be included and will be ignored or passed through as appropriate.

A prefix `dwc:` on column headings is stripped off.  Yes, I know this
is wrong.  Some input sources use this prefix, some don't.

`canonicalName` is an additional column used formerly by GBIF that is
not from Darwin Core but is important.

`managed_id` is a column used internally to the tools but is not
relevant to most users.  (It is used for distinguishing the case of
aligning checklists with a 'managed' identifier space (e.g. versions
of GBIF or NCBI Taxonomy) from checklists that could have accidental
identifier collisions.)


Here are Darwin Core headings used by one or more of the tools,
post-cleaning.  Usage of these columns by other software may vary.

 - `taxonID`: the record's primary key, uniquely specifying a record within a checklist
 - `scientificName`: the full taxonomic name with authorship information
 - `canonicalName`: the taxonomic name without authorship
 - `scientificNameAuthorship`: the authorship (e.g. `Smith, 1825`).
   Optional but the matcher is more precise if authorship is present
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


## Exemplar file format

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

