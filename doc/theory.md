
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

