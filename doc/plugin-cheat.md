# List tools "plugin" summary 2023-11-10

Jonathan Rees - work in progress

Scenario:
- Start with two checklists, call them A and B.
- If you're doing regression analysis, think of A as the 'baseline' 
  checklist and B as 'proposed successor' (this is not the only use case)
- Checklists must be in Darwin Core format.

Semantics of checklists:
- Each record of each checklist has an associated "taxon concept"
  determined by the fields of the record with the checklist as context.
  We may not know what the taxon concept is and it would be challenging to find 
  out.  However, that doesn't prevent us from reasoning about it
  using information in the checklist.
- Each record of each checklist has a name ('scientific' if it includes authority 
  information; 'canonical' if not),
  and the name is assumed to have an associated type specimen under
  control of the rules of biological nomenclature (again, we may not have any
  details at hand on the type specimen - that's OK).  If it's not
  clear whether two records have the same type specimen or not,
  information beyond just the name may allow this determination.
- Therefore each record has an associated type specimen.
- We know that the type specimen associated with a record falls under
  the taxon concept associated with that record.

'Exemplar' determination:
- Find sets of records from the two checklists such that
    . in each set, all records have the same associated type specimen
    . all records known to have that type specimen are in the set
    . each set contains at least one record from each checklist
- For each such set, call the common type specimen an "exemplar"
- (Again, we don't know much about the exemplar, but that doesn't matter
  since algorithmically we'll only be concerned with the associated set.)
- Typically an exemplar's record set would be two records, one in each
  checklist, with the same accepted species name.  However each name could
  be synonym, a subspecies name, etc. and the set could have multiple
  records from either or both checklists.
- See `exemplar` in [guide.md](guide.md#exemplar).

RCC-5 inference:
- Input: exemplars CSV file
- Exemplars play a role similar to protonyms in the Pyle/Remsen formulation.
- For each taxon concept, we consider the set of exemplars for it and all 
  taxon concepts hierarchically inferior to it in its checklist - that is, 
  the set of all exemplars that fall under the taxon concept.
- We take exemplar sets to be approximations to the (unknown) concepts.
- We can do set operations on exemplar sets to estimate RCC-5
  relationships between concepts.
- Occasionally there are multiple taxon concepts 
  in a checklist with the same exemplar set, requiring
  use of other information (hierarchy, names, ranks, and so on)
  to determine its relationships.
- See `plugin` in [guide.md](guide.md#plugin).

Inference report columns:
1. Taxon id for a record of A.
1. Canonical name for the record (redundant with taxonid).
1. If the A record is for a species, a list of relationships (semicolon 
   separated) for species
   whose taxon concepts are inferred to intersect that A taxon concept.
   Each relationship is given as the RCC-5
   relationship of the A concept to the B concept, along
   with the taxon id and canonical name of the B record.
   '-' means there may be intersecting species but the list was not computed 
   because the A record was not for a species.
1. The relationship to the smallest concept in the B checklist
   that contains the A record's concept.  The relationship 
   is given as above: RCC-5 relationship, taxon id, canonical name.
   As long as the A
   and B names are accepted, the A concept is either the same (RCC-5 =)
   as the B concept or larger (RCC-5 >) than it.  For synonyms it may
   be hard to tell.

Probably of most interest for understanding the impact of changes in
taxonomy going from A to B are the rows with multiple species given in
the intersecting species column.  These are situations where an A
species concept corresponds to multiple B species concepts,
potentially creating mislabeled data or requiring re-curation to
replace each use of an A species name with the appropriate B species
name.

The canonical names in the output are present for human readability;
for a more compact report they might be omitted, and obtained from the
checklists as needed.

Report rows for synonyms matching synonyms or matching nothing are suppressed.

[Installation](guide#installation):
- Assumes GNU bash and python3
- Assumes gnparser shell command
- Assumes python regex module
- There is no build step

Use:
- See user guide (guide.md)[guide.md] and example (example.md)[example.md]

