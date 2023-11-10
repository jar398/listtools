# List tools "plugin" summary 2023-11-10

Jonathan Rees - work in progress

Scenario:
- Start with two checklists, call them A and B.
- If you're doing regression analysis, think of A as the 'baseline' 
  checklist and B as 'proposed successor' (this is not the only use case)
- Checklists must be in Darwin Core format.

Semantics of checklists:
- Each row of each checklist has an associated "taxon concept"
  determined by the fields of the record with the checklist as context.
  We may not know what the taxon concept is and it would be challenging to find 
  out.  However, that doesn't prevent us from reasoning about it
  using information in the checklist.
- Each row of each checklist has a name ('scientific' if in includes authority 
  information; 'canonical' if not),
  and the name is assumed to have an associated type specimen under
  control of the rules of biological nomenclature (again, we may not have any
  details at hand on the type specimen - that's OK).  Again, in case of 
  ambiguity, the checklist can be consulted for context.
- Therefore each row has an associated type specimen, via the name.
- We know that the type specimen associated with a row falls under
  the taxon concept associated with that row.

'Exemplar' determination:
- Find sets of rows from the two checklists such that
    . in each set, all rows have the same associated type specimen
    . all rows known to have that type specimen are in the set
    . each set contains at least one row from each checklist
- For each such set, call the common type specimen an "exemplar"
- (Again, we don't know much about the exemplar, but that doesn't matter
  since algorithmically we'll only be concerned with the associated set.)
- Typically an exemplar's record set would be two records, one in each
  checklist, with the same accepted species name.  However each name could
  be synonym, a subspecies name, etc. and the set could have multiple
  rows from either or both checklists.
- Output: exemplars CSV file.  Rows are: A-or-B, taxonid, exemplarid
    - A-or-B is 0 for A, 1 for B
    - exemplarid is locally unique to this analysis (not global)
    - A row says that the row given by taxonid in checklist A-or-B belongs
      to the record set associated with the exemplar given by exemplarid

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

Inference report columns:
1. Taxon id for a row of A.
1. Canonical name for the row (redundant with taxonid).
1. If the A row is for a species, a list of relationships (semicolon 
   separated) for species
   whose taxon concepts are inferred to intersect that A taxon concept.
   Each relationship is given as the RCC-5
   relationship of the A concept to the B concept is given, along
   with the taxon id and canonical name of the B row.
   '-' means there may be intersecting species but the list was not computed 
   because the A row was not for a species.
1. The relationship to the smallest concept associated
   with a B row that contains the A row's concept.  The relationship 
   is given as above: RCC-5 relationship, taxon id, canonical name.
   As long as the A
   and B names are accepted, the A concept is either the same (RCC-5 =)
   as the B concept or larger (RCC-5 >) than it.  For synonyms it may
   be hard to tell.

The canonical names in the output are present for human readability;
for a more compact report they might be omitted, and obtained from the
checklists as needed.

Report rows for synonyms matching synonyms or matching nothing are suppressed.

Installation:
- Assumes GNU bash and python3
- Assumes gnparser shell command
- Assumes python regex module
- There is no build step

Use:
- See user guide (guide.md)[guide.md]
