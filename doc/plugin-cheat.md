List tools "plugin" summary 2023-11-09

Jonathan Rees - work in progress

Scenario and theory:
- Start with two checklists, call them A and B.
- If you're doing regression analysis, think of A as the 'baseline' 
  checklist and B as 'proposed successor' (this is not the only use case)
- Each row of each checklist has an associated "taxon concept", but we
  don't know exactly what it is (and it would be challenging to find 
  out).  That doesn't prevent us from reasoning about the taxon concept.
- Each row of each checklist has a name (scientific or canonical),
  and the name is assumed to have an associated type specimen under
  control of the rules of nomenclature (even if we don't have any
  details on the specimen - that's OK)
- Therefore each row has an associated type [specimen].
- We do know that the type associated with a row falls under
  taxon concept for that row.

Name matching method:
- Find sets of rows from the two checklists such that
  . in each set, all rows have the same associated type
  . all rows known to have that type are in the set
  . each set contains at least one row from each checklist
- For each such set call the associated type specimen an "exemplar"
- (again, we don't "know what the exemplar is", but that doesn't matter
  since algorithmically we'll only be concerned with the associated set)
- Typically an exemplar's record set would be two records, one in each
  checklist, with the same accepted species name, but each name could
  also be synonym, a subspecies name, etc. and the set could have multiple
  rows from either or both checklists
- Output: CSV file of A-or-B, taxonid, exemplarid
    - says that the row given by taxonid in checklist A-or-B belongs
      to the record set given by exemplarid
    - A-or-B is 0 for A, 1 for B
    - exemplarid is locally unique to this analysis (not global)

RCC-5 inference:
- Exemplars play a role similar to protonyms in the Pyle/Remsen formulation
- N.b. taxon concepts may be distinct (different circumscriptions) even 
  when the exemplar and rank are the same.
- No matter what the implied taxon concept is for any row in 
  either checklist, we know which exemplars the concept contains.
  So we take exemplar sets to be proxies for the (unknown) concepts.
- We can do set operations on exemplar sets to estimate RCC-5
  relationships between concepts.

Inference report:
- First column gives taxon id for a row of A.
- Second column, if the A row is for a species, lists the B taxon ids for
  all B concepts (also for species) that intersect that A concept
  (according to the exemplar sets).  For each B concept the RCC-5
  relationship of the A concept to the B concept is given.
- Third column gives the taxon id for the smallest concept associated
  with a B row that contains the A row's concept.  As long as the A
  and B names are accepted, the A concept is either the same (RCC-5 =)
  as the B concept or larger (RCC-5 >) than it.  (For synonyms it may
  be impossible to tell.)

Installation:
- Assumes GNU make, bash, python3
- Assumes gnparser
- There is no build step

Use:
- A hodgepodge currently... I use the Makefile for some operations
  (such as 'start', 'use-gnparse', and downloading GBIF backbone
  snapshots) and go direct to the shell for others (like 'plugin')
  ... just have not yet made the Makefile catch up to the code
- doc files are out of date
- It's possible to separate name matching from RCC-5 analysis, to be
  documented

Question:
- How to script access to checklistbank
