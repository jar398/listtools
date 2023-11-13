DEPRECATED


# List tools "plugin" summary 2023-11-10

Jonathan Rees - work in progress

Start with two checklists, call them A and B.

'Exemplar' determination:
- Find sets of records from the two checklists such that
    . in each set, all records have the same associated type specimen
    . all records known to have that type specimen are in the set
    . each set contains at least one record from each checklist
- For each such set, call the common type specimen an "exemplar"
  with respect to A and B.
- Typically an exemplar's record set would be two records, one in each
  checklist, with the same accepted species name.  However each name could
  be synonym, a subspecies name, etc. and the set could have multiple
  records from either or both checklists.
- See `exemplar` in [guide.md](guide.md#exemplar).

RCC-5 inference:
- Input: exemplars CSV file
- See `plugin` in [guide.md](guide.md#plugin).

See `plugin` in [guide.md](guide.md#plugin) for how to interpret the report.


[Installation](guide#installation):
- Assumes GNU bash and python3
- Assumes gnparser shell command
- Assumes python regex module
- There is no build step

Use:
- See user guide (guide.md)[guide.md] and example (example.md)[example.md]

