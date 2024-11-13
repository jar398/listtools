from util import log

from checklist import get_redundant, set_redundant, blurb, \
  all_records, get_parts, get_rank, is_accepted, \
  get_primary_key
from proximity import get_family

def check_for_redundant(checklist):
  index_by_name = {}    # (epithet, token, year)
  for x in all_records(checklist):
    p = get_parts(x)
    fam = get_family(x)
    if True:
      key = (fam, p.genus, p.epithet, p.middle, p.token, p.year, get_rank(x, None))
    else:
      if p.epithet:
        key = (fam, p.epithet, p.middle, p.token, p.year)
      else:
        key = (fam, p.genus, p.middle, p.token, p.year)
    if key in index_by_name:
      index_by_name[key].append(x)
    else:
      index_by_name[key] = [x]
  for (key, xs) in index_by_name.items():
    if len(xs) > 1:
      accept = []
      reject = []
      for x in xs:
        if is_accepted(x):
          accept.append(x)
        else:
          reject.append(x)
      # TBD: eliminate all synonym dups of accepted records, and report.
      # TBD: eliminate all of each all-synonym dup set, and report.
      # Only issue: when 2+ accepteds are dups of one another.
      if len(reject) > 0:
        for x in reject:
          set_redundant(x, xs[0])
      if len(accept) > 1:
        mess = "Accepted duplicates, keeping one"
        accept.sort(key=lambda a: (-len(get_children(a, ())),
                                   get_primary_key(a)))
        for x in accept[1:]:
          set_redundant(x, accept[0])
      elif len(accept) == 1:
        mess = "Redundant synonym(s), discarding"
      else:
        mess = "Ambiguous synonyms, discarding"
      log("# %s: %s (%s) %s | %s" %
          (mess,
           blurb(xs[0]),
           key[0],
           ", ".join(map(lambda x:get_primary_key(x), accept)),
           ", ".join(map(lambda x:get_primary_key(x), reject))))


