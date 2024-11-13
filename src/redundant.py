from util import log

from checklist import get_redundant, set_redundant, blurb, \
  all_records, get_parts, get_rank, is_accepted, \
  get_primary_key, get_tag, get_children
from proximity import get_family

def check_for_redundant(checklist):
  tag = get_tag(checklist)
  log("# Finding redundancies in %s checklist" % tag)
  index_by_name = {}    # (epithet, token, year)
  count = 0
  for x in all_records(checklist):
    p = get_parts(x)
    fam = get_family(x)
    # This doesn't seem so good.  Very limited
    key = (fam, p.genus, p.epithet, p.middle, p.token,
           p.year, get_rank(x, None))
    if key in index_by_name:
      index_by_name[key].append(x)
    else:
      index_by_name[key] = [x]
  for (key, xs) in index_by_name.items():
    if len(xs) > 1:
      xs.sort(key=redundancy_sort_key)
      best = xs[0]
      accept = []
      synos = []
      for x in xs:
        if is_accepted(x):
          accept.append(x)
        else:
          synos.append(x)
      # TBD: eliminate all synonym dups of accepted records, and report.
      # TBD: eliminate all of each all-synonym dup set, and report.
      # Only issue: when 2+ accepteds are dups of one another.
      for x in xs[1:]:
        if len(get_children(x, ())) == 0:
          count += 1
          set_redundant(x, best)
        else:
          log("# Redundant node has children: %s" % blurb(x))
      if len(accept) > 1:
        mess = "Redundant accepted(s), discarding some"
      elif len(synos) > 0:
        mess = "Redundant synonym(s), discarding"
      log("# %s: %s (%s) %s | %s" %
          (mess,
           blurb(best),
           key[0],
           ", ".join(map(lambda x:get_primary_key(x), accept)),
           ", ".join(map(lambda x:get_primary_key(x), synos))))
  if count > 0:
    log("# Redundant in %s: %s" % (tag, count))

def redundancy_sort_key(x):
  return (not is_accepted(x),
          -len(get_children(x, ())),
          get_primary_key(x)) 
