from util import log

from checklist import get_redundant, set_redundant, blurb, \
  all_records, get_parts, get_rank, is_accepted, \
  get_primary_key, get_tag, get_children, \
  get_canonical, get_superior, get_source_tag, \
  monitor
from proximity import get_family

def check_for_redundant(checklist):
  if True:
    return
  index_by_name = {}    # (epithet, token, year)
  count = 0
  tag = get_tag(checklist)
  for x in all_records(checklist): # not incl top
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
          if monitor(x):
            log("# set_redundant(%s, %s)" % (blurb(x), blurb(best)))
          #set_redundant(x, best)
        else:
          log("# Redundant node has children: %s" % blurb(x))
      mess = None
      assert len(accept) + len(synos) >= 2
      if len(accept) > 1:
        mess = ("Redundant accepted(s) in %s, discarding some" %
                tag)
      else:    # 0 accept and 2+ syn, or 1 accept and 1+ syn
        mess = ("Redundant synonym(s) in %s, discarding" %
                tag)
      if mess:
        log("# %s: %s (%s) %s | %s" %
            (mess,
             blurb(best),
             key[0],
             ", ".join(map(lambda x:get_primary_key(x), accept)),
             ", ".join(map(lambda x:get_primary_key(x), synos))))
  if count > 0:
    log("# Redundant in %s: %s" % (tag, count))

# TBD: prefer child with name most similar to parent (usu. genus)
# get_canonical(x)
# get_superior(x)
# mismatch_point(x, y)

def redundancy_sort_key(x):
  sup = get_superior(x, None)   # relation
  supcan = get_canonical(sup.record, '') if sup else ''
  return (not is_accepted(x),
          -len(get_children(x, ())),
          mismatch_point(get_canonical(x, ''),
                         supcan) or 0,
          get_primary_key(x)) 

# -1 means no common suffix
# Smaller result means better match (earlier in sort sequence)
# Maybe use Levenstein distance instead??

def mismatch_point(s1, s2):
  return (-mismatch_point_from_start(s1, s2) +
          mismatch_point_from_end(s1, s2))

# 10 = good, 3 = bad

def mismatch_point_from_start(s1, s2):
  for i in range(0, min(len(s1), len(s2)), 1):
    if s1[i] != s2[i]:
      if i > 10:
        log("# %s compared to %s is %s" % (s1, s2, i))
      return i
  return 0

# -10 = good, -3 = bad

def mismatch_point_from_end(s1, s2):
  for i in range(-1, -min(len(s1), len(s2)), -1):
    if s1[i] != s2[i]:
      if i < -10:
        log("# %s compared to %s is %s" % (s1, s2, i))
      return i
  return -1
