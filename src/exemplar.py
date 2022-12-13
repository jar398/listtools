#!/usr/bin/env python3

from util import log
from checklist import * 
from workspace import * 
import property as prop

# Exemplars (representing individual specimens or occurrences with
# known classification in BOTH checklists) are chosen heuristically
# via name matches.  Not every name match qualifies; only those that
# are 'tipward' do.  'Tipward' means that at least one of the two taxa
# has no tipward accepted descendant taxon.

# TBD: Be able to write and read the exemplars.

def choose_exemplars(AB):
  analyze_tipwards(AB)
  exemplar_records = []
  seen = prop.mep()
  j = [0]
  def do_exemplars(AB, inB):
    def traverse(x):
      for c in get_inferiors(x):
        traverse(c)
      v = AB.in_left(x)
      rel = get_tipward(v, None)
      j[0] += 1
      if rel:
        if inB:
          w = v
          v = rel.record
        else:
          w = rel.record
        if prop.mep_get(seen, v, False): return
        prop.mep_set(seen, v, True)
        prop.mep_set(seen, w, True)
        assert separated(v, w)
        id = len(exemplar_records)
        # Should save reason (in rel)??
        ex = (id, v, w)
        exemplar_records.append(ex)
    traverse(AB.A.top)
  do_exemplars(AB, False)
  do_exemplars(swap(AB), True)
  log("# visited %s records" % j[0])
  log("# %s exemplar records" % len(exemplar_records))
  return exemplar_records

def exemplar_id(ex): return ex[0]


# Find tipward record matches (TRMs)

(get_tipward, set_tipward) = prop.get_set(prop.declare_property("tipward"))

def analyze_tipwards(AB):
  find_tipwards(AB)
  find_tipwards(AB.swap())

# Postorder iteration over one of the two summands.
# Sets the 'tipward' property meaning we want to use the node to represent its
# exemplar.

def find_tipwards(AB):
  def traverse(AB, x):
    seen = seen2 = 0
    for c in get_children(x, ()):
      seen += traverse(AB, c)
    for c in get_synonyms(x, ()):
      seen2 += traverse(AB, c)
    if seen == 0:
      # not seen means that this node could be tipward
      v = AB.in_left(x)            # v is tipward...
      rel = get_matched(v)  # rel is a Relative...
      if monitor(x): log("tipwards %s %s" % (blurb(x), blurb(rel)))
      if rel:
        assert separated(v, rel.record)
        # v is tipward and has a match, not necessarily tipward ...
        set_tipward(v, rel)
        if True:
          w = rel.record
          rel2 = get_matched(w)
          assert rel2           # symmetric
          assert rel2.record == v
          set_tipward(w, rel2)
        seen = 1
    return seen + seen2
  traverse(AB, AB.A.top)
  traverse(swap(AB), AB.B.top)

