
from util import log
from checklist import * 
from workspace import * 

# Exemplars (representing individual specimens or occurrences with
# known classification in BOTH checklists) are chosen heuristically
# via name matches.  Not every name match qualifies; only those that
# are 'tipward' do.  'Tipward' means that at least one of the two taxa
# has no tipward accepted descendant taxon.

# TBD: Be able to write and read the exemplars.

def choose_exemplars(AB):
  analyze_tipwards(AB)
  exemplar_records = []
  def do_exemplars(AB, inB):
    def traverse(x):
      v = AB.in_left(x)
      for c in get_inferiors(v):
        traverse(c, inB)
      probe = get_exemplar(v, None) # (id, v, w)
      if probe: return
      rel = get_tipward(v, None)
      if rel:
        if inB:
          w = v
          v = rel.record
        else:
          w = rel.record
        assert separated(v, w)
        id = len(exemplar_records)
        ex = (id, v, w)
        exemplar_records.append(ex)
        set_exemplar(v, ex)
        set_exemplar(w, ex)
        return id                   # might be 0
    traverse(AB.A.top)
  do_exemplars(AB, False)
  do_exemplars(swap(AB), True)
  AB.exemplar_records = exemplar_records

def exemplar_id(ex): return ex[0]

(get_exemplar, set_exemplar) = prop.get_set(prop.declare_property("exemplar"))


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

