#!/usr/bin/env python3

import sys, argparse, regex
import property as prop
import util, workspace
import rows

from util import log
from property import mep_get, mep_set
from checklist import * 
from workspace import * 
from match_records import match_records

# Exemplars (representing individual specimens or occurrences with
# known classification in BOTH checklists) are chosen heuristically
# via name matches.  Not every name match qualifies; only those that
# are 'tipward' do.  'Tipward' means that at least one of the two taxa
# has no tipward accepted descendant taxon.

def generate_exemplars(AB):
  yield ("exemplar", "A_taxonID", "B_taxonID")
  count = 0
  for (id, vs, rels, pnym) in choose_exemplars(AB):
    count += 1
    yield (id,
           get_primary_key(get_outject(v)),
           get_primary_key(get_outject(rel.record)),
           pnym)
  log("# exemplar report rows: %s" % count)

# Returns generator of (id, v, w, protonymic).
# Called from 

def choose_exemplars(AB):
  analyze_tipwards(AB)
  seen = prop.mep()             # set of Record
  counter = [0]
  def do_exemplars(AB, inB):
    def traverse(x, species_exem):      # x in A
      # log("# considering %s" % (blurb(x)))
      x_exem = get_proto(x)
      inf_exem = (x_exem
                  if get_rank(x, None) == 'species'
                  else species_exem)
      for c in get_inferiors(x): # Generator
        # One child of c, together with its homotypic synonyms, will have the same protonym
        yield from traverse(c, inf_exem)
      v = AB.in_left(x)
      rel = get_tipward(v, None)
      log("# tipward? %s %s" % (blurb(v), blurb(rel)))
      if rel:
        w = rel.record          # in other checklist. not nec. tipward
        y = get_outject(w)
        if x_exem:              # Passed down from ancestor species
          (_, _, _, x_proto) = x_exem
          proto = parse.unify_protos(x_proto, get_proto(y))
        else:
          proto = x_proto
        nymic = parse.unparse_proto(proto) if proto else MISSING
        if have == None:
          id = counter[0]; counter[0] += 1
          have = [id, None, None, nymic]
          log("# exemplar %s" % have)
        [_, v, w, _] = have
        if not mep_get(seen, v, False):
          have[1] = v
          mep_set(seen, v, have)
        if not mep_get(seen, w, False):
          have[2] = w
          mep_set(seen, w, have)
        if v and w:
          yield have
    yield from traverse(AB.A.top, None, None)
  yield from do_exemplars(AB, False)
  yield from do_exemplars(swap(AB), True)
  log("# exemplars: %s exemplar records" % counter[0])

def unify_protonymics(x, y):          # x and y match uniquely
  e1 = get_epithet(x)
  e2 = get_epithet(y)
  if e1 and e2:
    return min(e1, e2)
  if e1 or e2:
    return e1 or e2
  else:
    assert get_canonical(x) and get_canonical(y)
    return min(get_canonical(x), get_canonical(y)) # ??

def exemplar_id(ex): return ex[0]

# Selecting the 'best' name for an epithet.

def has_better_epithet(x, y):
  if y == None: return True
  if is_accepted(x) and not is_accepted(y): return True
  if rank(x) > rank(y): return True
  else: return False

def rank(x):
  r = get_rank(x)
  if r == 'species': return 5
  if r == 'subspecies': return 4
  if r == 'variety': return 3
  else: return 0

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
      rel = get_match(v, None)  # rel is a Relative...
      if monitor(x): log("# exemplars: tipwards %s %s" % (blurb(x), blurb(rel)))
      if rel and rel.relationship == EQ:
        assert separated(v, rel.record)
        # v is tipward and has a match, not necessarily tipward ...
        w = rel.record
        # v in get_matches(w) ... ??? no
        rel2 = get_matched(w)
        if not rel2 or rel2.record != v:
          pass
          log("# Nonreciprocal tipwards %s %s" % (blurb(w), blurb(rel2)))
          # diagnose_match(v)
        else:
          set_tipward(v, rel)
          set_tipward(w, rel2)
        seen = 1
        # log("# Tipwards %s %s" % (blurb(v), blurb(w)))
      else:
        pass
        # log("# No tipwards %s" % (blurb(v),))
    return seen + seen2
  traverse(AB, AB.A.top)
  traverse(swap(AB), AB.B.top)

# exemplar_records

def read_exemplars(inpath, AB):
  exemplars = {}
  with open(inpath) as infile:
    reader = csv.reader(infile)
    header = next(reader)
    id_col = windex(header, "exemplar")
    a_col = windex(header, "A_taxonID")
    b_col = windex(header, "B_taxonID")
    for row in reader:
      exemplars[id_col] = \
        (row[id_col],
         checklist.look_up_record(AB.A, row[a_col]),
         checklist.look_up_record(AB,B, row[b_col]))
  return exemplars

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Generate list of exemplars proposed for two checklists
    """)
  parser.add_argument('--A', help="the A checklist, as path name or -",
                      default='-')
  parser.add_argument('--B', help="the B checklist, as path name or -",
                      default='-')
  parser.add_argument('--matches', help="the A-B matches list")
  args=parser.parse_args()

  a_name = "A"; b_name = "B"
  a_path = args.A
  b_path = args.B
  m_path = args.matches
  with rows.open(a_path) as a_rows:
    with rows.open(b_path) as b_rows:
      if m_path:
        with rows.open(m_path) as m_rows:
          AB = ingest_workspace(a_rows.rows(), a_name, b_rows.rows(), b_name, m_rows.rows())
      else:
        # compute name matches afresh
        AB = ingest_workspace(a_rows.rows(), a_name, b_rows.rows(), b_name, None)
      x_gen = generate_exemplars(AB)
      util.write_rows(x_gen, sys.stdout)
