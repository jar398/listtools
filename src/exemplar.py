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

# Returns generator of (id, v, w, parts).      !!!!!!!!!
# Called from theory.theorize.
# Each exemplar has a representative in A and one in B.

def choose_exemplars(AB):
  analyze_tipwards(AB)
  seen = prop.mep()             # set of Record
  counter = [0]
  def do_exemplars(AB, inB):
    def traverse(x, species):      # x in A
      # log("# considering %s" % (blurb(x)))
      if get_rank(x, None) == 'species':
        species = x
      for c in get_inferiors(x): # Generator
        # One child of c, together with its homotypic synonyms,
        # will have the same exemplar
        yield from traverse(c, species)
      u = AB.in_left(x)
      rel = get_tipward(u, None)
      log("# tipward? %s %s" % (blurb(u), blurb(rel)))
      if rel:
        v = rel.record          # in other checklist. not nec. tipward
        have = mep_get(seen, u)
        if have:              # Passed down from ancestor species
          assert have[1] == v
        else:
          id = counter[0]; counter[0] += 1
          have = [id, u, v]
          log("# exemplar %s" % (id, blurb(u), blurb(w)))
          mep_set(seen, v, have)
          yield have
      elif (species and
            are_homotypic(u, species)):
        # Homotypic synonyms all have same exemplar


    yield from traverse(AB.A.top, None)
  yield from do_exemplars(AB, False)
  yield from do_exemplars(swap(AB), True)
  log("# exemplars: %s exemplar records" % counter[0])


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
