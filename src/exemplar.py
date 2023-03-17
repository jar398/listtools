#!/usr/bin/env python3

import sys, argparse, regex
import property as prop
import util, workspace, linkage
import rows

from util import log, MISSING
from property import mep_get, mep_set
from checklist import * 
from workspace import * 

# Exemplars (representing individual specimens or occurrences with
# known classification in BOTH checklists) are chosen heuristically
# via name matches.  Not every name match qualifies; only those that
# are 'tipward' do.  'Tipward' means that at least one of the two taxa
# has no tipward accepted descendant taxon.

def generate_exemplars(AB):
  yield ("exemplar", "A_rep", "B_rep") # representatives
  count = 0
  for (id, u, v) in choose_exemplars(AB):
    count += 1
    yield (id,
           get_primary_key(get_outject(u)) if u else MISSING,
           get_primary_key(get_outject(v)) if v else MISSING,
           blurb(u) if u else blurb(v))
  log("# exemplar report rows: %s" % count)

# Returns list of (id, v, w) ...
# Called from theory.theorize.
# Each exemplar has a representative in A and one in B.

def choose_exemplars(AB):
  analyze_tipwards(AB)

  generateme = []
  def emit(id, u, v):
    if (isinA(AB, u) if u else isinB(AB, v)):
      generateme.append((id, u, v))
    else:
      generateme.append((id, v, u))

  id_counter = [1]
  def do_exemplars(AB, inB):
    def traverse(x, species, demand):      # x in A
      # log("# considering %s" % (blurb(x)))
      u = AB.in_left(x)
      if get_rank(x, None) == 'species':
        species = u
        demand = [None]
      def get_species_xid():
        if demand[0] == None:
          demand[0] = id_counter[0]; id_counter[0] += 1
        return demand[0]

      hom = het = None
      for c in get_inferiors(x): # Generator
        # One child of c, together with its homotypic synonyms,
        # will have the same exemplar
        (more_hom, more_het) = traverse(c, species, demand)
        hom = hom or more_hom
        het = het or more_het

      # Invent a record-less exemplar for species that need them
      v = get_tipward(u, None)
      if v:
        if not get_exemplar(v, None):
          # Not seen in other checklist
          if hom:
            # Need a species for all the homotypics
            id = get_species_xid()
          else:
            id = id_counter[0]; id_counter[0] += 1
          set_exemplar(v, id)
          set_exemplar(u, id)
          emit(id, u, v)
        hom = (species and homotypic_synonyms(u, species))
        het = not hom

      # Species need special dispensation ... ?
      # Really, do this for every good match perhaps other than synonyms?
      elif (get_rank(u, None) == 'species' and
            het and not hom):
        set_exemplar(u, get_species_xid())
        emit(id, u, None)

      return (hom, het)

    traverse(AB.A.top, None, None)
  do_exemplars(AB, False)
  do_exemplars(swap(AB), True)

  return generateme

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
      u = AB.in_left(x)            # u is tipward...
      v = get_link(u, None)
      if monitor(v): log("# exemplars: tipwards %s %s" % (blurb(u), blurb(v)))
      if v:
        assert separated(u, v)
        # u is tipward and has a match, not necessarily tipward ...
        u2 = get_link(v)
        if not u2 or u2 != u:
          pass
          log("# Nonreciprocal tipwards %s %s" % (blurb(u), blurb(v)))
          # diagnose_match(u)
        else:
          set_tipward(u, v)
          set_tipward(v, v)
        seen = 1
        # log("# Tipwards %s %s" % (blurb(u), blurb(v)))

    return seen + seen2
  traverse(AB, AB.A.top)
  traverse(swap(AB), AB.B.top)

# exemplars

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
  args=parser.parse_args()

  a_name = 'A'; b_name = 'B'
  a_path = args.A
  b_path = args.B
  with rows.open(a_path) as a_rows:
    with rows.open(b_path) as b_rows:
      # compute name matches afresh
      AB = ingest_workspace(a_rows.rows(), b_rows.rows(),
                            A_name=a_name, B_name=b_name)
      linkage.find_links(AB)
      util.write_rows(util.write_rows(x_gen, sys.stdout),
                      sys.stdout)
