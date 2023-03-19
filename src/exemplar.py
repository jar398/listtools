#!/usr/bin/env python3

import sys, argparse, regex
import property as prop
import util, workspace, linkage, simple
import rows

from util import log, MISSING, UnionFindable
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

# See AB.exemplars for map {id: [id, v, w], ...}
# Called from theory.theorize.
# Each exemplar has a representative in A and one in B.

def choose_exemplars(AB):
  analyze_tipwards(AB)
  AB.exemplar_counter = 0
  AB.exemplars = {}

  def do_exemplars(CD, species):
    def traverse(x, species_uf):      # x in A
      # log("# considering %s" % (blurb(x)))
      u = CD.in_left(x)
      
      if get_rank(u, None) == 'species':
        species = u

      if species:
        if homotypic(u, species):
          equate_exemplars(species, u)

      v = get_tipward(u, None)
      if v:
        equate_exemplars(u, v)

      for c in get_inferiors(x):
        traverse(c, species)

  do_exemplars(AB, None)
  do_exemplars(swap(AB), None)
  equate_exemplars(AB.in_left(AB.A.top),
                   AB.in_right(AB.B.top))

# Issues: 
#   1. choice of taxa (both sides)
#   2. choice of name

def equate_exemplars(u, v):     # opposite checklists. u might be species
  uf = get_exemplar_info(u)
  vf = get_exemplar_info(v)
  merge_representatives(uf.find(), vf.find())
  uf.absorb(vf)

def merge_representatives(uf, vf): # absorb into uf
  exem = uf.id()
  (_, u1, v1) = exem
  (_, u2, v2) = vf.id()
  exem[1] = pick_preferred(u1, u2)
  exem[2] = pick_preferred(v1, v2)

def pick_preferred(u1, u2):
  if not u2: return u1
  if not u1: return u2
  rel = compare_per_checklist(u1, u2)
  if rel.relationship == GT:
    return u2                   # Prefer more tipward
  if (rel.relationship == EQ and
      get_parts(u1).moved and not get_parts(u2).moved):
    log("# preferring protonym %s %s" % (blurb(u1), blurb(u2)))
    return u2                   # Prefer protonym
  else:
    return u1

def get_exemplar_info(u):
  probe = really_get_exemplar_info(u, None)
  if probe: return probe
  AB = get_workspace(u)
  uf = UnionFindable([None, u, None] if isinA(AB, u) else [None, None, u])
  set_exemplar_info(u, uf)
  return uf

(really_get_exemplar_info, set_exemplar_info) = \
  prop.get_set(prop.declare_property("exemplar_info"))

# Returns exemplar id or None

def get_exemplar(u):
  uf = really_get_exemplar_info(u)
  if uf:
    (id, u, v) = uf.id()
    if u and v:
      if not id:
        # Create id on demand
        ws = get_workspace(u)
        ws.exemplar_counter += 1
        uf[0] = ws.exemplar_counter
        ws.exemplars[id] = uf
        log("# Exemplar %s %s %s" % (id, blurb(uf[1]), blurb(uf[2])))
        id = ws.exemplar_counter
      return uf
  return None

# -----------------------------------------------------------------------------
# Find tipward record matches (TRMs)

(get_tipward, set_tipward) = prop.get_set(prop.declare_property("tipward"))

# Postorder iteration over one of the two summands.
# Sets the 'tipward' property meaning we want to use the node to represent its
# exemplar.

def analyze_tipwards(AB):
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
        # u is tipward and has link, but link not necessarily tipward ...
        u2 = get_link(v)
        if u2:
          assert u2 == u
          log("# Nonreciprocal tipwards %s %s" % (blurb(u), blurb(v)))
          # diagnose_match(u)
          set_tipward(u, v)
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
      choose_exemplars(AB)
      writer = csv.writer(sys.stdout)
      writer.writerow(("exemplar_id", "A_taxonID", "B_taxonID"))
      for (id, exem) in AB.exemplars.items():
        (_, u, v) = exem
        writer.writerow((id, get_primary_key(u), get_primary_key(v)))
