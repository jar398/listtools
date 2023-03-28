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

# See AB.exemplars for map {id: [id, v, w], ...}
# Called from theory.theorize.
# Each exemplar has a representative in A and one in B.

def analyze_exemplars(AB):
  analyze_tipwards(AB)

  def do_exemplars(CD):
    def traverse(x, species):      # x in A
      u = CD.in_left(x)
      
      if species:
        if linkage.homotypic(u, species):
          equate_exemplars(species, u)

      if get_rank(u, None) == 'species':
        species = u

      v = get_tipward(u, None)
      if v:                     # not None or False
        u_back = get_tipward(v, None)
        if u_back is u:         # Mutual?
          if monitor(u) or monitor(v):
            log("# exemplar: Mutual: %s -> %s" % (blurb(u), blurb(v)))
          equate_exemplars(u, v)
        elif u_back == None:
          if monitor(u) or monitor(v):
            log("# exemplar: Nonmutual: %s -> %s -> (none)" % 
                (blurb(u), blurb(v)))
          equate_exemplars(u, v)
        elif u_back == False:
          if monitor(u) or monitor(v):
            log("# exemplar: Ambiguous: %s -> %s -> (2 or more)" % 
                (blurb(u), blurb(v)))
          equate_exemplars(u, v)
        else: # u_back != None:     doesn't seem to happen?
          log("# exemplar: Returns to different tipward %s -> %s -> %s" % 
              (blurb(u), blurb(v), blurb(u_back)))

      for c in get_inferiors(x):
        traverse(c, species)
    traverse(CD.A.top, None)

  do_exemplars(AB)
  do_exemplars(swap(AB))
  equate_exemplars(AB.in_left(AB.A.top),
                   AB.in_right(AB.B.top))
  report_on_exemplars(AB)

def report_on_exemplars(AB):
  exemplars = set()
  for x in preorder_records(AB.A):
    u = AB.in_left(x)
    ex = get_exemplar(u)
    if ex:
      exemplars.add(ex[0])
  log("-- Exemplars: %s " % len(exemplars))

# Issues: 
#   1. choice of taxa (both sides)
#   2. choice of name

def equate_exemplars(u, v):     # opposite checklists. u might be species
  if u != v:
    uf = get_exemplar_uf(u)
    vf = get_exemplar_uf(v)
    merge_representatives(uf.find(), vf.find())
    uf.absorb(vf)

def merge_representatives(uf, vf): # absorb into uf
  r = uf.payload()
  (i1, u1, v1) = r
  (i2, u2, v2) = vf.payload()
  r[1] = unify_records(u1, u2)
  r[2] = unify_records(v1, v2)
  if r[1] and r[2]:
    #log("# Same exemplar: e(%s) = e(%s)" % (blurb(r[1]), blurb(r[2])))
    pass
  if i1 and i2:
    r[0] = min(i1, i2)
  else:
    r[0] = i1 or i2

def unify_records(u1, u2):
  if not u2: return u1
  if not u1: return u2
  rel = simple.compare_per_checklist(u1, u2)
  if rel.relationship == GT:
    return u2                   # Prefer more tipward
  if (rel.relationship == EQ and
      not get_parts(u1).protonymp and get_parts(u2).protonymp):
    log("# preferring protonym %s to %s" % (blurb(u2), blurb(u1)))
    return u2                   # Prefer protonym
  else:
    return u1

def get_exemplar_uf(u):
  #log("# Thinking about exemplar for %s" % blurb(u))
  probe = really_get_exemplar_uf(u, None)
  if probe: return probe
  AB = get_workspace(u)
  uf = UnionFindable([None, u, None] if isinA(AB, u) else [None, None, u])
  set_exemplar_uf(u, uf)
  return uf

(really_get_exemplar_uf, set_exemplar_uf) = \
  prop.get_set(prop.declare_property("exemplar_uf"))


# Apply this to an exemplar id to obtain an exemplar union/find node,
# and return the associated taxon record that's in same checklist as z.

def xid_to_record(AB, id, z):
  uf = AB.exemplars[id]
  (_, u, v) = uf.payload()
  return u if isinA(AB, z) else v

def xid_to_opposite_record(AB, id, z):
  uf = AB.exemplars[id]
  (_, u, v) = uf.payload()
  return v if isinA(AB, z) else u


# Returns exemplar record or None.  TO BE USED ONLY AFTER
# analyze_exemplars HAS FINISHED ITS WORK.

def get_exemplar(u):
  uf = really_get_exemplar_uf(u, None)
  if uf:
    r = uf.payload()
    (id, u, v) = r
    if u and v:
      if not id:
        # Create id (for set operations) on demand
        ws = get_workspace(u)
        id = fresh_exemplar_id(ws)
        r[0] = id
        ws.exemplars[id] = uf
        #log("# Exemplar %s: e(%s) = (%s)" % (id, blurb(u), blurb(v)))
      return r
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
      v = get_link(u, None)        # same name -> same exemplar(s?)
      if monitor(v):
        log("# exemplar: Tipward link %s -> %s" % (blurb(u), blurb(v)))
      if v != None:                     # not None, not False = ambiguous
        if v != False:
          assert separated(u, v)
          seen = 1
        # u is tipward and has link, but link not necessarily tipward ...
        set_tipward(u, v)
        if monitor(u):
          log("# exemplar: Set tipward %s -> %s" % (blurb(u), blurb(v)))
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

# could do this as a generator + write_rows

def write_exemplar_list(AB, out=sys.stdout):
  util.write_rows(generate_exemplars(AB), out)

def generate_exemplars(AB):
  yield ("exemplar", "A_taxonID", "B_taxonID", "A_blurb", "B_blurb") # representatives
  count = 0
  seen = {}
  for x in preorder_records(AB.A):
    z = AB.in_left(x)
    r = get_exemplar(z)
    if r:
      (id, u, v) = r
      if u == z:
        u_key = get_primary_key(get_outject(u))
        v_key = get_primary_key(get_outject(v))
        yield (id, u_key, v_key, blurb(u), blurb(v))
        count += 1
  log("# %s rows in exemplar report" % count)

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
      analyze_exemplars(AB)
      write_exemplar_list(AB)
