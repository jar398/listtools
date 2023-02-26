#!/usr/bin/env python3

import sys, argparse
import property as prop
import util, workspace

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

def exemplars(A_iter, B_iter, m_iter):
  A = rows_to_checklist(A_iter, {'tag': "A"})
  B = rows_to_checklist(B_iter, {'tag': "B"})
  AB = workspace.make_workspace(A, B, {'tag': "AB"})
  if not m_iter:
    m_iter = match_records(checklist_to_rows(A), checklist_to_rows(B))
  checklist.load_matches(m_iter, AB)
  yield ("exemplar", "A_taxonID", "B_taxonID")
  for (id, v, w) in choose_exemplars(AB):
    yield (id, get_primary_key(get_outject(v)), get_primary_key(get_outject(w)))

# Returns list of (id, v, w)

def choose_exemplars(AB):
  analyze_tipwards(AB)
  exemplar_records = []
  seen = prop.mep()             # set
  def do_exemplars(AB, inB):
    def traverse(x, emap):
      if emap == None and get_rank(x, None) == 'species':
        emap = {}
      for c in get_inferiors(x): # Generator
        traverse(c, emap)
      v = AB.in_left(x)
      rel = get_tipward(v, None)
      if rel:
        if emap == None: emap = {}
        w = rel.record          # in other checklist. not nec. tipward
        y = get_outject(w)
        e = epithetlike(x, y)
        ex = emap.get(e, None)
        if ex == None:
          id = len(exemplar_records)
          # Should save reason (in rel)??
          ex = (id, [], [])
          exemplar_records.append(ex)
          emap[e] = ex
        (_, vs, ws) = ex
        if not mep_get(seen, v, False):
          vs.append(v)
          mep_set(seen, v, True)
        if not mep_get(seen, w, False):
          ws.append(w)
          mep_set(seen, w, True)
    traverse(AB.A.top, None)
  do_exemplars(AB, False)
  do_exemplars(swap(AB), True)
  log("# %s exemplar records" % len(exemplar_records))
  return exemplar_records

def epithetlike(x, y):          # x and y match uniquely
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

# Stemmed last epithet.

def get_epithet(x):
  # Stemmed name will always be present if it's been run through gnparse,
  # absent otherwise.
  stemmed = get_stemmed(x, None)
  return stemmed.split(' ')[-1] if stemmed else None

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
      if monitor(x): log("# tipwards %s %s" % (blurb(x), blurb(rel)))
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

  a_path = args.A
  b_path = args.B
  m_path = args.matches
  with util.stdopen(a_path) as a_file:
    with util.stdopen(b_path) as b_file:
      def doit(m_iter):
        rows = exemplars(csv.reader(a_file),
                         csv.reader(b_file),
                         m_iter)
        writer = csv.writer(sys.stdout)
        for row in rows: writer.writerow(row)
      if m_path:
        with open(m_path) as m_file:
          doit(csv.reader(m_file))
      else:
        doit(None)
