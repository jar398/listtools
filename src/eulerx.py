# Returns generator of lines (strings)

import sys
import theory

from rcc5 import *
from checklist import *

def generate_eulerx(AB, al):
  yield from generate_eulerx_checklist(AB.A)
  yield from generate_eulerx_checklist(AB.B)
  yield from eulerx_alignment(AB, al)

"""
taxonomy 2015 Prum
(Aves Gall_Neoa_Clade Palaeognathae)
(Gall_Neoa_Clade Galloanserae Neoaves)

taxonomy 2014 Jarvis
(Aves Neognathae Paleognathae)
(Paleognathae Struthioniformes Tinamiformes)

articulation 2015-2014 Prum-Jarvis
[2015.Aves is_included_in 2014.Aves]
[2015.Gall_Neoa_Clade equals 2014.Neognathae]
[2015.Palaeognathae is_included_in 2014.Paleognathae]
[2015.Galloanserae equals 2014.Struthioniformes]
[2015.Neoaves equals 2014.Tinamiformes]
"""

# Unique name of the sort Euler/X likes
# TBD: ensure name is unique (deal with multiple taxa with same canonical)

def get_eulerx_name(x, C=None):
  if not C: C = get_source(x)
  which = C.eulerx_which.get(get_primary_key(x))
  if not which: return None
  (e, i, cell) = which
  if cell[0] <= 1:
    return e
  else:
    if False:                   # lots of these in ncbi mammals
      tag = get_source_tag(x)
      print("# Discriminating: %s %s/%s, %s in %s" % (e, i, cell[0], get_primary_key(x), tag),
            file=sys.stderr)
    return "%s_%s" % (e, i)

def assign_eulerx_names(C):
  spin = {}    # maps eulerx base name to integer
  eulerx_which = {}
  C.eulerx_which = eulerx_which
  def traverse(x):
    # x is accepted
    e = get_eulerx_base_name(x)
    cell = spin.get(e, 0)
    if not cell:
      cell = [0]
      spin[e] = cell
    cell[0] += 1
    key = get_primary_key(x)
    eulerx_which[key] = (e, cell[0], cell)
    if get_rank(x, None) != 'species':
      for c in get_children(x, ()):
        traverse(c)
  traverse(C.top)

def get_eulerx_base_name(x):
  e = blurb(x)
  e = e.replace(' ', '_')
  e = e.replace('.', '_')
  # Nico's cosmetic preference
  e = e.replace('__', '_')
  return e

# x is a record (list)

def get_eulerx_qualified_name(x):
  src = get_source_tag(x)       # e.g. "A"
  return "%s.%s" % (src, get_eulerx_name(x))

def generate_eulerx_checklist(C):
  assign_eulerx_names(C)
  tag = get_tag(C)
  descr = checklist_description(C)
  yield ("taxonomy %s %s" % (tag, descr))
  names = []
  for rec in preorder_records(C):
    if rec != C.top:
      children = get_children(rec, None) # not the synonyms
      if children:
        sup_name = get_eulerx_name(rec, C)
        if sup_name and get_rank(rec, None) != 'species':
          names = list(map(lambda x:get_eulerx_name(x, C), children))
          names.sort()
          yield ("(%s %s)" % (sup_name, ' '.join(names)))
  yield ''

def checklist_description(C):
  return get_tag(C) + "_checklist"

def eulerx_alignment(AB, al):
  A = AB.A; B = AB.B
  yield ("articulation %s-%s %s-%s" %
         (get_tag(A), get_tag(B),
          checklist_description(A), checklist_description(B)))
  for (v, rel) in al:
    rel = theory.maybe_graft(v, rel) # Mark grafts with ! per Nico's suggestion
    w = rel.record
    x = get_outject(v); y = get_outject(w)
    if not is_top(x) and not is_top(y):
      yield eulerx_articulation(x, rel.relationship, y, rel.note)

def eulerx_articulation(x, ship, y, note):
  sym = rcc5_eulerx(ship)
  if ok_for_eulerx(ship):
    return "[%s %s %s]" % (get_eulerx_qualified_name(x),
                           sym,
                           get_eulerx_qualified_name(y))
  else:
    return ("#[%s %s %s] %s" %
            (get_eulerx_qualified_name(x),
             sym,
             get_eulerx_qualified_name(y),
             note or ''))

def ok_for_eulerx(ship):
  if False:
    if ship == EQ: return '='
    elif ship == LT: return '<'
    elif ship == GT: return '>'
    elif ship == OVERLAP: return '><'
    elif ship == DISJOINT: return '!'
    else: return False
  return True

