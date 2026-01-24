# Returns generator of lines (strings)

import sys
import theory

from property import *
from rcc5 import *
from checklist import *
from workspace import get_workspace, swap, is_accepted_locally

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

def generate_eulerx(AB, al):
  yield from generate_eulerx_checklist(AB)
  yield from generate_eulerx_checklist(swap(AB))
  yield from eulerx_alignment(AB, al)

def generate_eulerx_checklist(C):
  assert C.workspace
  assign_eulerx_names(C)
  tag = get_tag(C.A)
  descr = checklist_description(C.A)
  yield ("taxonomy %s %s" % (tag, descr))
  names = []
  for x in preorder_records(C.A):
    if x != C.A.top:
      u = C.in_left(x)
      children = get_children(x, None) # not the synonyms
      if children:
        sup_name = get_eulerx_name(get_outject(u))
        if sup_name and get_rank(x, None) != 'species':
          names = list(map(lambda c:get_eulerx_name(c),
                           children))
          names.sort()
          yield ("(%s %s)" % (sup_name, ' '.join(names)))
  yield ''

def checklist_description(C):
  return get_tag(C) + "_checklist"

def eulerx_alignment(AB, al):
  assert AB.workspace
  A = AB.A; B = AB.B
  yield ("articulation %s-%s %s-%s" %
         (get_tag(A), get_tag(B),
          checklist_description(A), checklist_description(B)))
  for (u, rel, kind) in al:
    # rel = theory.maybe_graft(u, rel) # Mark grafts with ! per Nico's suggestion
    w = rel.record
    if u and w:
      assert get_workspace(u)
      assert is_accepted_locally(AB, u)
      assert is_accepted_locally(AB, w), blurb(w)
      if not is_top(get_outject(u)) and not is_top(get_outject(w)):
        yield eulerx_articulation(u, rel.relation, w, rel.note)

def eulerx_articulation(u, ship, v, note):
  sym = rcc5_eulerx(ship)
  name1 = get_eulerx_qualified_name(u)
  name2 = get_eulerx_qualified_name(v)
  return "[%s %s %s]" % (name1, sym, name2)


# Unique name of the sort Euler/X likes
# TBD: ensure name is unique (deal with multiple taxa with same canonical)

def get_eulerx_name(x):
  A = get_source(x)
  which = mep_get(A.eulerx_which, x, None) # (ecanon, i, cell)
  if not is_accepted(x):
    log("# %s not accepted in articulation" % blurb(x))
    assert False
  if not which:
    log("# No euler/x name for %s" % blurb(x))
    assert False
  (ecanon, i, cell) = which
  if cell[0] <= 1:
    return ecanon
  else:
    if False:                   # lots of these in ncbi mammals
      tag = get_source_tag(x)
      print("# Discriminating: %s %s/%s, %s in %s" % (ecanon, i, cell[0], get_primary_key(x), tag),
            file=sys.stderr)
    return "%s_%s" % (ecanon, i)

def assign_eulerx_names(C):
  def doit(C):
    eulerx_which = mep()
    spin = {}    # maps eulerx base name to integer
    def traverse(x):
      if is_accepted(x):
        ecanon = get_eulerx_base_name(x)  #NOT in workspace
        cell = spin.get(ecanon, None)
        if cell == None:
          cell = [0]
          spin[ecanon] = cell
        cell[0] += 1
        mep_set(eulerx_which, x, (ecanon, cell[0], cell))
        if get_rank(x, None) != 'species':
          for c in get_inferiors(x):  # children only ??
            traverse(c)
    traverse(C.A.top)
    C.A.eulerx_which = eulerx_which
  doit(C)
  doit(swap(C))

def get_eulerx_base_name(x):
  e = get_ok_name(x)
  e = e.replace(' ', '_')
  e = e.replace('.', '_')
  # Nico's cosmetic preference
  e = e.replace('__', '_')
  if not is_accepted(x): e = e + "_SYN"
  return e

# x is a record (list)

def get_eulerx_qualified_name(z):
  x = get_outject(z)
  src = get_source_tag(x)       # e.g. "A"
  ename = get_eulerx_name(x)
  assert ename                  # why why why
  return "%s.%s" % (src, ename)

