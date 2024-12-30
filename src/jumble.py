from util import log
from checklist import *
from workspace import *
from rcc5 import *

import theory
from estimate import get_estimate

def jumble_workspace(AB):
  for x in postorder_records(AB.A):
    z = AB.in_left(x)
    f = jumbled_superior(AB, z)
    if f: set_superior(z, f)
  for y in postorder_records(AB.B):
    z = AB.in_right(y)
    f = jumbled_superior(AB, z)
    if f: set_superior(z, f)
  AB.top = AB.in_left(AB.A.top)

# I recommend drawing a picture

def jumbled_superior(AB, u):
  sup = local_sup(AB, u)      # Relation
  cos = cosuperior(AB, u)
  if not sup: return cos
  if not cos: return sup
  assert separated(sup.record, cos.record)

  rel = theory.compare(AB, sup.record, cos.record)
  if rel.relationship == LT or rel.relationship == LE:
    return sup
  elif rel.relationship == GT or rel.relationship == GE:
    return cos
  elif rel.relationship == EQ:
    if isinA(AB, u):
      return sup                # sup is in A, cos is in B
    else:   # Do not add same concept twice
      return None               # cos is in A, sup is in B
  elif rel.relationship == OVERLAP:
    return sup
  else:  # rel.relationship == OVERLAP:
    log("# jumble: %s %s %s" % (blurb(sup), rcc5_symbol(rel.relationship), blurb(cos)))
    return sup

# Let's say u is in checklist, cos should be in checklist 2

def cosuperior(AB, u):
  est1 = get_estimate(u, None)  # v in 2
  if not est1: return None      # at top??
  v = est1.record
  est2 = get_estimate(v, None)  # u2 >= u in 1
  if not est2: return est1      # ????
  u2 = est2.record              # in 1
  if u2 is u:
    # v did not provide the 2nd superior.
    # Try v's parent (u -> v -> v3).
    sup3 = local_sup(AB, v)   # v3 in 2
    if not sup3: return None
    return compose_relations(est1, sup3) # in 2
  else:
    return est1                 # in 2

if __name__ == '__main__':
  import newick
  A = rows_to_checklist(parse_newick('a(b)'), {})
  B = rows_to_checklist(parse_newick('a(c)'), {})
  AB = make_workspace(A, B)
  jumble_workspace(AB)
  # checklist -> newick ?
