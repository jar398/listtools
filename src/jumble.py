from util import log
from checklist import *
from workspace import *
from rcc5 import *
from checklist import link_superior

import theory
from estimate import get_estimate, get_equivalent
from estimate import get_block  # for debugging
from ranks import ranks_dict

# link_superior(inf, sup)
#   inf is a Record, sup is a Relation
#   sup becomes the superior of inf
#   inf becomes a child or synonym of sup, depending 
#     on inf's taxonomicStatus


# Set superior and inferior links

def jumble_workspace(AB):
  log("# Jumbling")
  A_sups = B_sups = 0
  for x in preorder_records(AB.A): # includes AB.A.top
    z = AB.in_left(x)
    f = jumbled_superior(AB, z)
    if f:                       # should be accepted
      link_superior(z, f)
      A_sups += 1

  for y in preorder_records(AB.B):
    z = AB.in_right(y)
    f = jumbled_superior(AB, z)
    if f:
      link_superior(z, f)
      B_sups += 1
  log("# %s A superiors, %s B superiors" % (A_sups, B_sups))
  A_top = AB.in_left(AB.A.top)
  B_top = AB.in_right(AB.B.top)
  assert get_block(A_top) == get_block(B_top)
  A_roots = tuple(get_inferiors(A_top))
  B_roots = tuple(get_inferiors(B_top))
  assert A_roots or B_roots
  if A_roots and B_roots:
    log("# Residual roots from both A and B NYI: %s %s" %
        (map(blurb, A_roots), map(blurb, B_roots)))
    assert False

  if A_roots:                   # they're in AB, from either A or B
    log("# AB top is A top; roots = %s" %
        (tuple(map(blurb, A_roots)),))
    AB.top = A_top
  else:
    log("# AB top is B top; roots = %s" %
        (tuple(map(blurb, B_roots)),))
    AB.top = B_top

def set_workspace_top(AB):
  AB.top = get_workspace_top(AB)

# I recommend drawing a picture
# u could be in A or B in AB

def jumbled_superior(AB, u):
  # Suppress redundant nodes
  if is_redundant(AB, u):
    return None

  sup = local_sup(AB, u)      # Relation
  cos = cosuperior(AB, u)
  if not sup and not cos: return None
  elif not sup: prefer = cos        # u is A or B top in AB
  elif not cos: prefer = sup
  else:
    assert separated(sup.record, cos.record)
    assert local_accepted(AB, sup.record)
    assert local_accepted(AB, cos.record)

    rel = theory.compare(AB, sup.record, cos.record)
    # might be NOINFO
    if rel.relationship == GT or rel.relationship == GE:
      prefer = cos

    elif rel.relationship == LT or rel.relationship == LE:
      prefer = sup
    elif rel.relationship == EQ:
      prefer = sup                # sup is in A, cos is in B
    elif rel.relationship == DISJOINT:
      assert "should not happen"

    # OVERLAP NOINFO COMPARABLE INTERSECT ...
    # Superior should be MRCA ??

    else:
      # Way too many of these
      if local_accepted(AB, u):
        log("# Superiors of %s are:\n  %s %s %s" %
            (blurb(u), blurb(sup), rcc5_symbol(rel.relationship), blurb(cos)))
      prefer = sup
  assert local_accepted(AB, prefer.record)
  eq = is_redundant(AB, prefer.record)
  if eq:
    prefer = compose_relations(prefer, eq)
  return prefer

# Records to suppress

def is_redundant(AB, u):
  if isinA(AB, u): return None
  return get_equivalent(AB, u)
  #return (e_rel and
  #        get_rank(e_rel.record, None) == get_rank(u, None))

# Node's parent in opposite taxonomy.

def cosuperior(AB, u):
  est1 = get_estimate(u, None)  # v in other checklist
  if est1:                      # u is not top
    v = est1.record
    # u <= v
    # Cannot tolerate a synonym as parent
    if not is_accepted_locally(AB, v):
      # u <= v <= accepted-v
      est1 = compose_relations(est1, local_sup(AB, v))
      v = est1.record
    est2 = get_estimate(v, None)
    # u <= v <= u2
    if est2:      # u2 is not top
      u2 = est2.record              # in 1
      if u2 is u:                   # Congruent
        # u = v = u2
        # Try v's parent v3 (u = v < v3).
        sup3 = local_sup(AB, v)   # v <= v3 in 2
        if sup3:
          # u = v < v3
          answer = compose_relations(est1, sup3) # in 2
        else:
          return None           # u is at top
      else:
        # u < v <= u2 (cannot have u = v < u2)
        answer = est1           # v in checklist 2, noncongruent
    else:                       # u2 is top
      answer = est1
  else:
    return None                 # u is top
  # u < v
  # What if B.Balaenoptera physalus velifera <= A.Balaenoptera velifera*
  if not is_accepted_locally(AB, answer.record):
    # Shouldn't happen
    log("! nested synonyms: %s <= %s" % (blurb(u), blurb(answer.record)))
    assert False
  return answer

if __name__ == '__main__':
  import newick
  A = rows_to_checklist(parse_newick('a(b)'), {})
  B = rows_to_checklist(parse_newick('a(c)'), {})
  AB = make_workspace(A, B)
  jumble_workspace(AB)
  # checklist -> newick ?
