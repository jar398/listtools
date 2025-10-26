from util import log
from checklist import *
from workspace import *
from rcc5 import *
from checklist import link_superior

import util, theory
from estimate import get_estimate, get_equivalent

from estimate import get_block  # for debugging

# Set superior and inferior links

def jumble_workspace(AB):
  log("# Jumbling")
  util.log_allowance = 1000
  set_workspace_top(AB)
  links = 0
  for x in preorder_records(AB.A): # includes AB.A.top
    z = AB.in_left(x)
    f = jumbled_superior(AB, z)
    # f could be < parent of z in A
    if f:                       # should be accepted
      link_superior(z, f)
      if is_top(f.record) or is_top(z):
        log("# Linking near top: %s %s" %
            (blurb(z), blurb(f)))
      links += 1
  log("# made %s parent links for A records in AB" % links)
  links = 0
  for y in preorder_records(AB.B):
    z = AB.in_right(y)
    f = jumbled_superior(AB, z)
    if f:
      link_superior(z, f)
      links += 1
  log("# made %s parent links for B records in AB" % links)
  log("# workspace top is %s" % blurb(AB.top))
  infs = list(get_inferiors(AB.top))
  log("# inferiors of top are %s" % list(map(blurb, infs)))
  assert len(infs) > 0

def set_workspace_top(AB):
  AB.top = get_workspace_top(AB)

def get_workspace_top(AB):
  A_top = AB.in_left(AB.A.top)
  B_top = AB.in_right(AB.B.top)
  assert get_block(A_top) == get_block(B_top)
  rel = theory.compare(AB, A_top, B_top)
  log("# Top comparison: %s %s %s" %
      (blurb(A_top), rcc5_symbol(rel.relationship), blurb(B_top)))

  if rel.relationship == EQ or rel.relationship == GT or rel.relationship == GE:
    answer = A_top
    if get_equivalent(AB, B_top):
      if rel.relationship == EQ: assert A_top == get_equivalent(AB, B_top).record
      # could be DISJOINT or OVERLAP
    else:
      log("# losing B top: %s %s %s" %
          (blurb(A_top), rcc5_symbol(rel.relationship), blurb(B_top)))
  else:  # rel.relationship == LT or rel.relationship == LE:
    answer = B_top
    if not get_equivalent(AB, A_top):
      log("# losing A top: %s %s %s" %
          (blurb(A_top), rcc5_symbol(rel.relationship), blurb(B_top)))
  return answer

def supremum(AB, u, v):
  return None

# I recommend drawing a picture

# 'u' could be in either A or B

def jumbled_superior(AB, u):
  # Suppress nodes in B that have an equivalent in A
  if isinB(AB, u) and get_equivalent(AB, u):
    # Omit from jumbled hierarchy; redundant
    return None

  if isinA(AB, u):
    sup = local_sup(AB, u)      # Relation from u to another record in AB
    cos = cosuperior(AB, u)
  else:
    sup = cosuperior(AB, u)
    cos = local_sup(AB, u)
  # We want min(sup, cos).
  if not sup:
    return cos    # u has no parent in AB.A
  if not cos:
    return sup    # u has no parent in AB.B
  # sup and cos are both in AB
  assert separated(sup.record, cos.record)
  assert local_accepted(AB, sup.record)
  assert local_accepted(AB, cos.record)

  # Find the supremum.  sup in AB is derived from A.
  # Trouble case: theory doesn't know whether it's = or >

  while True:

    rel = theory.compare(AB, sup.record, cos.record)
    if rel.relationship == LT or rel.relationship == LE \
       or rel.relationship == EQ:
      return sup

    elif rel.relationship == GT or rel.relationship == GE:
      return cos

    else:
      sup2 = get_superior(sup.record, None)
      if rel.relationship != OVERLAP or not sup2:
        # rel.relationship == anything else: NOINFO, etc.
        # might be NOINFO
        log("# jumble at %s:\n  %s %s %s" %
            (blurb(u), blurb(sup), rcc5_symbol(rel.relationship), blurb(cos)))
      if not sup2:
        return sup              # ?

      # Loop around with new sup
      sup = sup2

# Let's say u is in checklist 1.  Cosuperior should be in checklist 2.
# Answer is None for record in checklist 2 whose sup is 
# "lifted" i.e. congruent to a checklist 1 node.

def cosuperior(AB, u):
  est1 = get_estimate(u, None)
  if not est1: return None
  v = est1.record

  assert separated(u, v)

  est2 = get_estimate(v, None)
  if not est2: return None
  u2 = est2.record

  assert separated(u2, v)

  if u2 is u:                   # Equivalent
    sup = local_sup(AB, v)
    if not sup: return None
    # u -> v
    answer = compose_relations(est1, sup)
  else:
    # u < u2, or u = u2 and some minor difference
    answer = est1

  assert separated(u, answer.record)

  return answer


def cosuperior_old(AB, u):
  # Find v with u != v and u <= v
  est1 = get_estimate(u, None)  # v in checklist 2
  if est1:                      # u is not top
    # u <= v
    v = est1.record
    # if v is a synonym, we would have u <= v, which is no good
    if not is_accepted_locally(AB, v):
      # u <= v <= accepted-v
      acc = local_accepted(AB, v)    # same thing I guess
      est1 = compose_relations(est1, acc)
      v = est1.record
      # Now we have u <= v and it's OK
    # est1: u -> v
    est2 = get_estimate(v, None)
    # u <= v <= u2
    if est2:      # u2 is not top
      u2 = est2.record              # in 1
      # Should eliminate synonym u <= u2 as we did above by lifting u2
      if not is_accepted_locally(AB, u2):
        acc = local_accepted(AB, u2)
        est2 = compose_relations(est2, acc)
        u2 = est2.record

      if u2 is u:                   # Congruent; co-estimates
        # v did not provide the 2nd superior.
        # Try v's parent (u -> v -> v3).
        sup3 = local_sup(AB, v)   # v <= v3 in 2
        if sup3:
          # u <= v <= v3
          # We want overall est1 -> u -> v -> u2 -> est2
          #   where est1: u -> v
          #   WHERE EST2: V -> U2
          answer = compose_relations(est1, sup3) # in 2
        else:
          return None           # v is top
      else:
        answer = est1           # v in checklist 2, noncongruent
    else:                       # u2 is top
      answer = est1
  else:
    return None                 # u is top
  # e.g. B.Balaenoptera physalus velifera <= A.Balaenoptera velifera*
  if not is_accepted_locally(AB, answer.record):
    # Shouldn't happen
    log("! nested synonyms: %s <= %s" % (blurb(u), blurb(answer.record)))
    assert False
  # assert u.record < answer.record
  return answer

if __name__ == '__main__':
  import newick
  A = rows_to_checklist(parse_newick('a(b)'), {})
  B = rows_to_checklist(parse_newick('a(c)'), {})
  AB = make_workspace(A, B)
  jumble_workspace(AB)
  # checklist -> newick ?
