from util import log
from checklist import *
from workspace import *
from rcc5 import *
from checklist import link_superior

import theory
from estimate import get_estimate, get_equivalent

from estimate import get_block  # for debugging

# Set superior and inferior links

def jumble_workspace(AB):
  log("# Jumbling")
  set_workspace_top(AB)
  i = eq = 0
  for x in preorder_records(AB.A): # includes AB.A.top
    z = AB.in_left(x)
    f = jumbled_superior(AB, z)
    if f:                       # should be accepted
      link_superior(z, f)
      if is_top(f.record) or is_top(z):
        log("# Linking near top: %s %s" %
            (blurb(z), blurb(f)))
    else: eq += 1
    i += 1
  log("# A: %s jumbled parent/child links, %s congruences" % (i, eq))
  i = eq = 0
  for y in preorder_records(AB.B):
    z = AB.in_right(y)
    f = jumbled_superior(AB, z)
    if f:
      link_superior(z, f)
    else: eq += 1
    i += 1
  log("# B: %s jumbled parent/child links, %s suppressed congruents" % (i, eq))
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

# I recommend drawing a picture

def jumbled_superior(AB, u):
  sup = local_sup(AB, u)      # Relation
  cos = cosuperior(AB, u)
  if not sup: return cos        # u is AB.top
  if not cos: return sup
  assert separated(sup.record, cos.record)
  assert local_accepted(AB, sup.record)
  assert local_accepted(AB, cos.record)

  # What if the result is a synonym?

  rel = theory.compare(AB, sup.record, cos.record)
  # might be NOINFO
  if rel.relationship == LT or rel.relationship == LE:
    answer = sup
  elif rel.relationship == GT or rel.relationship == GE:
    answer = cos
  elif rel.relationship == EQ:
    if isinA(AB, u):
      answer = sup                # sup is in A, cos is in B
    else:   # Do not add same concept twice
      answer = None               # cos is in A, sup is in B
  elif rel.relationship == OVERLAP:
    answer = sup
  else:  # rel.relationship == anything else: NOINFO, OVERLAP, etc.
    log("# jumble at %s:\n  %s %s %s" %
        (blurb(u), blurb(sup), rcc5_symbol(rel.relationship), blurb(cos)))
    answer = sup
  if answer: assert local_accepted(AB, answer.record)
  return answer

# Let's say u is in checklist 1.  Cosuperior should be in checklist 2.
# Answer is None for record in checklist 2 whose sup is 
# congruent to a checklist 1 node.

def cosuperior(AB, u):
  est1 = get_estimate(u, None)  # v in checklist 2
  if est1:                      # u is not top
    # u <= v
    v = est1.record
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
        # v did not provide the 2nd superior.
        # Try v's parent (u -> v -> v3).
        sup3 = local_sup(AB, v)   # v <= v3 in 2
        if sup3:
          # u <= v <= v3
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
  return answer

if __name__ == '__main__':
  import newick
  A = rows_to_checklist(parse_newick('a(b)'), {})
  B = rows_to_checklist(parse_newick('a(c)'), {})
  AB = make_workspace(A, B)
  jumble_workspace(AB)
  # checklist -> newick ?
