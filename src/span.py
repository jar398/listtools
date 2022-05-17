
# Create a spanning tree for A+B based on B priority

# For each equivalent pair, only the B record will have a superior and inferiors.


import theory
from util import log

def span(AB):
  def traverse(AB, B_priority):
    for y in preorder_records(AB.B):
      if B_priority or not get_equivalent(AB, y):
        p_rel = possible_insertion(AB, y, B_priority)
        if p_rel:
          p = p_rel.record
          z = in_left(p)
          sup = relation(p_rel.relationship, z, p_rel.status, p_rel.note)
          plug_superior(y, sup, B_priority)
        else:
          q_rel = get_superior(y)
          if q_rel:
            z = in_right(q_rel.record)
            sup = relation(q_rel.relationship, z, q_rel.status, q_rel.note)
            plug_superior(y, sup, B_priority)
  traverse(AB, True)
  traverse(AB.flip(), False)
  ensure_inferiors_indexed(AB)

def plug_superior(y, sup, B_priority):
  if B_priority or not get_equivalent(y):
    set_superior(y, sup)
    if sup.relationship == ACCEPTED:
      ch = get_children(sup.record, None)
      if ch != None:
        ch.append(y)
      else:
        set_children(sup.record, [y])
    else:
      ch = get_synonyms(sup.record, None)
      if ch != None:
        ch.append(y)
      else:
        set_synonyms(sup.record, [y])

# y is in B (possibly flipped)

def possible_insertion(AB, y, B_priority):
  # Candidate superiors in A+B are {get_superior(y), q}
  p = cross_superior(AB, y)
  rel = get_superior(y)         # default
  q = rel.record
  ship = cross_relation(p, q).relationship
  if ship == LT:
    return rel             #q > p > y
  if ship == GT or ship == EQ:
    return None            #p > q > y  or p = q > y
  # CONFLICT or OVERLAP case.  Priority wins.
  if B_priority:
    log("# making arbitrary parent choice (B, priority)")
    return None    #p is the default
  else:
    log("# making arbitrary parent choice (A, non-priority)")
    return rel
