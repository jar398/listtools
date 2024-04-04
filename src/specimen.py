
import util
import property as prop
import simple
import parse

from util import log, UnionFindable

from checklist import get_source, blurb, blorb, get_parts, monitor, \
  is_accepted, get_duplicate_from
from workspace import get_workspace, get_children, isinA
  

# Specimens per se

# Get a specimen id on demand (never returns falsish)

def get_specimen_id(uf):
  r = uf.payload()
  (sid, x, y) = r
  if sid == None:
    # Create specimen id (for set operations) on demand
    w = x or y
    assert w
    ws = get_workspace(w)
    sid = len(ws.specimen_ufs)
    r[0] = sid
    ws.specimen_ufs[sid] = uf
    #log("# Specimen %s: (%s) <-> (%s)" % (sid, blorb(x), blorb(y)))
  return sid

# Assert identity of specimens

def same_specimens(uf, vf):     # formerly: same_typification_ufs
  return uf.find() is vf.find()

def equate_specimens(uf, vf):
  uf = uf.find()
  vf = vf.find()
  if uf is vf: return

  (i1, u1, v1) = uf.payload()
  (i2, u2, v2) = vf.payload()

  ef = uf.absorb(vf)          # ef happens to be uf
  r = ef.payload()

  # Not sure this is necessary, but should at least be harmless
  if i2 == None: r[0] = i1
  elif i1 == None: r[0] = i2
  elif i1 != i2: r[0] = min(i1, i2)

  # Choose or improve a record for which this specimen is to be its type.
  if _better_type_taxon(u1, u2):
    r[1] = u1; r[2] = u2
  else:
    r[1] = u2; r[2] = u1
  return ef

# This helps find, given a specimen, the most tipward taxon in a
# checklist containing that specimen.

def _pick_type_taxon(x, y):
  # See whether y is an improvement over x (lower unimportance value)
  if x == None: return y
  if y == None: return x
  assert get_source(x) is get_source(y)
  if x is y: return x
  d1 = get_duplicate_from(x, None)
  d2 = get_duplicate_from(y, None)
  if d2 and not d1: return x
  if d1 and not d2: return y
  a1 = is_accepted(x)
  a2 = is_accepted(x2)
  if a1 and not a2: return x
  if a2 and not a1: return y
  if not a1:
    return x                   # Arbitrary
  assert not x is x2
  m = simple.mrca(x, x2)
  # Choose the most tipward taxon
  if m == x: return y
  if m == x2: return x
  # x and y are disjoint??  How can this happen?  Arbitrary I guess
  # Prefer the one that has children
  c1 = len(get_children(x, ()))
  c2 = len(get_children(x2, ()))
  if c1 == 0 and c2 > 0: return c2
  if c2 == 0 and c1 > 0: return c1
  log("# Which is better as a type? %s %s" % (blorb(x), blorb(y)))
  assert False
  return x

# Access to specimen records [sid, x, y]

def sid_to_specimen(source, sid):
  return source.specimen_ufs[sid]

def sid_to_record(source, sid, z): # in A checklist
  uf = sid_to_specimen(source, sid)
  return uf.payload()[1]

def sid_to_opposite_record(source, sid, z): # in B checklist
  uf = sid_to_specimen(source, sid)
  return uf.payload()[2]

def sid_to_epithet(source, sid):
  parts = get_parts(sid_to_record(source, sid))
  return parts.epithet or parts.genus 

# -----------------------------------------------------------------------------

# x = smallest taxon in A checklist containing this specimen
# y = smallest taxon in B checklist containing this specimen

# typification_uf = specimen that is a type of some taxon.
# We might just call these "types"

(maybe_get_typification, set_typification_uf) = \
  prop.get_set(prop.declare_property("typification_uf",
                                     inherit=False))

# Only workspace nodes have uf records

def get_typification(x):
  probe = maybe_get_typification(x, None)
  if probe: return probe
  AB = get_workspace(x)
  r = [None, x, None] if isinA(AB, x) else [None, None, x]
  uf = UnionFindable(r)
  set_typification_uf(x, uf)
  return uf

# Are x and y known to have the same typification (type specimen)?

def same_typifications(x, y):
  uf = maybe_get_typification(x, None)
  vf = maybe_get_typification(y, None)
  return uf and vf and same_specimens(uf, vf)

# x and y are in workspace but may or may not be from same checklist

def equate_typifications(x, y):     # opposite checklists. x might be species
  if x is not y:
    equate_specimens(get_typification(x), get_typification(y))
    if monitor(x) or monitor(y):
      log("# Unifying exemplar(%s) = exemplar(%s)" % (blorb(x), blorb(y)))
  return x

# Are x and y equatable?
# x and y are assumed to be in the same checklist.
# This implements transitivity of equality, I think?

def equatable_typifications(x, y):
  uf = maybe_get_typification(x, None)
  if uf:
    vf = maybe_get_typification(y, None)
    if vf:
      return uf.payload() is vf.payload()
  return True

# For debugging mainly?

def get_typifies(spec):
  (sid, x, y) = spec.payload()
  return y or x

def bridge(AB, x, y):
  x = AB.in_left(x)
  y = AB.in_right(y)
  equate_typifications(x, y)

# -----------------------------------------------------------------------------
# An exemplar is a specimen that's the typification of two or more taxa in distinct
# checklists.

# Returns typification record (sid, x, y) or None.
# z is in AB.

def is_exemplar(uf):
  if uf:
    (_, x, y) = uf.payload()
    return x and y
  return False

def get_exemplar(z):            # Returns a uf node
  uf = maybe_get_typification(z, None)
  if is_exemplar(uf): return uf
  return None
  
def get_exemplar_record(z):
  uf = get_exemplar(z)
  if is_exemplar(uf): return uf.payload()
  return None

def get_exemplar_id(uf):
  if is_exemplar(uf):
    return get_specimen_id(uf)
  else:
    return None
