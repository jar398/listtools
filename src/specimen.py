
import util
import property as prop
import simple
import parse

from util import log, UnionFindable

from checklist import get_source, blurb, blorb, get_parts, monitor, \
  is_accepted, get_duplicate_from
from workspace import get_workspace, get_children, \
  get_outject, isinA
  

# Specimens per se

# Get a specimen id on demand (never returns falsish)

def get_specimen_id(uf):
  r = uf.payload()
  (sid, u, v) = r
  if sid == None:
    # Create specimen id (for set operations) on demand
    w = u or v
    assert w
    ws = get_workspace(w)
    sid = len(ws.specimen_ufs)
    r[0] = sid
    ws.specimen_ufs[sid] = uf
    #log("# Specimen %s: (%s) <-> (%s)" % (sid, blorb(u), blorb(v)))
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
  r[1] = _pick_type_taxon(u1, u2)
  r[2] = _pick_type_taxon(v1, v2)
  return ef

# This helps find, given a specimen, the most tipward taxon in a
# checklist containing that specimen.

def _pick_type_taxon(v1, v2):
  # See whether v2 is an improvement over v1 (lower unimportance value)
  if v1 == None: return v2
  if v2 == None: return v1
  assert get_source(v1) is get_source(v2)
  assert not get_duplicate_from(v1, None)
  assert not get_duplicate_from(v2, None)
  if v1 is v2: return v1
  x1 = get_outject(v1)
  x2 = get_outject(v2)
  a1 = is_accepted(x1)
  a2 = is_accepted(x2)
  if a1 and not a2: return v1
  if a2 and not a1: return v2
  if not a1:
    return v1                   # Arbitrary
  assert not x1 is x2
  m = simple.mrca(x1, x2)
  # Choose the most tipward taxon
  if m == x1: return v2
  if m == x2: return v1
  # v1 and v2 are disjoint??  How can this happen?  Arbitrary I guess
  # Prefer the one that has children
  c1 = len(get_children(x1, ()))
  c2 = len(get_children(x2, ()))
  if c1 == 0 and c2 > 0: return c2
  if c2 == 0 and c1 > 0: return c1
  log("# Which is better as a type? %s %s" % (blorb(v1), blorb(v2)))
  assert False
  return v1

# Access to specimen records [sid, u, v]

def sid_to_specimen(AB, sid):
  return AB.specimen_ufs[sid]

def sid_to_record(AB, sid, z):
  uf = sid_to_specimen(AB, sid)
  (_, u, v) = uf.payload()
  return u if isinA(AB, z) else v

def sid_to_opposite_record(AB, sid, z):
  uf = sid_to_specimen(AB, sid)
  (_, u, v) = uf.payload()
  return v if isinA(AB, z) else u

def sid_to_epithet(AB, sid):
  uf = sid_to_specimen(AB, sid)
  (_, u, v) = uf.payload()
  if u:
    parts = get_parts(get_outject(u))
  else:
    parts = get_parts(get_outject(v))
  return parts.epithet or parts.genus 

# -----------------------------------------------------------------------------

# u = smallest taxon in A checklist containing this specimen
# v = smallest taxon in B checklist containing this specimen

# typification_uf = specimen that is a type of some taxon.
# We might just call these "types"

(maybe_get_typification, set_typification_uf) = \
  prop.get_set(prop.declare_property("typification_uf"))

# Only workspace nodes have uf records

def get_typification(u):
  probe = maybe_get_typification(u, None)
  if probe: return probe
  AB = get_workspace(u)
  r = [None, u, None] if isinA(AB, u) else [None, None, u]
  uf = UnionFindable(r)
  set_typification_uf(u, uf)
  return uf

# Are u and v known to have the same typification (type specimen)?

def same_typifications(u, v):
  uf = maybe_get_typification(u, None)
  vf = maybe_get_typification(v, None)
  return uf and vf and same_specimens(uf, vf)

# u and v are in workspace but may or may not be from same checklist

def equate_typifications(u, v):     # opposite checklists. u might be species
  if u is not v:
    equate_specimens(get_typification(u), get_typification(v))
    if monitor(u) or monitor(v):
      log("# Unifying exemplar(%s) = exemplar(%s)" % (blorb(u), blorb(v)))
  return u

# Are u and v equatable?
# u and v are assumed to be in the same checklist.
# This implements transitivity of equality, I think?

def equatable_typifications(u, v):
  uf = maybe_get_typification(u, None)
  if uf:
    vf = maybe_get_typification(v, None)
    if vf:
      return uf.payload() is vf.payload()
  return True

# For debugging mainly?

def get_typifies(spec):
  (sid, u, v) = spec.payload()
  return v or u

# -----------------------------------------------------------------------------
# An exemplar is a specimen that's the typification of two or more taxa in distinct
# checklists.

# Returns typification record (sid, u, v) or None.
# z is in AB.

def is_exemplar(uf):
  if uf:
    (_, u, v) = uf.payload()
    return u and v
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
