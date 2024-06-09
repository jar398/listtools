
import util
import property as prop
import simple
import parse

from util import log, UnionFindable

from checklist import get_source, blurb, blorb, get_parts, monitor, \
  is_accepted, get_suppressed
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

# Prefer information coming from first argument

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

  # Formerly this captured "better" information
  r[1] = u1 or u2; r[2] = v1 or v2
  return ef

# Access to specimen records [sid, x, y]

def sid_to_specimen(source, sid):
  return source.specimen_ufs[sid]

def sid_to_record(source, sid):   # result is in A checklist
  uf = sid_to_specimen(source, sid)
  return uf.payload()[1]

def sid_to_opposite_record(source, sid):   # in B checklist
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
