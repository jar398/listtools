
import util
import property as prop
import simple
import parse

from util import log, UnionFindable

from checklist import get_source, blurb, blorb, get_parts, monitor, \
  is_accepted, get_redundant, get_canonical, get_rank, get_inferiors
from workspace import get_workspace, get_children, \
  get_outject, isinA
from ranks import ranks_dict
from rcc5 import DISJOINT
  
# Specimens per se

# Get a specimen id on demand (never returns falsish)
# Either an exemplar, or a specimen for which an id is needed for 
# indexing purposes (typify).

def get_specimen_id(uf):
  r = uf.payload()
  (sid, u, v) = r
  if sid == None:
    # Create specimen id (for set operations) on demand
    w = u or v
    assert w
    ws = get_workspace(w)
    sid = ws.max_sid + 1
    ws.max_sid = sid
    r[0] = sid
    ws.specimens[sid] = uf
    #log("# Specimen %s: (%s) <-> (%s)" % (sid, blorb(u), blorb(v)))
  return sid

def sid_to_specimen(AB, sid):
  return AB.specimens[sid]


# Assert identity of specimens

def same_specimens(uf, vf):
  return uf.find() is vf.find()

def equate_specimens(uf, vf):
  uf = uf.find()
  vf = vf.find()
  if uf is vf: return

  (i1, u1, v1) = uf.payload()
  (i2, u2, v2) = vf.payload()

  # Choose or improve a record for which this specimen is to be its type.
  # Prefer protonyms.
  t1 = _pick_type_taxon(u1, u2)
  t2 = _pick_type_taxon(v1, v2)

  # Not sure this is necessary, but should at least be harmless
  if i2 == None: i = i1
  elif i1 == None: i = i2
  else: i = min(i1, i2)

  ef = uf.absorb(vf)          # ef happens to be uf
  r = ef.payload()
  r[0] = i
  r[1] = t1
  r[2] = t2
  return ef

# This chooses, given two homotypic taxa, the one that best 'represents' 
# the type specimen.

def _pick_type_taxon(v1, v2):
  # See whether v2 is an improvement over v1 (lower unimportance value)
  if v1 == None: return v2
  if v2 == None: return v1
  assert get_source(v1) is get_source(v2)
  if v1 is v2: return v1
  x1 = get_outject(v1)
  x2 = get_outject(v2)
  assert not x1 is x2

  # Prefer species (A b) to (subspecies A b c and even to A b b)
  # ...

  # Prefer higher rank (e.g. species over subspecies)
  k1 = ranks_dict.get(get_rank(x1, None))
  k2 = ranks_dict.get(get_rank(x2, None))
  # lower numbers mean lower (more tipward) rank
  if k1 and k2:
    if k1 < k2: return v2       # higher/bigger taxon
    if k2 < k1: return v1

  # Prefer the most rootward taxon (backup in case ranks missing)
  m = simple.mrca(x1, x2)
  if m == x1 and m != x2: return v1
  if m == x2 and m != x1: return v2

  # Prefer the protonym, if any
  p1 = get_parts(x1).protonymp
  p2 = get_parts(x2).protonymp
  if p1 and not p2: return v1
  if p2 and not p1: return v2

  # Prefer the accepted record, if any
  a1 = is_accepted(x1)
  a2 = is_accepted(x2)
  if a1 and not a2: return v1
  if a2 and not a1: return v2

  # Prefer one with more children
  c1 = len(get_children(x1, ()))
  c2 = len(get_children(x2, ()))
  if c1 < c2: return v2
  if c1 > c2: return v1

  # log("# Difficult choice: %s vs. %s" % (blorb(v1), blorb(v2)))

  # Prefer the older one
  y1 = get_parts(x1).year
  y2 = get_parts(x2).year
  if y1 and y2:
    if y1 < y2: return v1
    if y1 > y2: return v2

  # Prefer earlier in alphabet ... ugh.
  n1 = get_canonical(x1)
  n2 = get_canonical(x2)
  if n1 < n2: return v1
  if n1 > n2: return v2

  log("# Don't know how to choose: %s %s" % (blurb(v1), blurb(v2)))
  return v1


# 0 is "good"
# A c c < A c < A b c

def typelike(x):
  # TBD: Deal with subgenus
  parts = get_canonical(x).split(' ')
  if len(parts) > 2 and parts[-1] == parts[-2]:
    return 0                    # A c c, best
  if len(parts) < 2:
    return 1                    # A c, better
  return 10                     # bad

# Access to specimen records

# sid -> payload
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

# type_uf = specimen that is the type of a given taxon. ???
# We might just call these "types".

(maybe_get_type_uf, set_type_uf) = \
  prop.get_set(prop.declare_property("type_uf"))

# The type specimen for a given record.
# Not necessarily an exemplar, not necessarily identifier-carrying.
# Only workspace nodes have uf nodes.

def get_type_uf(u):
  uf = maybe_get_type_uf(u, None)
  if not uf:
    AB = get_workspace(u)
    r = [None, u, None] if isinA(AB, u) else [None, None, u]
    uf = UnionFindable(r)
    set_type_uf(u, uf)
  return uf

# Are u and v known to have the same type_uf?
#  I.e. are they known to be homotypic?

def same_type_ufs(u, v):
  uf = maybe_get_type_uf(u, None)
  vf = maybe_get_type_uf(v, None)
  return uf and vf and same_specimens(uf, vf)

# u and v are in workspace but may or may not be from same checklist

def equate_type_ufs(u, v):     # opposite checklists. u might be species
  equate_specimens(get_type_uf(u), get_type_uf(v))
  if monitor(u) or monitor(v):
    log("# Unifying specimen(%s) = specimen(%s)" % (blorb(u), blorb(v)))
  return u

# Are u and v equatable?
# u and v are assumed to be in the same checklist.
# This implements transitivity of equality, I think?

def equatable_type_ufs(u, v):
  uf = maybe_get_type_uf(u, None)
  if uf:
    vf = maybe_get_type_uf(v, None)
    if vf:
      return uf.payload() is vf.payload()
  return True

# For debugging mainly?

def get_typifies(spec_uf):
  (sid, u, v) = spec_uf.payload()
  return v or u

# -----------------------------------------------------------------------------
# An exemplar is a specimen that's the type_uf of two or more taxa in distinct
# checklists.

# Returns type_uf record (sid, u, v) or None.
# z is in AB.

def is_exemplar(uf):
  if uf:
    (_, u, v) = uf.payload()
    return u and v
  return False

# Cf. analyze_blocks

def get_exemplar(z):            # Returns a uf node
  uf = maybe_get_type_uf(z, None)
  if is_exemplar(uf): return uf
  return None
  
def get_exemplar_info(z):
  uf = get_exemplar(z)
  if is_exemplar(uf): return uf.payload()
  return None

def get_exemplar_id(uf):
  if is_exemplar(uf):
    return get_specimen_id(uf)
  else:
    return None
