# Identify homotypies that tie checklists together

import util
import simple
import homotypy

from parse import PROBE

from util import log, MISSING
from rcc5 import DISJOINT
from checklist import get_parts, monitor, \
  get_inferiors, \
  is_accepted, blurb, blorb, get_scientific, get_primary_key, get_rank, \
  get_nomenclatural_status
from workspace import separated, get_outject, get_workspace, local_sup
from workspace import isinA, isinB, local_accepted, \
  swap
from simple import compare_per_checklist
from specimen import sid_to_epithet, same_specimens, \
  equate_type_ufs, equate_specimens, \
  get_type_uf, get_specimen_id, \
  sid_to_specimen, get_typifies, same_specimens, \
  maybe_get_type_uf

from homotypy import compare_parts, MOTION, REVIEW, explain_classified
from homotypy import relate_records

# --------------------

# listtools's exemplar-finding procedure.  If there is some other way
# of finding exemplars, that's fine, don't need to use this.

def find_exemplars(AB):
  find_endohomotypics(AB)       # Within each checklist

  subproblems = find_subproblems(AB)
  log("* Finding type_ufs (single pass):")
  find_type_ufs(AB, subproblems, None, True)

  # maybe compute better estimates - see theory.py
  report_on_exemplars(AB)

# Find blocks/chunks, one per epithet

def find_subproblems(AB):
  log("* Finding subproblems:")
  (A_index, B_index) = \
    map(lambda CD: \
        index_by_some_key(CD,
                          # should use genus if epithet is missing
                          get_subproblem_key),
        (AB, swap(AB)))
  subprobs = {}
  for (key, us) in A_index.items():
    assert key != MISSING, blurb(us[0])
    vs = B_index.get(key, None)
    if vs != None:
      if any(map(monitor, us)) or any(map(monitor, vs)):
        log("* Found monitored subproblem %s / %s" %
            (map(blurb, us), map(blurb, vs)))
      us.sort(key=unimportance)
      vs.sort(key=unimportance)
      subprobs[key] = (us, vs)
      if (any(monitor(u) for u in us) or
          any(monitor(v) for v in vs)):
        log("* Added subproblem %s e.g. %s.  %s x %s" %
            (key, blurb(us[0]), len(list(map(blurb, us))), len(list(map(blurb, vs)))))
    else:
      if PROBE in key:
        log("# Null subproblem %s" % key)
  log("* There are %s subproblems." % len(subprobs))
  AB.subproblems = subprobs
  return subprobs

# More important -> lower number, earlier in sequence

def unimportance(u):
  parts = get_parts(u)
  if parts.epithet == MISSING: unimp = 4      # Foo
  elif parts.middle == parts.epithet: unimp = 2     # Foo bar bar
  elif parts.middle == None or parts.middle == '':  unimp = 0     # Foo bar
  else: unimp = 1                         # Foo bar baz
  x = get_outject(u)
  return (1 if is_accepted(x) else 2,
          unimp,
          # Prefer to match the duplicate that has children
          # (or more children)
          -len(get_children(x, ())),
          -len(get_synonyms(x, ())),
          get_scientific(x, None),
          get_primary_key(x))

# Returns dict value -> key
# fn is a function over AB records

def index_by_some_key(AB, fn):
  index = {}
  for x in postorder_records(AB.A):
    u = AB.in_left(x)
    key = fn(u)
    if False and monitor(u):
      log("# 1 Indexing %s, key %s, monitored? %s" % (blurb(u), key, monitor(u)))
      log("# 2 Indexing %s, key %s, monitored? %s" % (blurb(u), key, monitor(x)))
    #assert key  - MSW has scientificName = ?
    have = index.get(key, None) # part
    if have:
      have.append(u)
    else:
      index[key] = [u]
  return index

# Each subproblem covers a single epithet (or name, if higher taxon)
# z is in AB

def get_subproblem_key(z):
  x = get_outject(z)
  parts = get_parts(x)
  ep = parts.epithet            # stemmed
  key = ep if ep else parts.genus
  if key:
    if False and monitor(x):
      log("# Subproblem key is %s for %s" % (key, blurb(x)))
  else:
    log("** %s: Name missing or ill-formed: %s" %
        (get_primary_key(x), parts,))
    key = '?' + get_primary_key(x)
  return key

# ------

def report_on_exemplars(AB):
  count = ufcount = 0      # of nodes having exemplars?
  
  # but we could just look at AB.specimen_ufs, instead?
  for x in preorder_records(AB.A):
    u = AB.in_left(x)
    uf = maybe_get_type_uf(u, None)
    if uf:
      ufcount += 1
      b = get_exemplar(u)        # forces sid assignment, return (sid,u,v) ?
      if b:
        count += 1
        get_exemplar_id(uf)        # forces sid assignment  ??
  log("# Nodes with type specimens: %s, nodes with exemplars: %s, specimen id UFs: %s" %
      (ufcount, count, len(AB.specimen_ufs)))

# -----------------------------------------------------------------------------

# Problem perhaps: This creates specimen objects even when they aren't matched.

def find_endohomotypics(AB):
  homotypy.find_homotypics_in_checklist(AB)
  homotypy.find_homotypics_in_checklist(swap(AB))

# Unify the type specimens of A with the type specimens of B.

def find_type_ufs(AB, subprobs, get_estimate, last):
  # This sets the 'type_uf' property of ... some ... records.

  n = 1
  for (key, (us, vs)) in subprobs.items():  # For each subproblem
    if n % 1000 == 0 or PROBE in key:
      log("# Subproblem %s %s %s %s %s" %
          (n, len(us), len(vs), blurb(us[0]), blurb(vs[0])))
    # us and vs are sorted by 'unimportance' (i.e. most important first)
    n += 1

    (u_matches, v_matches) = find_matches(key, us, vs)

    # u_specs and v_specs are sorted by 'unimportance'
    # u_matches : u_sid -> (v_clas, v_specs)
    # clas means a classification
    for (u_sid, (v_clas, v_specs)) in u_matches.items():
      u_spec = sid_to_specimen(AB, u_sid)
      # if monitor(get_typifies(u_spec)): ...
      if v_clas == MOTION:
        # Maybe log AUTHORLESS and/or REVIEW
        # Put this in the report somehow ??
        u0 = get_typifies(u_spec)
        v0 = get_typifies(v_specs[0])
        log("# %s: %s -> %s" %  # or, make a note of it for review
            (explain_classified(v_clas), blorb(u0), blorb(v0)))
      if v_clas < MOTION:
        pass
      elif len(v_specs) > 1:
        # Problem here
        # But it seems to not occur (in MDD/COL).
        v0 = get_typifies(v_specs[0]) # record that has given specimen
        v1 = get_typifies(v_specs[1])
        if is_accepted(v0) and is_accepted(v1):
          if compare_per_checklist(get_outject(v0), get_outject(v1)) \
             .relationship & DISJOINT != 0:
            # Another kind of REVIEW
            # TBD: suppress if synonym?
            log("# B ambiguity: %s %s -> %s %s, %s %s" %
                (explain_classified(v_clas),
                 blorb(get_typifies(u_spec)),
                 get_primary_key(v0),
                 blorb(v0),
                 get_primary_key(v1),
                 blorb(v1)))
          else:
            # DISJOINT is impossible
            equate_specimens(u_spec, v_specs[0])
            equate_specimens(u_spec, v_specs[1])
      else:
        v_spec = v_specs[0]
        v_sid = get_specimen_id(v_spec) # invents new id on demand
        # v_matches : v_sid -> (u_clas, u_specs)
        results = v_matches.get(v_sid)
        if results:
          (u_clas, u_specs) = results
          if u_clas == MOTION:
            # also maybe AUTHORLESS, REVIEW ?
            log("# Review 2 %s -> %s" %         # or, make a note of it for review
                (blorb(get_typifies(u_spec)),
                 blorb(get_typifies(v_spec)),))
          if u_clas < MOTION:    # e.g. REVIEW
            pass
          elif len(u_specs) > 1:
            # Select type A b b from u_specs if possible...
            u0 = get_typifies(u_specs[0])
            u1 = get_typifies(u_specs[1])
            if compare_per_checklist(get_outject(u0), get_outject(u1)) \
               .relationship & DISJOINT != 0:
              # TBD: suppress if synonym?
              log("# A ambiguity %s: %s, %s -> %s" %
                  (explain_classified(u_clas),
                   blorb(u0),
                   blorb(u1),
                   blorb(get_typifies(v_spec))))
            else:
              # Can be simultaneously compatible with both.
              # (TBD: deal with 3-way or more ambiguity.)
              equate_specimens(u_specs[0], v_spec)
              equate_specimens(u_specs[1], v_spec)
          elif not same_specimens(u_specs[0], u_spec):
            # How can this happen ???
            log("# Vee %s ->\n  %s -> %s" %
                (blorb(get_typifies(u_spec)),
                 blorb(get_typifies(v_spec)),
                 blorb(get_typifies(u_specs[0])))),
          else:
            equate_specimens(u_spec, v_spec)
        else:
          log("# No return match: %s -> %s" %   # No return match for v - shouldn't happen
              (blorb(get_typifies(u_spec)),
               blorb(get_typifies(v_spec))))
      # End u_sid loop.
      # end subproblem loop

  equate_type_ufs(AB.in_left(AB.A.top),
                       AB.in_right(AB.B.top))

def find_matches(key, us, vs):
  if any(map(monitor, us)):
    log("# Doing subproblem %s" % key)
    log("#  us = %s" % list(map(blurb, us),))
    log("#  vs = %s" % list(map(blurb, vs),))

  u_matches = {}    # u_sid -> (class, [v_spec, ...])
  v_matches = {}    # v_sid -> (class, [u_spec, ...])
                    # spec = (sid, u, v)
                    # Didn't count on multiple specs per sid!
  for i in range(0, len(us)):
    u = us[i]
    if monitor(u): log("# Subproblem row: '%s' '%s'" % (key, blorb(u)))
    for j in range(0, len(vs)):
      v = vs[j]
      classified = relate_records(u, v) # shows!
      if classified >= MOTION:
        # relate_records checks proximity
        observe_match(u, v, u_matches, classified)
        observe_match(v, u, v_matches, classified)
        if monitor(u):
          log("# observe %s => %s = %s" %
              (blurb(u), blurb(v), explain_classified(classified)))
    # end j loop
  # end i loop
  return u_matches, v_matches

# u_matches : u_sid -> (v_clas, v_specs)
#   'spec' is short for 'specimen'

# Ideally, only one match per specimen id.
# 'Classified' might be a failure e.g. HETEROTYPIC

def observe_match(u, v, u_matches, classified):
  u_spec = get_type_uf(u) # there may be multiple u's with same specimen
  v_spec = get_type_uf(v)

  u_sid = get_specimen_id(u_spec)
  v_sid = get_specimen_id(v_spec)

  have = u_matches.get(u_sid)   # (v_clas, [v_speq, ...])
  if have:
    (v_clas, v_speqs) = have
    if classified > v_clas:
      # Discard previously seen lower value matches.
      if monitor(u):
        log("* Discarding previously seen lower value matches: %s" %
            blurb(u))
      u_matches[u_sid] = (classified, [v_spec])
    elif classified != v_clas:
      if monitor(u):
        log("* Discarding lower value match: %s" %
            blurb(u))
      pass
    elif any(map(lambda v_speq: same_specimens(v_spec, v_speq),
                 v_speqs)):
      # Don't add same specimen twice
      if monitor(u):
        log("* Discarding redundant match: %s" %
            blurb(u))
      pass
    else:
      if monitor(u):
        log("* Adding match: %s" % blurb(u))
      v_speqs.append(v_spec)    # ????
  else:
    u_matches[u_sid] = (classified, [v_spec])

def collate_by_sid(us):
  dict = {}
  for u in us:
    sid = get_specimen_id(u, None)
    seen = dict[sid]
    if sid != None:
      seen.append(u)
    else:
      seen = [u]
  return dict
