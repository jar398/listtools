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
  sid_to_specimen, get_typifies, same_specimens
from proximity import near_enough

from homotypy import compare_parts, MOTION, REVIEW, explain_classified

# Problem perhaps: This creates specimen objects even when they aren't matched.

def find_endohomotypics(AB):
  homotypy.find_homotypics_in_checklist(AB)
  homotypy.find_homotypics_in_checklist(swap(AB))

# Replaced by version in specimen.py
def find_homotypics_in_checklist(AB):
  def process(x):

    epithets = {}
    def consider(c):
      c_ep = get_parts(c).epithet
      if c_ep:
        if c_ep in epithets:
          equate_type_ufs(AB.in_left(c), epithets[c_ep])
        else:
          epithets[c_ep] = AB.in_left(c)

    # Match siblings
    for c in get_inferiors(x):
      process(c)
      classified = relate_records(AB.in_left(c), AB.in_left(x))
      if classified >= MOTION:
        consider(c)

    # Match parent/child
    consider(x)

  process(AB.A.top)

# This can be configured to run once or run twice.  'last' means we're
# on the last pass so if there was a first pass, proximity is available.

# Unify the types of A with the types of B.

def find_type_ufs(AB, subprobs, get_estimate, last):
  # This sets the 'type_uf' property of ... some ... records.

  n = 1
  for (key, (us, vs)) in subprobs.items():  # For each subproblem
    if n % 1000 == 0 or PROBE in key:
      log("# Subproblem %s %s %s %s %s" %
          (n, len(us), len(vs), blurb(us[0]), blurb(vs[0])))
    # us and vs are sorted by 'unimportance' (i.e. most important first)
    n += 1
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
        # relate_records checks proximity
        observe_match(u, v, u_matches, classified)
        observe_match(v, u, v_matches, classified)
        if monitor(u):
          log("# observe %s => %s = %s" %
              (blurb(u), blurb(v), explain_classified(classified)))
      # end j loop
    # end i loop
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
            equate_specimens(u_spec, v_specs[0])
            equate_specimens(u_spec, v_specs[1])
      else:
        v_spec = v_specs[0]
        v_sid = get_specimen_id(v_spec)
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
            # WRONG? - only do this for 'tipward' nodes?
            equate_specimens(u_spec, v_spec)
        else:
          log("# No return match: %s -> %s" %   # No return match for v - shouldn't happen
              (blorb(get_typifies(u_spec)),
               blorb(get_typifies(v_spec))))
      # End u_sid loop.
      # end subproblem loop

  equate_type_ufs(AB.in_left(AB.A.top),
                       AB.in_right(AB.B.top))

# u_matches : u_sid -> (v_clas, v_specs)
#   'spec' is short for 'specimen'

# Ideally, only one match per specimen id.
# 'Classified' might be a failure e.g. HETEROTYPIC

def observe_match(u, v, u_matches, classified):
  u_spec = get_type_uf(u)
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

# Compare potentially homotypic taxa in same workspace.

def relate_records(u, v):
  classified = compare_parts(u, v)

  if classified == MOTION:
    if not near_enough(u, v):
    # Proximity within the hierarchy is essential if no genus match
      return REVIEW

  if monitor(u) or monitor(v):
    log("# Compare %s, %s = %s" %
        (blurb(u), blurb(v), explain_classified(classified)))
  return classified

