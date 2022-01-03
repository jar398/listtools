#!/usr/bin/env python3

import sys, csv, argparse
import util, property as prop
import checklist, workspace, merge
import theory

from util import windex, MISSING
from property import mep_get, mep_set
from rcc5 import *
from checklist import *
from workspace import *

from merge import reverse_note

get_difference = prop.getter(merge.difference_prop)

#

def reports(merged_iter, mode):
  AB = rows_to_merged(merged_iter, {'name': 'AB'})
  resolve_merge_links(AB)
  for r in all_records(AB): assert get_superior(r)  # sanity
  ensure_inferiors_indexed(AB)
  #theory.ensure_levels(AB)
  return generate_reports(AB, mode)

readable_template = "%s"

# Full mode shows every concept in the sum.
# MDD mode only shows changed/new/removed rows.

def generate_reports(AB, mode):
  stats = {}

  def generate_report(AB, mode, stats):
    yield ("A name", "B name", "B accepted", "rank", "differences", "basis of match")

    def tick(stat, y, rank):
      if stat in stats:
        s = stats[stat]
      else:
        s = [0, 0]      # [all, species only]
        stats[stat] = s
      s[0] += 1
      if rank == "species":
        s[1] += 1

    for z in all_records(AB):
      (x, y) = personae(z)
      if x or y:   #(x and is_accepted(x)) or (y and is_accepted(y)):
        diff = get_difference(z, None)
        rank = ((y and get_rank(y, None)) or (x and get_rank(x, None)))
        tick(diff or 'same', z, rank)
        if mode == 'full':
          doit = True
        elif mode == 'mdd':
          doit = (diff and rank == 'species')
        else:
          assert False   # invalid mode
        if doit:
          accepted = ''
          if y:
            sup = get_superior(y)
            if sup and sup.relationship != ACCEPTED:
              accepted = blurb(sup.record)
          yield (blurb(x) if x else '',
                 blurb(y) if y else '',
                 accepted,
                 rank, diff,
                 get_basis_of_match(y, '') if y else '')

  def generate_summary(stats):
    for key in sorted(list(stats.keys())):
      s = stats[key]
      yield("%7d %7d %s" % (s[0], s[1], key))

  return (generate_report(AB, mode, stats), generate_summary(stats))

def is_accepted(x):
  sup = get_superior(x, None)
  return (not sup) or sup.relationship == ACCEPTED

def personae(z):
  if is_top(z):
    return (None, None)
  if get_equated(z, None):
    return (get_equated(z).record, z)  # B with equated A
  else:
    rs = get_superior(z)
    if rs.relationship == EQ:        # A equated to B - ignore
      return (None, None)
    else:
      diff = get_difference(z, '')
      # this is a kludge. there ought to be a more direct test
      if 'not in A' in diff:
        return (None, z)
      elif 'not in B' in diff:
        return (z, None)
      else:
        assert False

# -----------------------------------------------------------------------------

# Load a merged checklist, recovering equated, match, and
# maybe... left/right id links? ... ...
# N.b. this will be a mere checklist, not a workspace.

def rows_to_merged(rows, meta):
  M = checklist.rows_to_checklist(rows, meta)  # sets superiors
  resolve_merge_links(M)
  return M

# This is used for reporting (not for merging).
# When we get here, superiors have already been set, so the ghosts
# will be bypassed in enumerations (preorder and all_records).

def resolve_merge_links(M):
  conjured = {}
  def register(rec):
    pk = get_primary_key(rec)
    if pk in conjured: log("** duplicate: %s" % pk)
    conjured[pk] = rec
  for record in all_records(M):

    note = get_equated_note(record, None)
    key = get_equated_key(record, None)
    if key != None:
      # Only B records have equated's:
      #  have equation record EQ z, want synonym z EQ record
      z = look_up_record(M, key, record) or conjured.get(key, None)
      if not z:
        z = conjure_record(M, key, record, register)   # in A
        set_superior(z, relation(EQ, record, "equivalent", note))
        if False:
          # big fragile kludge
          sup_level = get_level(record, None)
          assert sup_level != None, blurb(record)
          set_level(z, sup_level+1)
      else:
        # If z already exsts, it will take care of its own synonymity
        pass
      set_equated(record, relation(EQ, z, "equivalent", note))  # B->A

    note = get_basis_of_match(record, None)
    key = get_match_key(record, None)
    if key != None:
      # B: record MATCHES z, A: z MATCHES record
      z = look_up_record(M, key, record) or conjured.get(key, None)
      if not z:
        z = conjure_record(M, key, record, register)
        set_match(z, relation(EQ, record, "nominal match",
                              reverse_note(note)))
      else:
        # If z already exists, it will take care of its own match
        pass
      set_match(record, relation(EQ, z, "nominal match", note))
    elif note:
      set_match(record, relation(NOINFO, None, None, note))
  return conjured

# Create new dummy records if needed so that we can link to them.
# key is for what record matches and is unresolved.  That means that the
# record was a B record, and key is for a suppressed A record that 
# we need to 'conjure'.

def conjure_record(C, key, record, register):
  #log("# Creating record %s by demand from %s" % (key, blurb(record)))
  rec = checklist.make_record(key, C)
  register(rec)
  return rec

#      x0 ---
#             \   match
#              \
#       x  <->  y    SYNONYM / EQ
#        \
#  match  \
#           --- y0


# -----------------------------------------------------------------------------

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    standard input = the merged checklist we're reporting on
    """)
  parser.add_argument('--mode',
                      default='mdd',
                      help="report type: 'mdd' or 'full'")
  parser.add_argument('--summary',
                      help="where to write the summary report (path)")
  args=parser.parse_args()

  (report, summary) = reports(csv.reader(sys.stdin), args.mode)
  util.write_rows(report, sys.stdout)

  summary = list(summary)  # need it 2x
  log("# Summary has %s rows" % len(summary))
  for line in summary: print(line, file=sys.stderr)
  if args.summary:
    with open(args.summary, "w") as sumfile:
      for line in summary: print(line, file=sumfile)
