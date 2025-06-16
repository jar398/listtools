#!/usr/bin/env python3

# Assumes cross-checklist type specimen matches have been made

from parse import PROBE

import sys, argparse
import util, rows

from util import log, windex
from workspace import *

from checklist import blorb, blurb

# Obsolete 
from specimen import get_exemplar, get_exemplar_id, sid_to_epithet
from specimen import equate_specimens, equate_type_ufs, \
  get_type_uf

import cluster

# Write exemplars to a file

def write_exemplar_list(AB, out=sys.stdout):
  util.write_rows(generate_exemplars(AB), out)

def generate_exemplars(AB):
  yield ("exemplar id", "epithet", "checklist", "taxonID", "canonicalName")
  count = [0]
  rows = []
  def doit(ws, which):
    rcount = ecount = 0
    for x in preorder_records(ws.A):
      rcount += 1
      if not is_top(x):               # ?
        u = ws.in_left(x)
        exem = get_exemplar(u)     # exemplar record [sid, u, v] or None
        if exem:
          ecount += 1
          sid = get_exemplar_id(exem)
          epithet = sid_to_epithet(AB, sid)
          rows.append((sid, epithet, which, get_primary_key(x), get_canonical(x)))
          count[0] += 1
    log("# preorder: %s, exemplars: %s" % (rcount, ecount)) 
  doit(swap(AB), 'B')
  doit(AB, 'A')
  rows.sort(key=lambda row:(row[0], row[2], row[4]))
  yield from rows
  log("# %s rows in exemplar report" % count[0])

# Read list of exemplars from file (given as a Rows)

def read_exemplars(in_rows, AB):
  assert len(AB.specimen_ufs) == 0
  equate_type_ufs(AB.in_left(AB.A.top), AB.in_right(AB.B.top))
  the_rows = in_rows.rows()     # caller will close in_rows
  header = next(the_rows)
  sid_col = windex(header, "exemplar id")
  which_col = windex(header, "checklist")
  taxonid_col = windex(header, "taxonID")
  for row in the_rows:
    taxonid = row[taxonid_col]
    which = row[which_col]
    sid = int(row[sid_col])
    if which == 'A':
      C = AB.A
    elif which == 'B':
      C = AB.B
    else:
      log("# Invalid checklist indicator %s" % which)
    x = checklist.look_up_record(C, taxonid)
    if not x:
      log("## read_exemplars: Record not found?! %s" % taxonid)
    else:
      if sid > AB.max_sid:
        AB.max_sid = sid
      u = AB.in_left(x) if which=='A' else AB.in_right(x)
      exem = get_type_uf(u)        # one exemplar for each taxon

      if sid in AB.specimen_ufs: # as a key
        uf2 = AB.specimen_ufs[sid]
        assert uf2.payload()[0] == sid
        exem = equate_specimens(exem, uf2)

      exem.payload()[0] = sid
      AB.specimen_ufs[sid] = exem # Replace


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Generate list of exemplars proposed for two checklists
    """)
  parser.add_argument('--A', help="the A checklist, as path name or -",
                      default='-')
  parser.add_argument('--B', help="the B checklist, as path name or -",
                      default='-')
  parser.add_argument('--Aname', help="short name of the A checklist",
                      default='A')
  parser.add_argument('--Bname', help="short name of the B checklist",
                      default='B')
  args=parser.parse_args()

  a_name = args.Aname
  b_name = args.Bname
  a_path = args.A
  b_path = args.B
  with rows.open(a_path) as a_rows:
    with rows.open(b_path) as b_rows:
      # compute name matches afresh
      AB = ingest_workspace(a_rows.rows(), b_rows.rows(),
                            A_name=a_name, B_name=b_name)
      find_exemplars(AB)
      write_exemplar_list(AB)
