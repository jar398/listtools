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
from cluster import find_exemplars

# Write exemplars to a file

def gather_exemplars(AB):
  return cluster.find_exemplars(AB.A, AB.B)  # list of exemplars

def write_exemplar_list(exemplars, AB, out=sys.stdout):
  util.write_rows(generate_exemplars(exemplars), out)

def generate_exemplars(exemplars):
  yield ("exemplar id", "epithet", "checklist", "taxonID", "canonicalName")
  count = 0

  for exemplar in exemplars:
    for record in left_cluster(examplar):
      yield exemplar_id, None, "A", get_primary_key(record), get_canonical(record)
      count += 1
    for record in right_cluster(examplar):
      yield exemplar_id, None, "B", get_primary_key(record), get_canonical(record)
      count += 1

  log("# %s rows in exemplar report" % count)

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
      write_exemplar_list(gather_exemplars(AB), AB)
