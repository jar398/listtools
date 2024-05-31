#!/usr/bin/env python3

from parse import PROBE

import sys, argparse
import util, rows

from util import log, windex
from workspace import *

from checklist import get_duplicate_from, set_duplicate_from, blorb, blurb
from specimen import get_exemplar, get_exemplar_id, sid_to_epithet
from specimen import equate_specimens, equate_typifications, \
  get_typification, maybe_get_typification

from estimate import find_estimates, get_estimate
from typify import find_typifications
from typify import basic_typifications, get_family

# listtools's exemplar-finding procedure.  If there is some other way
# of finding exemplars, that's fine, don't need to use this.

# This is invoked twice - two-pass method.  Purpose of first pass is
# to be able to compute distances on the second pass.

def find_exemplars(get_estimate, AB):
  basic_typifications(AB.A)
  basic_typifications(AB.B)
  match_typifications(AB)
  # maybe compute better estimates - see theory.py
  report_on_exemplars(AB)

# ------

def report_on_exemplars(AB):
  count = ufcount = 0      # of nodes having exemplars?
  
  # but we could just look at AB.specimen_ufs, instead?
  for x in preorder_records(AB.A):
    u = AB.in_left(x)
    uf = maybe_get_typification(u, None)
    if uf:
      ufcount += 1
      b = get_exemplar(u)        # forces sid assignment, return (sid,u,v) ?
      if b:
        count += 1
        get_exemplar_id(uf)        # forces sid assignment  ??
  log("* Nodes with typification: %s, nodes with exemplars: %s, specimen ids: %s" %
      (ufcount, count, len(AB.specimen_ufs)))

def write_exemplar_list(AB, out=sys.stdout):
  util.write_rows(generate_exemplars(AB), out)

def generate_exemplars(AB):
  yield ("exemplar id", "epithet", "checklist", "taxonID", "canonicalName", "duplicate of")
  count = [0]
  rows = []
  def doit(ws, which):
    rcount = ecount = 0
    for x in preorder_records(ws.A):
      rcount += 1
      if not is_top(x):               # ?
        u = ws.in_left(x)
        uf = get_exemplar(u)     # exemplar record [sid, u, v] or None
        if uf:
          ecount += 1
          sid = get_exemplar_id(uf)
          epithet = sid_to_epithet(AB, sid)
          dup = get_duplicate_from(x, None)
          rows.append((sid, epithet, which, get_primary_key(x), get_canonical(x),
                       get_primary_key(dup) if dup else MISSING))
          count[0] += 1
    log("# preorder: %s, exemplars: %s" % (rcount, ecount)) 
  doit(swap(AB), 'B')
  doit(AB, 'A')
  rows.sort(key=lambda row:(row[0], row[2], row[4]))
  yield from rows
  log("# %s rows in exemplar report" % count[0])

# Read list of exemplars from file (given as a Rows)

def read_exemplars(in_rows, AB):
  equate_typifications(AB.in_left(AB.A.top), AB.in_right(AB.B.top))
  the_rows = in_rows.rows()     # caller will close in_rows
  header = next(the_rows)
  sid_col = windex(header, "exemplar id")
  which_col = windex(header, "checklist")
  taxonid_col = windex(header, "taxonID")
  dup_col = windex(header, "duplicate of")
  for row in the_rows:
    taxonid = row[taxonid_col]
    which = row[which_col]
    if which == 'A':
      C = AB.A
    elif which == 'B':
      C = AB.B
    else:
      log("# Invalid checklist indicator %s" % which)
    dup_id = row[dup_col]
    x = checklist.look_up_record(C, taxonid)
    if not x:
      log("## read_exemplars: Record not found?! %s" % taxonid)
    else:
      u = AB.in_left(x) if which=='A' else AB.in_right(x)
      if monitor(u):
        log("# Found %s %s with dup id '%s'" % (taxonid, blurb(u), dup_id))

      if dup_id:
        drec = checklist.look_up_record(C, dup_id)
        assert drec
        set_duplicate_from(x, drec)
        if monitor(u):
          log("# %s %s is duplicated from %s %s" %
              (taxonid, blurb(u), dup_id, blurb(drec)))

      uf = get_typification(u)

      # row is (sid, which, taxonid)
      sid = int(row[sid_col])
      #log("# read exemplar #%s in %s = %s" % (sid, which, taxonid))
      if sid in AB.specimen_ufs: # as a key
        s = AB.specimen_ufs[sid]
        a = equate_specimens(s, uf)
      else:
        AB.specimen_ufs[sid] = uf
        a = uf
      a.payload()[0] = sid


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
      find_exemplars(get_estimate, AB)
      write_exemplar_list(AB)
