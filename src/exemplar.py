#!/usr/bin/env python3

import sys, argparse
import util, rows, linkage, some_exemplar

from workspace import *
from util import log, windex

from estimate import get_estimate, find_estimates
from some_exemplar import find_some_exemplars
from some_exemplar import get_exemplar, get_bare_exemplar, get_exemplar_uf
from some_exemplar import equate_exemplar_ufs, equate_exemplars

# listtools's exemplar-finding procedure.  If there is some other way
# of finding exemplars, that's fine, don't need to use this.

def find_exemplars(AB):
  log("# Finding subproblems")
  subproblems = linkage.find_subproblems(AB)
  log("#   There are %s subproblems" % len(subproblems))
  find_some_exemplars(AB, subproblems, lambda _u, _d: None, False)
  find_estimates(AB)
  find_some_exemplars(AB, subproblems, get_estimate, True)
  # maybe compute better estimates - see theory.py
  report_on_exemplars(AB)

def report_on_exemplars(AB):
  count = ufcount = 0      # of nodes having exemplars?
  
  # but we could just look at AB.exemplar_ufs, instead?
  for x in preorder_records(AB.A):
    u = AB.in_left(x)
    uf = some_exemplar.really_get_exemplar_uf(u, None)
    if uf:
      ufcount += 1
      b = get_bare_exemplar(u)        # forces xid assignment, return (xid,u,v)
      if b:
        count += 1
        get_exemplar(u)        # forces xid assignment, return (xid,u,v)
  log("# Nodes with union/find: %s, nodes with exemplars: %s, exemplars: %s" %
      (ufcount, count, len(AB.exemplar_ufs)))

def write_exemplar_list(AB, out=sys.stdout):
  util.write_rows(generate_exemplars(AB), out)

def generate_exemplars(AB):
  yield ("checklist", "taxonID", "exemplar id")
  count = [0]
  def doit(ws, which):
    rcount = ecount = lcount = 0
    for x in preorder_records(ws.A):
      rcount += 1
      if not is_top(x):               # ?
        if get_link(x, None): lcount += 1
        u = ws.in_left(x)
        r = get_exemplar(u)     # exemplar record [xid, u, v] or None
        if r:
          ecount += 1
          xid = r[0]
          yield (which, get_primary_key(x), xid)
          count[0] += 1
    log("# preorder: %s, links %s, exemplars: %s" % (rcount, lcount, ecount)) 
  yield from doit(AB, 0)
  yield from doit(swap(AB), 1)
  log("# %s rows in exemplar report" % count[0])

# Read list of exemplars from file (given as a Rows)

def read_exemplars(in_rows, AB):
  equate_exemplars(AB.in_left(AB.A.top), AB.in_right(AB.B.top))
  the_rows = in_rows.rows()     # caller will close in_rows
  header = next(the_rows)
  xid_col = windex(header, "exemplar id")
  which_col = windex(header, "checklist")
  taxonid_col = windex(header, "taxonID")
  for row in the_rows:
    which = int(row[which_col])
    taxonid = row[taxonid_col]
    C = AB.A if which==0 else AB.B
    x = checklist.look_up_record(C, taxonid)
    if x:
      u = AB.in_left(x) if which==0 else AB.in_right(x)
      uf = get_exemplar_uf(u)

      # row is xid, which, taxonid
      xid = int(row[xid_col])
      #log("# read exemplar #%s in %s = %s" % (xid, which, taxonid))
      if xid in AB.exemplar_ufs:
        a = equate_exemplar_ufs(AB.exemplar_ufs[xid], uf)
      else:
        AB.exemplar_ufs[xid] = uf
        a = uf
      a.payload()[0] = xid
  log("# %s exemplars" % len(AB.exemplar_ufs))


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Generate list of exemplars proposed for two checklists
    """)
  parser.add_argument('--A', help="the A checklist, as path name or -",
                      default='-')
  parser.add_argument('--B', help="the B checklist, as path name or -",
                      default='-')
  args=parser.parse_args()

  a_name = 'A'; b_name = 'B'
  a_path = args.A
  b_path = args.B
  with rows.open(a_path) as a_rows:
    with rows.open(b_path) as b_rows:
      # compute name matches afresh
      AB = ingest_workspace(a_rows.rows(), b_rows.rows(),
                            A_name=a_name, B_name=b_name)
      find_exemplars(AB)
      write_exemplar_list(AB)
