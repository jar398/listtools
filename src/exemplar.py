#!/usr/bin/env python3

from parse import PROBE

import sys, argparse
import util, rows

from workspace import *
from util import log, windex

from estimate import find_estimates, get_estimate

from typify import equate_typification_ufs
from typify import get_exemplar, get_typification_uf
from typify import get_typification_record, equate_typifications
from typify import really_get_typification_uf, find_typifications
from typify import unimportance

# listtools's exemplar-finding procedure.  If there is some other way
# of finding exemplars, that's fine, don't need to use this.

# This is invoked twice - two-pass method.  Purpose of first pass is
# to be able to compute distances on the second pass.

def find_exemplars(get_estimate, AB):
  log("# Finding subproblems")
  subproblems = find_subproblems(AB)
  log("#   There are %s subproblems" % len(subproblems))
  log("# Finding pass 1 typifications (for distance calculations)")
  find_typifications(AB, subproblems, None, False)
  find_estimates(AB)            # for distance calculations
  log("# Finding pass 2 typifications (using distance calculations)")
  find_typifications(AB, subproblems, get_estimate, True)
  # maybe compute better estimates - see theory.py
  report_on_exemplars(AB)

# Find blocks/chunks, one per epithet

def find_subproblems(AB):
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
      us.sort(key=unimportance)
      vs.sort(key=unimportance)
      subprobs[key] = (us, vs)
      if (any(monitor(u) for u in us) or
          any(monitor(v) for v in vs)):
        log("# Subproblem %s" % key)
        log("#  us = %s" % (list(map(blurb, us)),))
        log("#  vs = %s" % (list(map(blurb, vs)),))
    else:
      if PROBE in key:
        log("# Null subproblem %s" % key)
  AB.subproblems = subprobs
  return subprobs

# Returns dict value -> key
# fn is a function over AB records

def index_by_some_key(AB, fn):
  index = {}
  for x in postorder_records(AB.A):
    u = AB.in_left(x)
    key = fn(u)
    if monitor(u):
      log("# Indexing %s %s" % (key, blurb(u)))
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
  parts = get_parts(get_outject(z))
  ep = parts.epithet            # stemmed
  key = ep if ep else parts.genus
  if key:
    if monitor(z):
      log("# Subproblem key is %s (%s)" % (key, blurb(z)))
  else:
    log("## Falsish name: %s" % (parts,))
    key = '?'
  return key

# ------

def report_on_exemplars(AB):
  count = ufcount = 0      # of nodes having exemplars?
  
  # but we could just look at AB.exemplar_ufs, instead?
  for x in preorder_records(AB.A):
    u = AB.in_left(x)
    uf = really_get_typification_uf(u, None)
    if uf:
      ufcount += 1
      b = get_typification_record(u)        # forces xid assignment, return (xid,u,v)
      if b:
        count += 1
        get_exemplar(u)        # forces xid assignment, return (xid,u,v)
  log("# Nodes with typification: %s, nodes with exemplars: %s, exemplars: %s" %
      (ufcount, count, len(AB.exemplar_ufs)))

def write_exemplar_list(AB, out=sys.stdout):
  util.write_rows(generate_exemplars(AB), out)

def generate_exemplars(AB):
  yield ("checklist", "taxonID", "exemplar id", "canonicalName")
  count = [0]
  def doit(ws, which):
    rcount = ecount = 0
    for x in preorder_records(ws.A):
      rcount += 1
      if not is_top(x):               # ?
        u = ws.in_left(x)
        r = get_exemplar(u)     # exemplar record [xid, u, v] or None
        if r:
          ecount += 1
          xid = r[0]
          yield (which, get_primary_key(x), xid, get_canonical(x))
          count[0] += 1
    log("# preorder: %s, exemplars: %s" % (rcount, ecount)) 
  yield from doit(AB, 0)
  yield from doit(swap(AB), 1)
  log("# %s rows in exemplar report" % count[0])

# Read list of exemplars from file (given as a Rows)

def read_exemplars(in_rows, AB):
  equate_typifications(AB.in_left(AB.A.top), AB.in_right(AB.B.top))
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
      uf = get_typification_uf(u)

      # row is xid, which, taxonid
      xid = int(row[xid_col])
      #log("# read exemplar #%s in %s = %s" % (xid, which, taxonid))
      if xid in AB.exemplar_ufs:
        a = equate_typification_ufs(AB.exemplar_ufs[xid], uf)
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
