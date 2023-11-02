#!/usr/bin/env python3

import sys, argparse
import util, rows
from workspace import *
from util import log, windex

from estimate import get_estimate, find_estimates
from some_exemplar import find_some_exemplars, get_exemplar
from some_exemplar import init_exemplar

# listtools's exemplar-finding procedure.  If there is some other way
# of finding exemplars, that's fine, don't need to use this.

def find_exemplars(AB):
  find_some_exemplars(AB, get_estimate)
  find_estimates(AB)
  find_some_exemplars(AB, get_estimate)

# Read list of exemplars from file (given as a Rows)

def read_exemplars(in_rows, AB):
  exemplars = {}
  init_exemplar(AB, exemplars, 1, AB.A.top, AB.B.top)    # is there a better way?
  the_rows = in_rows.rows()     # caller will close in_rows
  header = next(the_rows)
  id_col = windex(header, "exemplar")
  a_col = windex(header, "A_taxonID")
  b_col = windex(header, "B_taxonID")
  for row in the_rows:
    id = int(row[id_col])
    log("# exemplar #%s = %s" % (id, row))
    x = checklist.look_up_record(AB.A, row[a_col])
    y = checklist.look_up_record(AB.B, row[b_col])
    init_exemplar(AB, exemplars, id, x, y)
  log("# %s exemplars" % len(exemplars))
  return exemplars

def write_exemplar_list(AB, out=sys.stdout):
  util.write_rows(generate_exemplars(AB), out)

def generate_exemplars(AB):
  yield ("exemplar", "A_taxonID", "B_taxonID", "A_blurb", "B_blurb") # representatives
  count = rcount = ecount = lcount = 0
  seen = {}
  for x in preorder_records(AB.A):
    if not is_top(x):               # ?
      rcount += 1
      z = AB.in_left(x)
      if get_link(x, None): lcount += 1
      r = get_exemplar(z)
      if r:
        ecount += 1
        (id, u, v) = r
        if u == z:
          u_key = get_primary_key(get_outject(u))
          v_key = get_primary_key(get_outject(v))
          yield (id, u_key, v_key, blurb(u), blurb(v))
          count += 1
  log("# %s rows in exemplar report" % count)
  log("# preorder: %s, links %s, exemplars: %s" % (rcount, lcount, ecount)) 


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
