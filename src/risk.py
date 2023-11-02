#!/usr/bin/env python3

import sys, csv, argparse
import rows, align, checklist, theory, exemplar
from workspace import ingest_workspace
from checklist import *

def generate_risk_report(AB):
  yield ("A name", "status", "kept", "split out", "lumped in")
  for x in preorder_records(AB.A):
    if get_rank(x, None) == 'species' and is_accepted(x):
      u = AB.in_left(x)
      v = get_link(u, None)
      if not v:
        yield (get_ok_name(x), "removed", "", "", "")
      else:
        b1 = theory.get_block(u)
        b2 = theory.get_block(v)
        keep  = b1.intersection(b2)
        flush = b1.difference(b2)
        add   = b2.difference(b1)
        k = not theory.is_empty_block(keep)
        f = not theory.is_empty_block(flush)
        a = not theory.is_empty_block(add)
        k_exemplars = exemplar_names(keep, AB, u, v)
        f_exemplars = exemplar_names(flush, AB, u, v)
        a_exemplars = exemplar_names(add, AB, u, v)
        if not k:
          #log("# No exemplars in common: %s" % blurb(u))
          yield (get_ok_name(x), "mismatch",
                 k_exemplars, f_exemplars, a_exemplars)
        elif f and a:
          log("# Reorg: %s" % blurb(u))
          yield (get_ok_name(x), "reorg",
                 k_exemplars, f_exemplars, a_exemplars)
        elif not f and not a:
          pass                  # No risk
        elif f:
          log("# Split: %s" % blurb(u))
          yield (get_ok_name(x), "split",
                 k_exemplars, f_exemplars, a_exemplars)
        else:
          log("# Lump: %s" % blurb(u))
          yield (get_ok_name(x), "lump",
                 k_exemplars, f_exemplars, a_exemplars)

def exemplar_names(xids, AB, u, v):
  def makename(xid):
    aname = get_ok_name(exemplar.xid_to_record(AB, xid, u))
    bname = get_ok_name(exemplar.xid_to_record(AB, xid, v))
    if aname == bname:
      return aname
    else:
      return "%s[%s]" % (aname, bname)
  return ";".join(map(makename, xids))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
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
  assert a_path != b_path
  d_path = '-'
  with rows.open(a_path) as a_rows:
    with rows.open(b_path) as b_rows:
      AB = ingest_workspace(a_rows.rows(), b_rows.rows(),
                            A_name=a_name, B_name=b_name)
      theory.theorize(AB, False)
      with rows.open(d_path, "w") as d_gen:
        d_gen.write_rows(generate_risk_report(AB))
