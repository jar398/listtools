#!/usr/bin/env python3

"""
 Makes a subset of a checklist based on one subtree of a taxonmoy.

 python3 subset_dwc.py [--taxonomy tax_dwc] source_dwc id --out out_dwc

 Assumption: every accepted record has a taxonID
"""

# This file was initiated on 20 April 2021 by making a copy of
# https://github.com/jar398/cldiff/blob/master/src/subset_dwc.py

debug = False

import sys, os, csv, argparse

from util import MISSING, csv_parameters

def main(infile, hier_path, root_id, outfile):
  topo = read_topology(hier_path)
  all = closure(topo, root_id)
  write_subset(infile, root_id, all, topo, outfile)

def write_subset(infile, root_id, all, topo, outfile):
  reader = csv.reader(infile)
  head = next(reader)

  tid_column = head.index("taxonID") 
  aid_column = head.index("acceptedNameUsageID")
  pid_column = head.index("parentNameUsageID")
  sid_column = head.index("taxonomicStatus")

  writer = csv.writer(outfile)
  writer.writerow(head)
  for row in reader:
    tid = row[tid_column]
    if tid in all:
      writer.writerow(row)

# Transitive closure of accepted records

def closure(topo, root_id):
  print("Computing transitive closure starting from %s" % root_id,
        flush=True, file=sys.stderr)
  all = {}
  empty = []
  def descend(id):
    if not id in all:
      all[id] = True
      if id in topo:
        (children, synonyms) = topo[id]
        for child in children:
          descend(child)
        for syn in synonyms:
          descend(syn)
  descend(root_id)
  assert len(all) > 1
  print("  %s nodes in transitive closure" % len(all), file=sys.stderr)
  return all

def read_topology(hier_path):
  # Keyed by taxon id
  topo = {}
  (delimiter, quotechar, mode) = csv_parameters(hier_path)
  counter = 0
  with open(hier_path, "r") as infile:
    print("Scanning %s to obtain hierarchy" % hier_path, flush=True, file=sys.stderr)
    reader = csv.reader(infile, delimiter=delimiter, quotechar=quotechar, quoting=mode)
    head = next(reader)

    tid_column = head.index("taxonID") 
    pid_column = head.index("parentNameUsageID")
    aid_column = head.index("acceptedNameUsageID")
    sid_column = head.index("taxonomicStatus")

    if tid_column == None:      # usually 0
      print("** No taxonID column found", file=sys.stderr)
    if pid_column == None:
      print("** No taxonID column found", file=sys.stderr)
    if aid_column == None:
      print("** No acceptedNameUsageID column found", file=sys.stderr)
    if sid_column == None:
      print("** No taxonomicStatus column found", file=sys.stderr)

    for row in reader:
      counter += 1
      tid = row[tid_column]
      parent_id = row[pid_column]
      accepted_id = row[aid_column]
      status = row[sid_column]
      # Not clear which part of the record is authoritative when there
      # is a conflict.
      # (accepted_id and not accepted_id == tid))
      if accepted_id != MISSING and accepted_id != tid:
        (_, syns) = get_topo_record(accepted_id, topo)
        syns.append(tid)
      if parent_id != MISSING and parent_id != tid:
        (children, _) = get_topo_record(parent_id, topo)
        children.append(tid)
    print("  %s nodes of which %s have children and/or synonyms" %
          (counter, len(topo)), file=sys.stderr)

  return topo

def get_topo_record(tid, topo):
  record = topo.get(tid)
  if not record:
    record = ([], [])
    topo[tid] = record
  return record

# main(checklist, taxonomy, root_id, outfile)

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--hierarchy',
                      help="file from which to extract complete hierarchy")
  parser.add_argument('--root',
                      help="taxonID of root of subtree to be extracted")
  args = parser.parse_args()
  main(sys.stdin, args.hierarchy, args.root, sys.stdout)
