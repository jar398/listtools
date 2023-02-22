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
import util

from util import MISSING, windex, csv_parameters

def extract_subset(infile, hier_file, root_id, outfile):
  (topo, root_tid) = read_topology(hier_file, root_id)
  assert root_tid
  all = closure(topo, root_tid)
  write_subset(infile, root_tid, all, topo, outfile)

def write_subset(infile, root_id, all, topo, outfile):
  reader = csv.reader(infile)
  head = next(reader)

  tid_column = head.index("taxonID") 
  aid_column = head.index("acceptedNameUsageID")
  pid_column = head.index("parentNameUsageID")

  writer = csv.writer(outfile)
  writer.writerow(head)
  for row in reader:
    tid = row[tid_column]
    if tid in all:
      aid = row[aid_column]
      if not aid in all:
        row[aid_column] = MISSING
      pid = row[pid_column]
      if not pid in all:
        row[pid_column] = MISSING
      writer.writerow(row)

# Transitive closure of accepted records

def closure(topo, root_id):
  print("-- subset: transitive closure starting with %s" % root_id, file=sys.stderr)
  (children, _) = topo[root_id]
  print("-- subset: root has %s children" % len(children), file=sys.stderr)
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
  print("-- subset: %s items in transitive closure" % len(all), file=sys.stderr)
  return all

def read_topology(hier_file, root_id):
  # Keyed by taxon id
  topo = {}
  # (delimiter, quotechar, mode) = csv_parameters(hier_path)
  counter = 0
  root_tid = None
  print("# subset: scanning to obtain hierarchy",
        flush=True,
        file=sys.stderr)
  reader = csv.reader(hier_file)  #, delimiter=delimiter, quotechar=quotechar, quoting=mode
  head = next(reader)

  tid_column = windex(head, "taxonID") 
  pid_column = windex(head, "parentNameUsageID")
  aid_column = windex(head, "acceptedNameUsageID")
  sid_column = windex(head, "taxonomicStatus")
  name_column = windex(head, "canonicalName")
  sci_column = windex(head, "scientificName")

  if tid_column == None:      # usually 0
    print("** No taxonID column found", file=sys.stderr)
  if pid_column == None:
    print("** No parentNameUsageID column found", file=sys.stderr)
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
    if accepted_id != MISSING:
      assert accepted_id != tid
      (_, syns) = get_topo_record(accepted_id, topo)
      syns.append(tid)
    if parent_id != MISSING and parent_id != tid:
      (children, _) = get_topo_record(parent_id, topo)
      children.append(tid)
    if (tid == root_id or
        (name_column != None and row[name_column] == root_id) or
        (sci_column != None and row[sci_column] == root_id)):
      root_tid = tid

  print("-- subset: %s hierarchy items of which %s have children and/or synonyms" %
        (counter, len(topo)), file=sys.stderr)

  if root_tid == None:
    print("*** Did not find root taxonID or canonicalName %s" % root_id)
  return (topo, root_tid or root_id)

def get_topo_record(tid, topo):
  record = topo.get(tid)
  if not record:
    record = ([], [])
    topo[tid] = record
  return record

# extract_subset(checklist, taxonomy, root_id, outfile)

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--hierarchy',
                      help="file from which to extract complete hierarchy")
  parser.add_argument('--root',
                      help="taxonID or name of root of subtree to be extracted")
  args = parser.parse_args()
  with open(args.hierarchy, 'r') as hier_file:
    extract_subset(sys.stdin, hier_file, args.root, sys.stdout)
