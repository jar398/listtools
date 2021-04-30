#!/usr/bin/env python3

# Write a separate file listing all names, and remove non-accepted
# name records (i.e. those for synonyms) from the taxon file

import sys, csv, argparse

# Compare map.py
accepted_taxon_id_label = "acceptedNameUsageID"
accepted_page_id_label = "acceptedEOLid"

def divert_names(inport, outport, synport):
  reader = csv.reader(inport)
  header = next(reader)
  id_label = "EOLid"
  if windex(header, id_label) == None: id_label = "taxonID"
  id_pos = windex(header, id_label)
  if id_pos == None:
    print("** names: No id column (%s) found in input header" % id_label,
          file=sys.stderr)
    print("Header: %s" % header, file=sys.stderr)
    assert False
  aid_label = accepted_page_id_label
  if aid_label == None: aid_label = accepted_taxon_id_label
  accepted_id_pos = windex(header, aid_label)
  if accepted_id_pos == None:
    print("** names: No accepted id column (%s) found in input header" % aid_label,
          file=sys.stderr)
    print("Header: %s" % header, file=sys.stderr)
    assert False
  columns_to_keep = list(range(len(header)))
  del columns_to_keep[accepted_id_pos]
  # Also delete?  taxonomic and nomenclatural status

  writer = csv.writer(outport)
  writer.writerow([header[i] for i in columns_to_keep])
  count = 0
  names_writer = csv.writer(synport)
  names_writer.writerow(header)
  syn_count = 0
  for row in reader:
    is_accepted = False
    id = row[id_pos]
    accepted_id = row[accepted_id_pos]
    if accepted_id == id:
      is_accepted = True
    elif accepted_id == None or accepted_id == '':    # can't remember
      is_accepted = True
      row[accepted_id_pos] = id
    if is_accepted:
      writer.writerow([row[i] for i in columns_to_keep])
      count += 1
    else:
      syn_count += 1
    names_writer.writerow(row)
  print("names: %s accepted, %s synonyms" % (count, syn_count),
        file=sys.stderr)

def windex(header, fieldname):
  if fieldname in header:
    return header.index(fieldname)
  else:
    return None

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Remove synonyms from input stream, while also creating
    a separate file containing all names.
    """)
  parser.add_argument('--names', 
                      help='path to file where names will be stored')
  args=parser.parse_args()
  with open(args.names, "w") as synport:
    divert_names(sys.stdin, sys.stdout, synport)
