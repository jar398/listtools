#!/usr/bin/env python3

# Write a separate file listing all names, and remove non-accepted
# name records (i.e. those for synonyms) from the taxon file

import sys, csv, argparse

# Compare map.py
accepted_taxon_id_label = "acceptedNameUsageID"
accepted_page_id_label = "acceptedEOLid"

def siphon_names(inport, outport, synport):
  reader = csv.reader(inport)
  header = next(reader)
  id_pos = (windex(header, "EOLid") or
            windex(header, "taxonID"))
  if id_pos == None:
    print("** names: No acceptedcolumn found in input header",
          file=sys.stderr)
    assert False
  accepted_id_pos = (windex(header, accepted_page_id_label) or
                     windex(header, accepted_taxon_id_label))
  if accepted_id_pos == None:
    print("** names: No accepted id column found in input header",
          file=sys.stderr)
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
    elif accepted_id == None:
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
                      help='name of file where names will be stored')
  args=parser.parse_args()
  with open(args.synonyms, "w") as synport:
    siphon_names(sys.stdin, sys.stdout, synport)
