#!/usr/bin/env python3

# Shunt records that have a non-null acceptedNameUsageID or
# acceptedEOLid field (i.e. those for synonyms) into a separate file

import sys, csv, argparse

# Compare map.py
page_id_label = "EOLid"
taxon_id_label = "taxonID"
accepted_page_id_label = "acceptedEOLid"
accepted_taxon_id_label = "acceptedNameUsageID"

def shunt_synonyms(inport, outport, synport):
  reader = csv.reader(inport)
  header = next(reader)

  # Find columns of primary key (id_pos) and accepted id (aid_pos)
  id_pos = windex(header, taxon_id_label)
  aid_pos = windex(header, accepted_taxon_id_label)
  if id_pos == None or aid_pos == None:
    id_pos = windex(header, page_id_label)
    aid_pos = windex(header, accepted_page_id_label)
    if id_pos == None or aid_pos == None:
      print("** shunt: Failed to find id and accepted id columns in input",
            file=sys.stderr)
      print("Header: %s" % header,
            file=sys.stderr)

  # Remove empty accepted-id column from accepted names ist
  columns_to_keep = list(range(len(header)))
  del columns_to_keep[aid_pos]

  acc_writer = csv.writer(outport)
  acc_writer.writerow([header[i] for i in columns_to_keep])
  acc_count = 0
  syn_writer = csv.writer(synport)
  syn_writer.writerow(header)
  syn_count = 0
  for row in reader:
    # Separate synonyms from accepteds
    is_accepted = False
    id = row[id_pos]
    accepted_id = row[aid_pos]
    if accepted_id == id:
      is_accepted = True
    elif accepted_id == None or accepted_id == '':    # can't remember
      is_accepted = True
      row[aid_pos] = id
    if is_accepted:
      acc_writer.writerow([row[i] for i in columns_to_keep])
      acc_count += 1
    else:
      syn_writer.writerow(row)
      syn_count += 1
  print("shunt: %s accepted, %s synonyms" % (acc_count, syn_count),
        file=sys.stderr)

def windex(header, fieldname):
  if fieldname in header:
    return header.index(fieldname)
  else:
    return None

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Move rows for taxonomic synonyms to a separate file, keeping
    all other rows.
    Rows for names are read from standard input and accepted (non-synonym)
    rows are written to standard output.
    """)
  parser.add_argument('--synonyms', 
                      help='name of file where synonyms will be stored')
  args=parser.parse_args()
  with open(args.synonyms, "w") as synport:
    shunt_synonyms(sys.stdin, sys.stdout, synport)
