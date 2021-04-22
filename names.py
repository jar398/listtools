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
  taxon_id_pos = windex(header, "taxonID")
  page_id_pos  = windex(header, "EOLid")
  accepted_taxon_id_pos = windex(header, accepted_taxon_id_label)
  accepted_page_id_pos  = windex(header, accepted_page_id_label)
  if (accepted_taxon_id_pos == None and accepted_page_id_pos == None):
    print("** names: No accepted-id column found in input header",
          file=sys.stderr)
    assert False
  columns_to_keep = [i
                     for i in list(range(len(header)))
                     if (i != accepted_taxon_id_pos and i != accepted_page_id_pos)]

  writer = csv.writer(outport)
  writer.writerow([header[i] for i in columns_to_keep])
  count = 0
  names_writer = csv.writer(synport)
  names_writer.writerow(header)
  syn_count = 0
  for row in reader:
    # Separate synonyms from accepteds
    is_accepted = False
    if accepted_page_id_pos != None:
      accepted_page_id = row[accepted_page_id_pos]
      page_id = row[page_id_pos] if (page_id_pos != None) else None
      if accepted_page_id == page_id:
        is_accepted = True
      elif accepted_page_id == None:
        is_accepted = True
        row[accepted_page_id_pos] = page_id
    elif accepted_taxon_id_pos != None:
      accepted_taxon_id = row[accepted_taxon_id_pos]
      taxon_id = row[taxon_id_pos] if (taxon_id_pos != None) else None
      if accepted_taxon_id == taxon_id:
        is_accepted = True
      elif accepted_taxon_id == None:
        is_accepted = True
        row[accepted_taxon_id_pos] = taxon_id
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
    Move rows for taxonomic synonyms to a separate file, leaving
    others unchanged.
    CSV rows are read from standard input and non-synonym 
    rows written to standard output.
    """)
  parser.add_argument('--names', 
                      help='name of file where names will be stored')
  args=parser.parse_args()
  with open(args.synonyms, "w") as synport:
    siphon_names(sys.stdin, sys.stdout, synport)
