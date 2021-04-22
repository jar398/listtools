#!/usr/bin/env python3

# Shunt records that have a non-null acceptedNameUsageID field
# (i.e. those for synonyms) into a separate file

import sys, csv, argparse

# Compare map.py
accepted_taxon_id_label = "acceptedNameUsageID"
accepted_page_id_label = "acceptedEOLid"

def shunt_synonyms(inport, outport, synport):
  reader = csv.reader(inport)
  header = next(reader)
  taxon_id_pos = windex(header, "taxonID")
  page_id_pos  = windex(header, "EOLid")
  accepted_taxon_id_pos = windex(header, accepted_taxon_id_label)
  accepted_page_id_pos  = windex(header, accepted_page_id_label)
  if (accepted_taxon_id_pos == None and accepted_page_id_pos == None):
    print("** No accepted-id column found in input header")
    assert False

  writer = csv.writer(outport)
  writer.writerow(header)
  count = 0
  syn_writer = csv.writer(synport)
  syn_writer.writerow(header)
  syn_count = 0
  for row in reader:
    # Separate synonyms from accepteds
    if accepted_page_id_pos != None:
      a = row[accepted_page_id_pos]
      syno = (a and (page_id_pos != None) and a != row[page_id_pos])
    elif accepted_taxon_id_pos != None:
      a = row[accepted_taxon_id_pos]
      syno = (a and (taxon_id_pos != None) and a != row[taxon_id_pos])
    if syno:
      syn_writer.writerow(row)
      syn_count += 1
    else:
      writer.writerow(row)
      count += 1
  print("shunt: %s non-synonyms, %s synonyms" % (count, syn_count),
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
  parser.add_argument('--synonyms', 
                      help='name of file where synonyms will be stored')
  args=parser.parse_args()
  with open(args.synonyms, "w") as synport:
    shunt_synonyms(sys.stdin, sys.stdout, synport)
