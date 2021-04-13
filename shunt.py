#!/usr/bin/env python3

# Shunt records that have a non-null acceptedNameUsageID field
# (i.e. those for synonyms) into a separate file

import sys, csv, argparse

# Compare map.py
accepted_page_id_label = "acceptedEOLid"

def shunt_synonyms(inport, outport, synport):
  reader = csv.reader(inport)
  header = next(reader)
  taxon_id_pos = windex(header, "taxonID")
  accepted_taxon_id_pos = windex(header, "acceptedNameUsageID")
  accepted_page_id_pos = windex(header, accepted_page_id_label)
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
    # DH 0.9 has 10 extra null columns at the end of synonym lines
    if len(row) > len(header):
      row = row[0:len(header)]
    if len(row) != len(header):
      print(("** shunt: Unexpected number of columns: have %s want %s" %
             (len(row), len(header))),
            file=sys.stderr)
      print(("** shunt: Row is %s" % (row,)), file=sys.stderr)
      assert False

    if ((accepted_taxon_id_pos and (row[accepted_taxon_id_pos] and
                                    row[accepted_taxon_id_pos] != row[taxon_id_pos])) or
        (accepted_page_id_pos and row[accepted_page_id_pos])):
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
