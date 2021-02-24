#!/bin/env python3

# Shunt records that have a non-null acceptedNameUsageID field
# (i.e. those for synonyms) into a separate file

import sys, csv

# Compare map.py
accepted_page_id_label = "acceptedEOLid"

def shunt_synonyms(csvp, inport, outport, synport):
  (d, q, g) = csv_parameters(csvp)
  reader = csv.reader(inport, delimiter=d, quotechar=q, quoting=g)
  header = next(reader)
  taxon_id_pos = windex(header, "taxonID")
  accepted_taxon_id_pos = windex(header, "acceptedNameUsageID")
  accepted_page_id_pos = windex(header, accepted_page_id_label)
  assert accepted_taxon_id_pos or accepted_page_id_pos

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

def csv_parameters(path):
  if ".csv" in path:
    return (",", '"', csv.QUOTE_MINIMAL)
  else:
    return ("\t", "\a", csv.QUOTE_NONE)

if __name__ == '__main__':
  synfile = sys.argv[1]
  csvp = sys.argv[2] if len(sys.argv) > 2 else "stdin.csv"
  with open(synfile, "w") as synport:
    shunt_synonyms(csvp, sys.stdin, sys.stdout, synport)
