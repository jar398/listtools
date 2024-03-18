#!/usr/bin/env python3

"""
 Converts an NCBI dump to a DwC TNU (taxon) file.
 Fold the scientific name (-> canonicalName) and authority (if it
 extends the scientific name) into the taxon record.

 python3 ncbi_to_dwc.py from >to
   'from' is a directory containing .dmp files from ncbi
   DwC file is written to standard output (in CSV format with header row)

 E.g.
  python3 ncbi_to_dwc.py ~/ncbi/2015-01-01/dump > ~/ncbi/2015-01-01.csv

 Copied from ../cldiff/src/ncbi_to_dwc.py on 2021-08-24
"""

import sys, os, csv, argparse, regex as re

from util import stable_hash

def ncbi_to_dwc(indir, outfile):
  assert os.path.exists(indir)
  accepteds = read_accepteds(os.path.join(indir, "nodes.dmp"))
  names = read_names(os.path.join(indir, "names.dmp"))
  merged = read_merged(os.path.join(indir, "merged.dmp"))
  (synonyms, scinames, authorities) = collate_names(names, accepteds)
  emit_dwc(accepteds, synonyms, scinames, authorities, merged, outfile)

def write_row(writer,
              taxonID, ncbi_id, parentNameUsageID, taxonRank,
              acceptedNameUsageID, scientificName, canonicalName,
              taxonomicStatus, nomenclaturalStatus):
  writer.writerow([taxonID, ncbi_id, parentNameUsageID, taxonRank,
                   acceptedNameUsageID, scientificName, canonicalName,
                   taxonomicStatus, nomenclaturalStatus])

def emit_dwc(accepteds, synonyms, scinames, authorities, merged, outfile):
  writer = csv.writer(outfile)
  write_row(writer,
            "taxonID", "managed_id", "parentNameUsageID", "taxonRank",
            "acceptedNameUsageID", "scientificName", "canonicalName",
            "taxonomicStatus", "nomenclaturalStatus")
  # Accepted names
  for (taxid, parent_id, rank) in accepteds:
    sci = scinames.get(taxid, None)
    can = authorities.get(taxid, None)
    write_row(writer,
              taxid, managed_id(taxid), parent_id, clean_rank(rank, can),
              None, sci, can,
              "accepted", None)
  # Synonyms.  taxonID is the id of the synonym; taxid of the accepted
  for (taxid, text, kind, spin) in synonyms:
    # synonym is a taxonomic status, not a nomenclatural status
    if kind == "synonym": kind = None
    taxonID = taxid + "." + spin      # primary key in this file
    rank = None
    if kind == "authority":    # shouldn't happen, but does
      can = None
      write_row(writer,
                taxonID, None, None, clean_rank(rank, can),
                taxid, text, can,
                "synonym", kind)
    else:
      can = text   # canonicalName
      write_row(writer,
                taxonID, None, None, clean_rank(rank, can),
                taxid, None, can,
                "synonym", kind)
  for (old_id, new_id) in merged:
    can = "deprecated taxon that had taxid %s" % old_id
    write_row(writer,
              old_id, managed_id(old_id), None, None,
              new_id, None, can,
              "synonym", "merged id")

def managed_id(taxid):
  return "ncbi:%s" % taxid

# Input: list of (id, text, kind, spin) from names.txt file
# Output: list of (id, text, kind, spin); 
#         dict: id -> text [canonical names];
#         dict: id -> text [authorities]

def collate_names(names, accepteds):
  keep = []
  scinames = {}
  for row in names:
    (id, text, kind, _) = row
    if kind == "scientific name":
      scinames[id] = text
    else:
      keep.append(row)
  # Remove the canonical names, keep the rest
  names2 = keep
  keep = None #GC
  print ("# %s canonicalNames (NCBI calls them 'scientific names')" %
         len(scinames),
         file=sys.stderr)
  synonyms = []
  authorities = {}
  for row in names2:
    (id, text, kind, _) = row
    if kind == "authority":
      probe = scinames.get(id, None)
      if probe and text.startswith(probe):
        authorities[id] = text
      else:
        synonyms.append(row)
    else:
      synonyms.append(row)
  # Remove authorities (= scientific names), keep the rest
  print ("# %s scientificNames (NCBI calls them 'authorities')" %
         len(authorities),
         file=sys.stderr)
  return (synonyms, scinames, authorities)

def read_accepteds(nodes_path):
  accepteds = []
  # Read the nodes file
  with open(nodes_path, "r") as infile:
    for row in csv.reader(infile,
                          delimiter="\t",
                          quotechar="\a",
                          quoting=csv.QUOTE_NONE):
      # tax_id, |, parent tax_id, |, rank, ... other stuff we don't use ...
      rank = row[4]
      if rank == "clade" or rank == "no rank":
        rank = None
      accepteds.append((row[0], row[2], rank)) # id, superior id, rank
  print ("# %s accepted names" % len(accepteds), file=sys.stderr)
  return accepteds

def read_names(names_path):
  names = []
  # Read the names file
  with open(names_path, "r") as infile:
    # Depends on names being grouped by taxa
    previous_id = None
    spins = {}
    for row in csv.reader(infile,
                          delimiter="\t",
                          quotechar="\a",
                          quoting=csv.QUOTE_NONE):
      # 0	tax_id					-- the id of node associated with this name
	    # 2 name_txt				-- name itself
	    # 4 unique name  		-- the unique variant of this name if name not unique
	    # 6 name class			-- (synonym, common name, ...)
      id = row[0]
      if str(id) != previous_id:
        previous_id = id
        spins = {}
      name = row[2]
      unique = row[4]
      nom_status = row[6]
      spin = stable_hash((name, nom_status, unique))
      out = (id, name, nom_status, spin)
      if not spin in spins:
        names.append(out)
        spins[spin] = True
      else:
        # 263 of these in 2020-08-01
        # print("# Ignoring redundant name record %s" % (out,), file=sys.stderr)
        pass
  print ("# %s names" % len(names), file=sys.stderr)
  return names

def read_merged(merged_path):
  merged = []
  # Read the merged file
  with open(merged_path, "r") as infile:
    for row in csv.reader(infile,
                          delimiter="\t",
                          quotechar="\a",
                          quoting=csv.QUOTE_NONE):
      # old_tax_id, |, new_tax_id
      merged.append((row[0], row[2]))
  print ("# %s merged taxa" % len(merged), file=sys.stderr)
  return merged

def clean_rank(rank, can):
  if not rank:

    # Lose.  This turns up too many common names.  The right thing would be to get the rank from the
    # record for the accepted name.

    if False and can and binomial_re.match(can):
      if (can.endswith('virus') or can.endswith('group') or # NCBI
          can.endswith('clade') or can.endswith('complex') or
          can.endswith('viruses') or can.endswith('isolates') or
          can.endswith('plasma') or can.endswith('phage')):
        pass
      else:
        print("# Assigning species rank to %s" % can,
              file=sys.stderr)
        return 'species'
  return rank

binomial_re = re.compile('\p{Uppercase_Letter}\p{Lowercase_Letter}+ \p{Lowercase_Letter}{2,}$')

# When invoked from command line:

if __name__ == '__main__':
  parser = argparse.ArgumentParser(prog='ncbi_to_dwc',
                                   description="Writes a CSV file in Darwin Core form to standard output")
  parser.add_argument('taxdump', help='directory containing NCBI dump i.e. names.txt, merged.txt, and so on')
  args = parser.parse_args()
  ncbi_to_dwc(args.taxdump, sys.stdout)
