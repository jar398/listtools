#!/usr/bin/env python3

# This should always be the first step in a processing pipeline.
#  - Normalizes tsv to csv
#  - Makes sure all rows have the same number of fields

import sys, csv, re, hashlib, argparse
from util import csv_parameters, windex

MISSING = ''

def start_csv(inport, params, outport, pk_col, cleanp):
  (d, q, g) = params
  reader = csv.reader(inport, delimiter=d, quotechar=q, quoting=g)
  in_header = next(reader)
  if len(in_header) == 1:
    if "," in in_header[0] or "\t" in in_header[0]:
      print("** start: Suspicious in_header", file=sys.stderr)
      print("** start: in_header is %s" % (in_header,), file=sys.stderr)
  pk_pos_in = windex(in_header, pk_col)

  can_pos = windex(in_header, "canonicalName")
  sci_pos = windex(in_header, "scientificName")
  source_pos = windex(in_header, "source")
  dataset_pos = windex(in_header, "datasetID")
  landmark_pos = windex(in_header, "Landmark")
  if landmark_pos != None: in_header[landmark_pos] = "landmark_status"
  taxon_id_pos = windex(in_header, "taxonID")
  accepted_pos = windex(in_header, "acceptedNameUsageID")

  must_affix_pk = (pk_col and pk_pos_in == None)
  if must_affix_pk:
    print("Prepending a %s column" % pk_col, file=sys.stderr)
    out_header = [pk_col] + in_header
  else:
    out_header = in_header
  pk_pos_out = windex(out_header, pk_col)
  print("# Output header: %s" % (out_header,), file=sys.stderr)

  writer = csv.writer(outport) # CSV not TSV
  writer.writerow(out_header)
  count = 0
  trimmed = 0
  names_cleaned = 0
  accepteds_cleaned = 0
  minted = 0
  seen_pks = {}
  previous_pk = 0
  for row in reader:

    # Deal with raggedness if any
    if len(row) > len(in_header):
      row = row[0:len(in_header)]
      trimmed += 1
    elif len(row) < len(in_header):
      print(("** start: Unexpected number of columns: have %s want %s" %
             (len(row), len(in_header))),
            file=sys.stderr)
      print(("** start: Row is %s" % (row,)), file=sys.stderr)
      assert False

    # Clean up if wrong values in canonical and/or scientific name columns
    if cleanp:
      if clean_name(row, can_pos, sci_pos):
        names_cleaned += 1
      if clean_accepted(row, accepted_pos, taxon_id_pos):
        accepteds_cleaned += 1

    # landmark_status is specific to EOL
    if landmark_pos != None: 
      l = row[landmark_pos]
      if l != MISSING:
        e = int(l)
        # enum landmark: %i[no_landmark minimal abbreviated extended full]
        if   e == 1: row[landmark_pos] = 'minimal'
        elif e == 2: row[landmark_pos] = 'abbreviated'
        elif e == 3: row[landmark_pos] = 'extended'
        elif e == 4: row[landmark_pos] = 'full'
        else: row[landmark_pos] = MISSING

    # If multiple sources (smasher output), use only the first
    if source_pos != None and row[source_pos] != MISSING:
      row[source_pos] = row[source_pos].split(',', 1)[0]

    # Mint a primary key if none provided
    if must_affix_pk:
      out_row = [MISSING] + row
    else:
      out_row = row
    pk = out_row[pk_pos_out]
    if pk == MISSING:
      text = "^".join(row)
      pk = hashlib.sha1(text.encode('utf-8')).hexdigest()[0:8]
      minted += 1
      out_row[pk_pos_out] = pk
    assert pk != MISSING
    if pk in seen_pks:
      print("%s is not a good primary key column.  Two or more rows with %s = %s\n" %
            (pk_col, pk_col, pk),
            file=sys.stderr)
      assert not (pk in seen_pks)
    seen_pks[pk] = True

    writer.writerow(out_row)
    count += 1
  print("-- start: %s rows, %s columns, %s ids minted, %s names cleaned, %s accepted cleaned" %
        (count, len(in_header), minted, names_cleaned, accepteds_cleaned),
        file=sys.stderr)
  if trimmed > 0:
    # Ignoring extra values is appropriate behavior for DH 0.9.  But
    # elsewhere we might want ragged input to be treated as an error.
    print("start: trimmed extra values from %s rows" % (trimmed,),
          file=sys.stderr)
    
"""
Let c = canonicalName from csv, s = scientificName from csv,
sci = satisfies scientific name regex.
Case analysis:
  c        s
  empty    empty     Leave.
  sci      empty     Maybe copy c to s ?
  not-sci  empty     Leave.
  empty    sci       Leave; but should use gnparse.
  sci      sci       Leave; but should use gnparse.
  not-sci  sci       Leave.
  empty    not-sci   Swap.
  sci      not-sci   Swap if s is a prefix of c, otherwise leave.
  not-sci  not-sci   Remove s if =.  Otherwise leave.
"""

def clean_accepted(row, accepted_pos, taxon_id_pos):
  if accepted_pos != None and row[accepted_pos] == MISSING:
    row[accepted_pos] = row[taxon_id_pos]
    return True
  return False

# Returns True if a change was made

def clean_name(row, can_pos, sci_pos):
  if can_pos != None and sci_pos != None:
    c = row[can_pos]
    s = row[sci_pos]
    if s == MISSING:
      if is_scientific(c):
        row[sci_pos] = c
        row[can_pos] = MISSING
        return True
      return False
    if is_scientific(s):
      return False
    # s is nonnull and not 'scientific'
    if c == MISSING:
      # swap
      row[sci_pos] = None
      row[can_pos] = s
      # print("start: c := s", file=sys.stderr) - frequent in DH 1.1
      return True
    if c == s:
      if is_scientific(c):
        # should gnparse!
        row[can_pos] = None
      else:
        row[sci_pos] = None
      # print("start: flush s", file=sys.stderr) - happens all the time in 1.1
      return True
  return False

sci_re = re.compile(" [1-2][0-9]{3}\\b")

def is_scientific(name):
  return sci_re.search(name)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Normalize to csv and check that number of columns is same in every row.
    CSV rows are written to standard output.
    """)
  parser.add_argument('--input', default=None,
                      help='name of input file.  TSV assumed unless name contains ".csv"')
  parser.add_argument('--pk', default=None,
                      help='name of column containing primary key')
  parser.add_argument('--clean', dest='clean', action='store_true',
                      help='clean up scientificName and canonicalName a little bit')
  parser.add_argument('--no-clean', dest='clean', action='store_false')
  parser.set_defaults(clean=True)

  args=parser.parse_args()
  inpath = args.input
  params = csv_parameters(inpath)
  if inpath == None:
    start_csv(sys.stdin, params, sys.stdout, args.pk, args.clean)
  else:
    with open(args.input, "r") as inport:
      start_csv(inport, params, sys.stdout, args.pk, args.clean)

"""
      # Assign ids (primary keys) to any nodes that don't have them
      pk = None
      if pk_pos != None and row[pk_pos] != MISSING:
        pk = row[pk_pos]
      if pk == None and taxon_pk_pos != None and row[taxon_pk_pos] != MISSING:
        pk = row[taxon_pk_pos]
      if pk == None:
        pk = count
      if pk in seen_pks:
        spin = 1
        while True:
          dodge = "%s..%s" % (pk, spin)
          if not dodge in seen_pks:
            pk = dodge
            break

      if pk_pos == None:
        row = row + [pk]
      else:
        row[out_pk_pos] = pk

    # Every table needs to have a column of unique primary keys
    pk_pos = windex(header, pk_col)
    if pk_pos != None:
      out_pk_pos = pk_pos
    else:
      out_header = header + [pk_col]    # Add primary key
      out_pk_pos = len(header)

"""
