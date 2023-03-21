#!/usr/bin/env python3

# This should always be the first step in a processing pipeline.
#  - Normalizes tsv to csv
#  - Makes sure all rows have the same number of fields

import sys, csv, regex, hashlib, argparse
from util import csv_parameters, windex, stable_hash, log, MISSING

def start_csv(inport, params, outport, args):
  cleanp = args.clean

  (d, q, g) = params
  reader = csv.reader(inport, delimiter=d, quotechar=q, quoting=g)
  in_header = next(reader)
  # Use built-in python csv sniffer ??
  if len(in_header) == 1:
    if "," in in_header[0] or "\t" in in_header[0]:
      log("** start: Suspicious in_header")
      log("** start: in_header is %s" % (in_header,))

  def fix(h):
    if h[0] == '\ufeff': h = h[1:]
    if h.startswith('dwc:'): h = h[4:]    # CoL 2021
    if h == 'taxonId': return 'taxonID'   # artsdatabanken 2023
    else: return h
  in_header = [fix(h) for h in in_header]

  pk_col = args.pk or 'taxonID'
  pk_pos_in = windex(in_header, pk_col)

  can_pos = windex(in_header, "canonicalName")
  sci_pos = windex(in_header, "scientificName")
  rank_pos = windex(in_header, "taxonRank")
  source_pos = windex(in_header, "source")
  landmark_pos = windex(in_header, "Landmark")
  if landmark_pos != None: in_header[landmark_pos] = "landmark_status"
  accepted_pos = windex(in_header, "acceptedNameUsageID")
  parent_pos = windex(in_header, "parentNameUsageID")
  tax_status_pos = windex(in_header, "taxonomicStatus")
  taxon_id_pos = windex(in_header, "taxonID")

  if taxon_id_pos == None:
    log("** No taxonID column. Header is %s" % (in_header,))

  out_header = in_header

  must_affix_pk = (pk_pos_in == None)
  if must_affix_pk:
    log("-- Appending a %s column" % pk_col)
    out_header = in_header + [pk_col]

  # --managed taxonID --prefix GBIF:   must always be used together
  if args.managed:
    assert not must_affix_pk
    (prefix, managed_col) = args.managed.split(':')     # GBIF:taxonID

    # Position of source id (usually taxonID) to be copied to managed_id column
    managed_col_pos = windex(in_header, managed_col)
    assert managed_col_pos != None
    assert windex(in_header, "managed_id") == None

    log("-- Appending a managed_id column")
    out_header = out_header + ["managed_id"]

  # log("# Output header: %s" % (out_header,))

  # Do these after header modifications
  pk_pos_out = windex(out_header, pk_col)

  writer = csv.writer(outport) # CSV not TSV
  writer.writerow(out_header)
  count = 0
  trimmed = 0
  names_cleaned = 0
  ranks_cleaned = 0
  accepteds_normalized = 0
  minted = 0
  seen_pks = {}
  previous_pk = 0
  conflicts = 0
  senior = 0
  managed = 0
  for row in reader:

    # Deal with raggedness if any
    if len(row) > len(in_header):
      row = row[0:len(in_header)]
      trimmed += 1
    elif len(row) < len(in_header):
      log(("** start: Unexpected number of columns: have %s want %s" %
             (len(row), len(in_header))))
      log(("** start: Row is %s" % (row,)))
      assert False

    # Filter out senior synonyms
    if tax_status_pos != None and row[tax_status_pos] == "senior synonym":
      senior += 1
      continue

    # Clean up if wrong values in canonical and/or scientific name columns
    if cleanp:
      if clean_name(row, can_pos, sci_pos):
        names_cleaned += 1
      if clean_rank(row, rank_pos, can_pos):
        ranks_cleaned += 1

    if normalize_accepted(row, taxon_id_pos, parent_pos, accepted_pos):
      accepteds_normalized += 1

    # Shouldn't have both accepted and parent
    if False:
      if accepted_pos and parent_pos and row[accepted_pos] and row[parent_pos]:
        row[parent_pos] = MISSING

    stat = row[tax_status_pos]
    indication_2 = (stat.startswith("accepted") or
                    stat.startswith("valid") or
                    stat.startswith("dubious"))

    # Two ways to test whether a usage is accepted/dubious
    usage_id = row[pk_pos_in] if pk_pos_in != None else MISSING
    au = row[accepted_pos] if accepted_pos != None else MISSING
    indication_1 = (au == MISSING or au == usage_id)

    if indication_1 != indication_2 and conflicts < 10:
      log("-- %s has accepted %s but taxstatus %s" %
             (usage_id, au, stat))
    conflicts += 1

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

    # Extend for new primary key column if necessary
    if must_affix_pk:
      out_row = row + [MISSING]
    else:
      out_row = row
    pk = out_row[pk_pos_out]

    # Add managed_id if necessary
    if args.managed:
      if prefix == "mdd" and (('synonym' in stat) or
                              row[rank_pos] != 'species'):
        # id added by Prashant, flush
        managed_id = MISSING
      elif pk:
        id = row[managed_col_pos]
        managed_id = "%s:%s" % (prefix, id) if id else MISSING
        managed += 1
      else:
        managed_id = MISSING
      out_row = out_row + [managed_id]

    # Set primary key if duplicate or mising
    if pk in seen_pks:
      log("** %s is not a good primary key column.  Two or more rows with %s = %s\n" %
            (pk_col, pk_col, pk))
      pk = fresh_pk(row, out_header)
      out_row[pk_pos_out] = pk
    elif pk == MISSING:
      pk = fresh_pk(row, out_header)
      minted += 1
      out_row[pk_pos_out] = pk
    assert pk != MISSING
    seen_pks[pk] = True

    assert len(out_row) == len(out_header)
    writer.writerow(out_row)
    count += 1
    if count % 250000 == 0:
      log("row %s id %s" % (count, row[taxon_id_pos]))
  log("-- start: %s rows, %s columns, %s ids minted, %s accepteds normalized" %
        (count, len(in_header), minted, accepteds_normalized))
  log("-- start: %s names cleaned, %s ranks cleaned" %
        (names_cleaned, ranks_cleaned))
  if senior > 0:
    log("-- start: filtered out %s senior synonyms" % senior)
  if managed > 0:
    log("-- start: %s managed ids" % managed)
  if trimmed > 0:
    # Ignoring extra values is appropriate behavior for DH 0.9.  But
    # elsewhere we might want ragged input to be treated as an error.
    log("-- start: trimmed extra values from %s rows" % (trimmed,))
    
def fresh_pk(row, out_header):
  try:
    return stable_hash("^".join(row))
  except:
    assert False, row

def normalize_accepted(row, taxon_id_pos, parent_pos, accepted_pos):
  if accepted_pos != None and row[accepted_pos] != MISSING:
    if (taxon_id_pos != None and
        row[accepted_pos] == row[taxon_id_pos]):
      row[accepted_pos] = MISSING
      return True
    elif False and parent_pos != None and row[accepted_pos] != row[parent_pos]:
      # Norwegian national checklist has bogus parent pointers for synonyms
      row[parent_pos] = MISSING
      return True
  return False

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

# Returns True if a change was made

def clean_name(row, can_pos, sci_pos):
  mod = False
  if can_pos != None and sci_pos != None:
    c = row[can_pos]
    s = row[sci_pos]
    if s == MISSING:
      if is_scientific(c):
        row[sci_pos] = c
        row[can_pos] = MISSING    # ?
        mod = True
    elif is_scientific(s):
      pass
    # s is nonnull and not 'scientific'
    elif c == MISSING:
      # swap
      row[sci_pos] = MISSING
      row[can_pos] = s
      # log("start: c := s") - frequent in DH 1.1
      mod = True
    elif c == s:
      if is_scientific(c):
        # should gnparse!
        row[can_pos] = MISSING
      else:
        row[sci_pos] = MISSING
      # log("start: flush s") - frequent in DH 1.1
      mod = True
  if sci_pos != None and row[sci_pos]:
    repl = row[sci_pos].replace(' and ', ' & ') # frequent in DH 1.1
    repl = repl.replace(',,', ',')    # kludge for MDD 1.0
    if repl != row[sci_pos]:
      row[sci_pos] = repl
      mod = True
  return mod

has_auth_re = regex.compile(u" (\()\p{Uppercase_Letter}[\p{Letter}-]+")

def is_scientific(name):
  return has_auth_re.search(name)

def clean_rank(row, rank_pos, can_pos):
  if rank_pos != None and row[rank_pos] == MISSING and \
     can_pos and binomial_re.match(row[can_pos]):
    row[rank_pos] = 'species'
    return True
  return False

binomial_re = regex.compile("[A-Z][a-z]+ [a-z]{2,}")


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Normalize to csv and check that number of columns is same in every row.
    CSV rows are written to standard output.
    """)
  parser.add_argument('--input', default=None,
                      help="""name of input file.  TSV assumed unless
                              name contains ".csv".  Default to CSV from stdin""")
  parser.add_argument('--pk', default=None,
                      help='name of column containing primary key')
  parser.add_argument('--clean', dest='clean', action='store_true',
                      help='clean up scientificName and canonicalName a little bit')
  parser.add_argument('--no-clean', dest='clean', action='store_false')
  parser.add_argument('--managed',
                      help='prefix and source column for managed record ids, e.g. eol:EOLid')
  parser.set_defaults(clean=True)

  args=parser.parse_args()
  inpath = args.input
  if inpath == None:
    params = csv_parameters("foo.csv")
    start_csv(sys.stdin, params, sys.stdout, args)
  else:
    params = csv_parameters(args.input)
    with open(args.input, "r") as inport:
      start_csv(inport, params, sys.stdout, args)

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
