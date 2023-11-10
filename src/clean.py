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
      log("** clean: Suspicious in_header")
      log("** clean: in_header is %s" % (in_header,))

  def fix(h):
    if h[0] == '\ufeff': h = h[1:]
    if h.startswith('dwc:'): h = h[4:]    # CoL 2021
    if h == 'taxonId': return 'taxonID'   # artsdatabanken 2023
    else: return h
  in_header = [fix(h) for h in in_header]

  pk_col = args.pk or 'taxonID'
  pk_pos_in = windex(in_header, pk_col)
  assert pk_pos_in != None

  can_pos = windex(in_header, "canonicalName")
  sci_pos = windex(in_header, "scientificName")
  rank_pos = windex(in_header, "taxonRank")
  source_pos = windex(in_header, "source")
  landmark_pos = windex(in_header, "Landmark")
  if landmark_pos != None: in_header[landmark_pos] = "landmark_status"
  accepted_pos = windex(in_header, "acceptedNameUsageID")
  parent_pos = windex(in_header, "parentNameUsageID")
  tax_status_pos = windex(in_header, "taxonomicStatus")
  auth_pos = windex(in_header, "Authorship")
  taxon_id_pos = windex(in_header, "taxonID")

  if taxon_id_pos == None:
    log("** No taxonID column. Header is %s" % (in_header,))

  out_header = in_header

  if can_pos == None:
    can_pos = len(out_header)
    out_header = out_header + ["canonicalName"]
    log("** Adding canonicalName column")
  if sci_pos == None:
    sci_pos = len(out_header)
    out_header = out_header + ["scientificName"]
    log("** Adding scienticicName column")

  # --managed taxonID --prefix GBIF:   must always be used together
  if args.managed:
    (prefix, managed_col) = args.managed.split(':')     # GBIF:taxonID

    # Position of source id (usually taxonID) to be copied to managed_id column
    managed_col_pos = windex(in_header, managed_col)
    assert managed_col_pos != None
    managed_id_pos = windex(in_header, "managed_id")
    if managed_id_pos == None:
      log("-- Appending a managed_id column")
      out_header = out_header + ["managed_id"]
      managed_id_pos = windex(out_header, "managed_id")

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
  # In COL 2023 we get:
  # _csv.Error: field larger than field limit (131072)
  csv.field_size_limit(131072 * 4)
  for in_row in reader:

    for field in in_row:
      if len(field) > 10000:
        log("# Really long row: %s %s" % (count, in_row[0]))

    # Deal with raggedness if any
    if len(in_row) > len(in_header):
      in_row = in_row[0:len(in_header)]
      trimmed += 1
    elif len(in_row) < len(in_header):
      log(("** clean: Unexpected number of columns: have %s want %s" %
             (len(in_row), len(in_header))))
      log(("** clean: in_row is %s" % (in_row,)))
      log(("** clean: in_header is %s" % (in_header,)))
      assert False

    # Filter out senior synonyms
    if tax_status_pos != None and in_row[tax_status_pos] == "senior synonym":
      senior += 1
      continue

    if normalize_accepted(in_row, taxon_id_pos, parent_pos, accepted_pos):
      accepteds_normalized += 1

    out_row = in_row + [MISSING] * (len(out_header) - len(in_row))

    # Shouldn't have both accepted and parent
    if False:                   # why disabled?
      if accepted_pos and parent_pos and in_row[accepted_pos] and in_row[parent_pos]:
        out_row[parent_pos] = MISSING

    stat = in_row[tax_status_pos]
    indication_2 = (stat.startswith("accepted") or
                    stat.startswith("valid") or
                    stat.startswith("dubious") or
                    stat == "provisionally accepted")  #seen in GBIF

    # Two ways to test whether a usage is accepted/dubious
    pk = in_row[pk_pos_in]
    au = in_row[accepted_pos] if accepted_pos != None else MISSING
    indication_1 = (au == MISSING or au == pk)

    if indication_1 != indication_2 and conflicts < 10:
      log("-- %s has accepted '%s' but taxstatus '%s'" %
          (pk, au, stat))
    conflicts += 1

    # landmark_status is specific to EOL
    if landmark_pos != None: 
      l = in_row[landmark_pos]
      if l != MISSING:
        e = int(l)
        # enum landmark: %i[no_landmark minimal abbreviated extended full]
        if   e == 1: out_row[landmark_pos] = 'minimal'
        elif e == 2: out_row[landmark_pos] = 'abbreviated'
        elif e == 3: out_row[landmark_pos] = 'extended'
        elif e == 4: out_row[landmark_pos] = 'full'
        else: out_row[landmark_pos] = MISSING

    # If multiple sources (smasher output), use only the first
    if source_pos != None and in_row[source_pos] != MISSING:
      out_row[source_pos] = in_row[source_pos].split(',', 1)[0]

    # Clean up if wrong values in canonical and/or scientific name columns
    if cleanp:
      if clean_name(out_row, can_pos, sci_pos):
        names_cleaned += 1
      if clean_rank(out_row, rank_pos, can_pos):
        ranks_cleaned += 1
      if auth_pos != None:      # clean_auth ...
        a = in_row[auth_pos]
        a = a.strip()
        if a.endswith(').'): a = a[0:-1] # for MDD
        out_row[auth_pos] = a

    # Add primary key if duplicate(?) or missing
    if pk == MISSING:
      pk = fresh_pk(in_row, out_header)
      minted += 1
      out_row[pk_pos_out] = pk
    elif pk in seen_pks:
      log("** %s Two or more rows have %s = %s\n" %
            (pk_col, pk_col, pk))
      pk = fresh_pk(in_row, out_header)
      out_row[pk_pos_out] = pk
    assert pk != MISSING
    seen_pks[pk] = True

    # Add managed_id if necessary
    if args.managed:            # --managed prefix:column
      if managed_col_pos < len(in_row):
        id = in_row[managed_col_pos]
        managed_id = "%s:%s" % (prefix, id) if id else MISSING
        managed += 1
      else:
        managed_id = MISSING
      out_row[managed_id_pos] = managed_id

    assert len(out_row) == len(out_header)
    writer.writerow(out_row)
    count += 1
    if count % 250000 == 0:
      log("row %s id %s" % (count, in_row[taxon_id_pos]))
  log("-- clean: %s rows, %s columns, %s ids minted, %s accepteds normalized" %
        (count, len(in_header), minted, accepteds_normalized))
  log("-- clean: %s names cleaned, %s ranks cleaned" %
        (names_cleaned, ranks_cleaned))
  if senior > 0:
    log("-- clean: filtered out %s senior synonyms" % senior)
  if managed > 0:
    log("-- clean: %s managed ids" % managed)
  if trimmed > 0:
    # Ignoring extra values is appropriate behavior for DH 0.9.  But
    # elsewhere we might want ragged input to be treated as an error.
    log("-- clean: trimmed extra values from %s rows" % (trimmed,))
    
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
  c = row[can_pos].strip() if can_pos != None else MISSING
  s = row[sci_pos].strip() if sci_pos != None else MISSING
  if s == MISSING:
    if is_scientific(c):
      s = c
      c = MISSING    # ?
  elif is_scientific(s):
    pass
  # s is nonnull and not 'scientific'
  elif c == MISSING:
    # swap
    c = s
    s = MISSING
    # log("clean: c := s") - frequent in DH 1.1
  elif c == s:
    if is_scientific(c):
      c = MISSING
    else:
      s = MISSING
    # log("clean: flush s") - frequent in DH 1.1
  # Remove subgenus: Foo (Bar) -> Foo
  s = remove_subgenus(s)
  s = s.replace(' and ', ' & ') # frequent in DH 1.1 ?
  s = s.replace(',,', ',')    # kludge for MDD 1.0
  s = s.replace('  ', ' ')    # kludge for MSW 3
  if s.endswith(').'): s = s[0:-1] # for MSW 3
  if s != row[sci_pos]:
    row[sci_pos] = s
    mod = True
  if can_pos and c != row[can_pos]:
    row[can_pos] = c
    mod = True
  return mod

has_subgenus_re = regex.compile(u'(\p{Uppercase_Letter}[\p{Letter}-]+)( \(\p{Uppercase_Letter}[\p{Letter}-]+\))(.*)')

def remove_subgenus(name):
  m = has_subgenus_re.match(name)
  if m:
    return m[1] + m[3]
  else:
    return name

# This may be too liberal... insist on there being a year?
has_auth_re = regex.compile(u' (\(?)\p{Uppercase_Letter}[\p{Letter}-]+')

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
