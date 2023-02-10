#!/usr/bin/env python3

# Make use of the output of gnparse, by adding new columns to the checklist

import sys, csv, argparse
import regex
import util
from util import windex, MISSING

# was auth_re = re.compile("[A-Z][A-Za-z'-]+")
auth_re = regex.compile(u"\p{Uppercase_Letter}[\p{Letter}-]+")
year_re = regex.compile(' ([12][0-9]{3})\)?$')

def use_parse(gn_iter, check_iter):

  # gnparse output: gn_iter / gn_header
  gn_header = next(gn_iter)
  # Id,Verbatim,Cardinality,CanonicalStem,CanonicalSimple,CanonicalFull,Authorship,Year,Quality
  if len(gn_header) != 9:
    print("** Unexpected number of columns in %s" % (gn_header,),
          file=sys.stderr)
    assert False
  verbatim_pos = windex(gn_header, "Verbatim")
  cardinality_pos = windex(gn_header, "Cardinality")
  stem_pos = windex(gn_header, "CanonicalStem")
  canonical_full_pos = windex(gn_header, "CanonicalFull")
  auth_pos = windex(gn_header, "Authorship")
  year_pos = windex(gn_header, "Year")
  quality_pos = windex(gn_header, "Quality")

  checklist_header = next(check_iter)
  # May need to consult the source record too
  add_canon = not "canonicalName" in checklist_header
  add_year = not "namePublishedInYear" in checklist_header
  status_pos = windex(checklist_header, "nomenclaturalStatus")

  out_header = checklist_header
  if add_canon:
    out_header = out_header + ["canonicalName"]
  if add_year:
    out_header = out_header + ["namePublishedInYear"]
  out_header = out_header + ["canonicalStem", "type"]

  out_canonical_pos = windex(out_header, "canonicalName")
  out_year_pos = windex(out_header, "namePublishedInYear")

  row_count = 0
  trim_count = 0
  year_count = 0
  canon_count = 0

  yield out_header
  for checklist_row in check_iter:
    assert len(checklist_row) == len(checklist_header)
    row_count += 1
    gn_row = next(gn_iter)
    out_row = checklist_row
    if add_canon: out_row = out_row + [MISSING]
    if add_year:  out_row = out_row + [MISSING]

    # Fill in year if it's missing from source
    year = out_row[out_year_pos]
    if year == MISSING:
      year = gn_row[year_pos]     # possibly missing
      if year:     # MISSING is falsish
        out_row[out_year_pos] = year
        year_count += 1

    # Fill in canonical if it's missing from source (checklist_row)
    canon = out_row[out_canonical_pos]
    if canon == MISSING:
      # Use quality metric to rule out bad results from gnparse?  E.g.
      #   Hygrophorus limacinus, s.auct., non (Scop.) Fr.  -> Hygrophorus
      canon = gn_row[canonical_full_pos]
      if canon:
        out_row[out_canonical_pos] = canon
        if canon_count < CANON_SAMPLE_LIMIT:
          print("# canonical := '%s' bc '%s'" % (canon_full, gn_row[verbatim_pos]),
                file=sys.stderr)
        canon_count += 1

    # gnparser tries to parse Verbatim into CanonicalFull [+ Authorship] [+ Year].
    # The original may have junk after the CanonicalFull.
    # We further parse CanonicalStem (if provided) to Pre_epithet + epithet
    # and discard pre_epithet, otherwise just use all of Verbatim in key.
    # If Authorship is present parse it into First + rest and discard
    # the rest.
    # If Year is absent try to get it from Verbatim, and if that fails
    # use 9999 (so that the record sorts after those with years).

    stemmed = MISSING
    tipe = MISSING   # default, overridden if possible
    quality = int(gn_row[quality_pos])

    if quality == 1 or quality == 2:
      # Figure out epithet or some substitute
      # do not trim non-epithet if year or auth is missing
      stemmed = gn_row[stem_pos]
      card = int(gn_row[cardinality_pos] or '0')

      # Figure out auth part of tipe... trim off authors after first
      m = auth_re.search(gn_row[auth_pos])
      auth = m[0] if m else MISSING

      if stemmed and year and auth and card > 1:
        # Figure out epithet part of type
        # e.g. 'speciosus' in 'Hygrophorus lucorum var. speciosus'
        epithet = stemmed.split(' ')[-1]
        trim_count += 1
        # Put them together
        tipe = "TS|%s|%s|%s" % (year, epithet, auth)

    else:
      status = (checklist_row[status_pos] if status_pos != None else '')
      if False and (status == 'accepted' or status == 'valid' or status == ''):
        # Happens very frequenctly with NCBI
        print("# Poor quality name: '%s' quality %s" %
              (gn_row[verbatim_pos], quality),
              file=sys.stderr)

    # Add extra columns to the original input
    out_row = out_row + [stemmed, tipe]
    if len(out_row) != len(out_header):
      print("! %s %s" % (len(out_header), out_header,), file=sys.stderr)
      print("! %s %s" % (len(out_row), out_row,), file=sys.stderr)
      assert False
    yield out_row

  print("# use_gnparse: of %s rows, got epithet for %s, got year for %s, added canonical for %s" %
        (row_count, trim_count, year_count, canon_count),
        file=sys.stderr)

CANON_SAMPLE_LIMIT = 0

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    TBD
    """)
  parser.add_argument('--source',
                      help="name of file containing checklist")
  args=parser.parse_args()
  assert args.source
  gn_file = sys.stdin
  with open(args.source, "r") as source_file:
    better = use_parse(csv.reader(gn_file),
                       csv.reader(source_file))
    util.write_rows(better, sys.stdout)
