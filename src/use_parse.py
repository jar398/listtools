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

  # gnparse output: gn_iter / header1
  header1 = next(gn_iter)
  # Id,Verbatim,Cardinality,CanonicalStem,CanonicalSimple,CanonicalFull,Authorship,Year,Quality
  if len(header1) != 9:
    print("** Unexpected number of columns in %s" % (header1,),
          file=sys.stderr)
    assert False
  verbatim_pos = windex(header1, "Verbatim")
  cardinality_pos = windex(header1, "Cardinality")
  stem_pos = windex(header1, "CanonicalStem")
  canonical_full_pos = windex(header1, "CanonicalFull")
  auth_pos = windex(header1, "Authorship")
  year_pos = windex(header1, "Year")
  quality_pos = windex(header1, "Quality")

  header2 = next(check_iter)
  # May need to consult the source record too
  add_canon = not "canonicalName" in header2
  if add_canon:
    header2 = header2 + ["canonicalName"]
  canonical_pos = windex(header2, "canonicalName")
  header2 = header2 + ["canonicalStem", "year", "type"]

  row_count = 0
  trim_count = 0
  year_count = 0
  canon_count = 0

  yield header2
  for checklist_row in check_iter:
    row_count += 1
    gn_row = next(gn_iter)

    # gnparser tries to parse Verbatim into CanonicalFull [+ Authorship] [+ Year].
    # The original may have junk after the CanonicalFull.
    # We further parse CanonicalStem (if provided) to Pre_epithet + epithet
    # and discard pre_epithet, otherwise just use all of Verbatim in key.
    # If Authorship is present parse it into First + rest and discard
    # the rest.
    # If Year is absent try to get it from Verbatim, and if that fails
    # get it use 9999 (so that the record sorts after those with years).

    # Figure out year part of tipe
    year = gn_row[year_pos]     # possibly missing
    if year:     # MISSING is falsish
      year_count += 1

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

      if stemmed and year and auth:
        if card > 1:
          # Figure out epithet part of type
          pos = stemmed.index(' ', card-1)     # location of space
          epithet = stemmed[pos+1:]         # species or subspecies epithet 
          trim_count += 1
        else:
          assert card == 1
          epithet = stemmed     # genus
        # Put them together
        tipe = "TS|%s|%s|%s" % (year, epithet, auth)

      # Extra benefit: fill in canonical if it's missing from source (checklist_row)
      full = gn_row[canonical_full_pos]
      if add_canon: full = full + [MISSING]
      if full:
        have = checklist_row[canonical_pos]
        if have == MISSING:
          quality = int(gn_row[quality_pos])
          verb = gn_row[verbatim_pos]
          if quality <= 2:
            if canon_count < CANON_SAMPLE_LIMIT:
              print("# canonical := '%s' bc '%s'" % (full, verb),
                    file=sys.stderr)
            checklist_row[canonical_pos] = full
            canon_count += 1

    # Add extra columns to the original input
    yield checklist_row + [stemmed, year, tipe]

  print("# use_parse: of %s rows, got epithet for %s, got year for %s, fixed canonical for %s" %
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
