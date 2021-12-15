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
  header1 = next(gn_iter)
  # Id,Verbatim,Cardinality,CanonicalStem,CanonicalSimple,CanonicalFull,Authorship,Year,Quality
  if len(header1) != 9:
    print("** Unexpected number of columns in %s" % (header1,),
          file=sys.stderr)
    assert False
  cardinality_pos = windex(header1, "Cardinality")
  stem_pos = windex(header1, "CanonicalStem")
  verbatim_pos = windex(header1, "Verbatim")
  canonical_full_pos = windex(header1, "CanonicalFull")
  auth_pos = windex(header1, "Authorship")
  year_pos = windex(header1, "Year")
  quality_pos = windex(header1, "Quality")

  # May need to consult the source record too
  header2 = next(check_iter) + ["canonicalStem", "year", "tipe"]
  id_pos = windex(header2, "taxonID")
  canonical_pos = windex(header2, "canonicalName")
  scientific_pos = windex(header2, "scientificName")
  record_pos = windex(header2, "record_id")

  row_count = 0
  trim_count = 0
  year_count = 0
  canon_count = 0

  yield header2
  for row2 in check_iter:
    row_count += 1
    gn_row = next(gn_iter)

    # gnparser tries to parse Verbatim into CanonicalFull [+ Authorship] [+ Year].
    # We further parse CanonicalFull (if provided) to Pre_epithet + epithet
    # and discard pre_epithet, otherwise just use all of Verbatim in key.
    # If Authorship is present parse it into First + rest and discard
    # the rest.
    # If Year is absent try to get it from Verbatim, and if that fails
    # get it use 9999 (so that the record sorts after those with years).

    # Figure out year part of tipe
    year = gn_row[year_pos]
    if year:
      year_count += 1
      ersatz_year = year
    else:
      m = year_re.search(year)  # rejected by gnparse
      if m:
        year_count += 1
        ersatz_year = m[1]
      else:
        ersatz_year = '9999'

    # Figure out auth part of tipe
    m = auth_re.search(gn_row[auth_pos])
    auth = m[0] if m else MISSING

    # Figure out epithet or some substitute
    # do not trim non-epithet if year or auth is missing
    stem = gn_row[stem_pos]
    card = int(gn_row[cardinality_pos] or '0')
    if stem and card > 1 and year and auth:
      # Figure out epithet part of type
      pos = stem.index(' ', card-1)     # location of space
      epithet = stem[pos+1:]         # species or subspecies epithet 
      trim_count += 1
    else:  # rejected by gnparse
      # Use other information
      epithet = stem
      if not epithet: epithet = gn_row[verbatim_pos]
      if not epithet and canonical_pos:  epithet = row2[canonical_pos]
      if not epithet and scientific_pos: epithet = row2[id_pos]
      if not epithet and record_pos:     epithet = row2[record_pos]
      if not epithet: epithet = row2[id_pos]
      assert epithet
      
    # Put them together
    tipe = "%s|%s|%s" % (ersatz_year, epithet, auth,)

    # Extra benefit: fill in canonical if it's missing from source (row2)
    full = gn_row[canonical_full_pos]
    if full and canonical_pos:
      have = row2[canonical_pos]
      if have == MISSING:
        quality = int(gn_row[quality_pos])
        verb = gn_row[verbatim_pos]
        if quality <= 2:
          if canon_count < CANON_SAMPLE_LIMIT:
            print("# canonical := '%s' bc '%s'" % (full, verb),
                  file=sys.stderr)
          row2[canonical_pos] = full
          canon_count += 1

    yield row2 + [stem, year, tipe]

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
    util.write_generated(better, sys.stdout)
