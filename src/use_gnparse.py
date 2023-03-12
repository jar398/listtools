#!/usr/bin/env python3

# Make use of the output of gnparse, by adding new columns to the checklist

import sys, csv, argparse
import regex
import util, parse
from util import windex, MISSING, log

# was auth_re = re.compile("[A-Z][A-Za-z'-]+")
auth_re = regex.compile(u"\p{Uppercase_Letter}[\p{Letter}-]+")
year_re = regex.compile(' ([12][0-9]{3})\)?$')

CANON_SAMPLE_LIMIT = 0    # for debugging
ADD_TIPES = True

def use_parse(gn_iter, check_iter):

  # gnparse output: gn_iter / gn_header
  gn_header = next(gn_iter)
  # Id,Verbatim,Cardinality,CanonicalStem,CanonicalSimple,CanonicalFull,Authorship,Year,Quality
  if len(gn_header) != 9:
    print("** Expected 9 columns in gnparse output, but got %s" % (len(gn_header),),
          file=sys.stderr)
    assert False
  verbatim_pos = windex(gn_header, "Verbatim")
  cardinality_pos = windex(gn_header, "Cardinality")
  stem_pos = windex(gn_header, "CanonicalStem")
  canonical_full_pos = windex(gn_header, "CanonicalFull")
  canonical_simple_pos = windex(gn_header, "CanonicalSimple")
  auth_pos = windex(gn_header, "Authorship")
  year_pos = windex(gn_header, "Year")
  quality_pos = windex(gn_header, "Quality")

  checklist_header = next(check_iter)
  out_header = checklist_header + []

  def ensure_column(col):
    pos = windex(out_header, col)
    if pos != None:
      return (pos, False)
    else:
      pos = len(out_header)
      out_header.append(col)
      return (pos, True)      

  (out_canon_pos, add_canon) = ensure_column("canonicalName")
  (out_year_pos, add_year) = ensure_column("namePublishedInYear")
  (out_stem_pos, add_stem) = ensure_column("canonicalStem")
  (out_auth_pos, add_auth) = ensure_column("scientificNameAuthorship")
  if ADD_TIPES:
    (out_tipe_pos, _) = ensure_column("tipe")

  # May need to consult the source record too
  status_pos = windex(checklist_header, "nomenclaturalStatus")

  out_scientific_pos = windex(out_header, "scientificName")
  out_year_pos = windex(out_header, "namePublishedInYear")

  row_count = 0
  tipe_count = 0
  year_count = 0
  canon_count = 0
  stemmed_count = 0
  lose_count = 0

  n_added_columns = (len(out_header) - len(checklist_header))

  yield out_header
  for checklist_row in check_iter:
    assert len(checklist_row) == len(checklist_header)
    row_count += 1
    gn_row = next(gn_iter)
    out_row = checklist_row + n_added_columns * [MISSING]

    sci_name = row[out_scientific_pos]
    if sci_name.startswith("Nil"):
      sci_name = "?" + sci_name[3:]
    (complete, authorship) = parse.split_protonymic(sci_name,
                                                    gn_row[canonical_full_pos],
                                                    gn_row[auth_pos])
    (token, year, _) = parse.analyze_authorship(authorship)
    assert year == gn_year

    # Fill in year if it's missing from source
    year = out_row[out_year_pos]
    if year == MISSING:
      year = gn_row[year_pos]     # possibly missing
      if year:     # MISSING is falsish
        out_row[out_year_pos] = year
        year_count += 1

    # gnparser tries to parse Verbatim into
    #   CanonicalFull [+ Junk] [+ Authorship].
    # where Junk is whatever sits between the CanonicalFull and the year.
    # We further parse CanonicalStem (if gnparse provides it) to Pre_epithets + epithets
    # and discard pre_epithets, otherwise just use all of Verbatim in key.
    # If Authorship is present parse it into Most + last, keep last, and discard
    # Most.  (The last epithet is most informative.)
    # Challenge:
    #   gnparser "Foo bar baz quux var. yow Jones, 2008"
    # which yields an empty Authorship and Year and Quality = 5.

    stemmed = MISSING
    tipe = MISSING   # default, overridden if possible
    quality = int(gn_row[quality_pos])

    # Use quality metric to rule out bad results from gnparse.  E.g.
    # if row has scientific "Hygrophorus limacinus, s.auct., non (Scop.) Fr."
    # then gnparser will say the canonicalName is "Hygrophorus" - bad.

    if quality < 4:  # was quality == 1 or quality == 2:

      # gnparser canonical.
      # Fill in canonical if it's missing from source (checklist_row)
      canon = out_row[out_canon_pos]
      if canon == MISSING:
        canon = gn_row[canonical_full_pos]
        if canon:
          out_row[out_canon_pos] = canon
          if canon_count < CANON_SAMPLE_LIMIT:
            print("# canonical := '%s' bc '%s'" % (canon_full, gn_row[verbatim_pos]),
                  file=sys.stderr)
          canon_count += 1

      # Figure out epithet or some substitute
      # do not trim non-epithet if year or auth is missing
      if gn_row[canonical_full_pos] == gn_row[canonical_simple_pos]:
        stemmed = gn_row[stem_pos]
        (_, _, epithet) = parse.analyze_complete(complete)   # but want the epithet to be stemmed...
        if not epithet.startswith(stemmed):
          # There's extra stuff after the end of the canonicalFull and canonicalStemmed,
          # so the end of canonicalStemmed is not going to give a solid match key.
          log("# Substituting nonstandard %s for gnparse's %s" %
              (epithet, stemmed))
          stemmed = epithet

        # See https://github.com/gnames/gnparser/issues/238
        if stemmed.endswith('i') and canon.endswith('ii'):
          stemmed = stemmed[0:-1]

      if ADD_TIPES:
        # Number of space-separated parts of name (e.g. Genus epithet = 2)
        card = int(gn_row[cardinality_pos] or '0')

        # Figure out auth part of tipe... trim off authors after first
        authorship = gn_row[auth_pos]
        m = auth_re.search(authorship)
        token = m[0] if m else MISSING

        if stemmed and year and token and card > 1:
          # Figure out epithet part of type
          # e.g. 'specios' in 'Hygrophorus lucorum var. speciosus'
          epithet = stemmed.split(' ')[-1]
          tipe_count += 1
          # Put them together
          if ';' in epithet or ';' in token:
            print("!!! Warning: semicolon in tipe should be quoted somehow",
                  file=sys.stderr)
          assert not ' ' in epithet
          assert not ' ' in token
          tipe = "[? %s %s, %s]" % (epithet, token, year)

        out_row[out_tipe_pos] = tipe

    else:
      status = (checklist_row[status_pos] if status_pos != None else '')
      if False and (status.startswith('accepted') or status == 'valid' or status == ''):
        # Happens very frequently with NCBI.  Do our own parse?
        print("# Poor quality name: '%s' quality %s" %
              (gn_row[verbatim_pos], quality),
              file=sys.stderr)

    # Add extra columns to the original input
    if stemmed: stemmed_count += 1
    out_row[out_stem_pos] = stemmed
    out_row[out_auth_pos] = authorship
    if len(out_row) != len(out_header):
      print("! %s %s" % (len(out_header), out_header,), file=sys.stderr)
      print("! %s %s" % (len(out_row), out_row,), file=sys.stderr)
      assert False
    yield out_row

    # Kill scientificNames that are redundant with their canonicalName
    if (out_row[out_canon_pos] != MISSING and
        (out_row[out_canon_pos] == sci_name)):
      out_row[out_scientific_pos] = MISSING
      lose_count += 1

  print("# use_gnparse: %s rows, added canonical %s, stemmed %s, year %s, tipe %s, lose %s" %
        (row_count, canon_count, stemmed_count, year_count,
         tipe_count, lose_count),
        file=sys.stderr)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Standard input is output generated by gnparse -s for a set of
    'scientific names.'
    Standard output is formed by combining the gnparse output with the
    parallel set of rows from the SOURCE Darwin Core checklist to 
    form a new checklist with a few additional columns.
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
