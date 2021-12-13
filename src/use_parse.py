#!/usr/bin/env python3

# Make use of the output of gnparse, by adding new columns to the checklist

import sys, csv, argparse
import regex
import util
from util import windex, MISSING

# auth_re = re.compile("[A-Z][A-Za-z'-]+")
auth_re = regex.compile(u"\p{Uppercase_Letter}[\p{Letter}-]+")

def use_parse(gn_iter, check_iter):
  header1 = next(gn_iter)
  # Id,Verbatim,Cardinality,CanonicalStem,CanonicalSimple,CanonicalFull,Authorship,Year,Quality
  if len(header1) != 9:
    print("** Incorrect number of columns in %s" % (header1,),
          file=sys.stderr)
    assert False
  stem_pos = windex(header1, "CanonicalStem")
  assert stem_pos
  canonical_pos = windex(header1, "CanonicalFull")
  auth_pos = windex(header1, "Authorship")
  year_pos = windex(header1, "Year")
  quality_pos = windex(header1, "Quality")
  header2 = next(check_iter) + ["canonicalStem", "tipe", "year"]
  yield header2
  row_count = 0
  stem_count = 0
  tipe_count = 0
  canon_count = 0
  have_pos = windex(header2, "canonical")
  for row2 in check_iter:
    row_count += 1
    row1 = next(gn_iter)
    stem = MISSING
    tipe = MISSING
    year = row1[year_pos]
    # https://github.com/gnames/gnparser/blob/master/quality.md
    if row1[quality_pos] == "1" or row1[quality_pos] == "2":
      stem = row1[stem_pos]
      if stem != MISSING:
        stem_count += 1
        # Get species or subspecies epithet stem
        epithet = stem.split(" ")[-1]
        part = auth_re.search(row1[auth_pos])
        if part: auth = part[0]
        if epithet != MISSING and auth and year != MISSING:
          tipe = "tipe[%s %s %s]" % (epithet, auth, year,)
          tipe_count += 1

      if (have_pos and row2[have_pos] == MISSING and
          row1[canonical_pos] != MISSING):
        row2[have_pos] = row1[canonical_pos]
        canon_count += 1

    yield row2 + [stem, tipe, year]
  print("# use_parse: added %s alt keys, %s stems, %s canonicals for %s rows" %
        (tipe_count, stem_count, canon_count, row_count),
        file=sys.stderr)

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
