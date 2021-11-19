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
  auth_pos = windex(header1, "Authorship")
  year_pos = windex(header1, "Year")
  quality_pos = windex(header1, "Quality")
  header2 = next(check_iter) + ["canonicalStem", "altKey"]
  yield header2
  row_count = 0
  stem_count = 0
  altkey_count = 0
  for row2 in check_iter:
    row_count += 1
    row1 = next(gn_iter)
    stem = MISSING
    altkey = MISSING
    # https://github.com/gnames/gnparser/blob/master/quality.md
    if row1[quality_pos] == "1" or row1[quality_pos] == "2":
      stem = row1[stem_pos]
      if stem != MISSING:
        stem_count += 1
        # Get species or subspecies epithet stem
        epithet = stem.split(" ")[-1]
        part = auth_re.search(row1[auth_pos])
        if part: auth = part[0]
        year = row1[year_pos]
        if epithet != MISSING and part and year != MISSING:
          altkey = "%s.%s.%s" % (epithet, auth, year)
          altkey_count += 1
    yield row2 + [stem, altkey]
  print("# use_parse: added %s alt keys and %s stems for %s rows" %
        (altkey_count, stem_count, row_count),
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
