#!/usr/bin/env python3

# Make use of the output of gnparse, by adding new columns to the checklist

import sys, csv, re, argparse
import util
from util import windex, MISSING

#r = re.compile("\p{Uppercase_Letter}\p{Lowercase_Letter}+")

auth_re = re.compile("[A-Z][A-Za-z'-]+")

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
    if row1[quality_pos] == "1":
      stem = row1[stem_pos]
      # Get species or subspecies epithet stem
      space = str.rfind(stem, " ")
      if space > 0:
        epithet = stem[space+1:]
      else:
        epithet = stem
      part = auth_re.search(row1[auth_pos])
      if part: part = part[0]
      year = row1[year_pos]
      if epithet != MISSING and part != MISSING and year != MISSING:
        altkey = "%s.%s.%s" % (epithet, part, year)
        altkey_count += 1
      if stem != MISSING:
        stem_count += 1
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
