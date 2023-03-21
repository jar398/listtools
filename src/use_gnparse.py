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
  canonical_stem_pos = windex(gn_header, "CanonicalStem")
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

  # Add gn_full, gn_stem, gn_auth
  (out_gn_full_pos, add_gn_full) = ensure_column("gn_canonical_full")
  (out_gn_stem_pos, add_gn_stem) = ensure_column("gn_canonical_stem")
  (out_gn_auth_pos, add_gn_auth) = ensure_column("gn_authorship")

  # May need to consult the source record too
  scientific_pos = windex(checklist_header, "scientificName")
  canonical_pos = windex(checklist_header, "canonicalName")
  status_pos = windex(checklist_header, "nomenclaturalStatus")

  n_added_columns = (len(out_header) - len(checklist_header))

  row_count = 0

  yield out_header
  for checklist_row in check_iter:
    assert len(checklist_row) == len(checklist_header)
    row_count += 1
    gn_row = next(gn_iter)
    out_row = checklist_row + n_added_columns*[MISSING]

    if scientific_pos:
      sci_name = out_row[scientific_pos]
    else:
      sci_name = out_row[canonical_pos]

    # Kludge to allow gnparse to parse names of the form '? foo'
    gn_full = gn_row[canonical_full_pos]
    if gn_full.startswith('Wildcard '):
      gn_full = '?' + gn_full[8:]

    gn_stem = gn_row[canonical_stem_pos]
    if not ' ' in gn_stem:
      gn_stem = MISSING
    if gn_stem.startswith('Wildcard '):
      gn_stem = '?' + gn_stem[8:]
    if gn_full.endswith('ii') and gn_stem.endswith('i'):
      # Fixed in gnparse, but I'm still using older version
      gn_stem = gn_stem[0:-1]

    parts = parse.parse_name(sci_name,
                             gn_full=gn_full,
                             gn_stem=gn_stem,
                             gn_authorship=gn_row[auth_pos])
    token = parts.token
    year = parts.year
    #assert year == gn_year

    out_row[out_gn_full_pos] = gn_full
    out_row[out_gn_stem_pos] = gn_stem
    out_row[out_gn_auth_pos] = gn_row[auth_pos]

    yield out_row

  print("-- use_gnparse: %s rows" % (row_count,),
        file=sys.stderr)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Standard input is the output generated by `gnparse -s` for a set of
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
