#!/usr/bin/env python3

# This should always be the first step in a processing pipeline.
#  - Normalizes tsv to csv
#  - Makes sure all rows have the same number of fields

import sys, csv, re, argparse

def start_csv(filename, outport, cleanp):
  with open(filename, "r") as inport:
    (d, q, g) = csv_parameters(filename)
    reader = csv.reader(inport, delimiter=d, quotechar=q, quoting=g)
    header = next(reader)
    if len(header) == 1:
      if "," in header[0] or "\t" in header[0]:
        print("** start: Suspicious header", file=sys.stderr)
        print("** start: Header is %s" % (row,), file=sys.stderr)
    can_pos = windex(header, "canonicalName")
    sci_pos = windex(header, "scientificName")
    source_pos = windex(header, "source")
    landmark_pos = windex(header, "Landmark")
    taxon_id_pos = windex(header, "taxonID")
    if landmark_pos != None: header[landmark_pos] = "landmark_status"
    writer = csv.writer(outport) # CSV not TSV
    writer.writerow(header)
    count = 0
    trimmed = 0
    clean_count = 0
    for row in reader:
      if len(row) > len(header):
        row = row[0:len(header)]
        trimmed += 1
      if len(row) != len(header):
        print(("** start: Unexpected number of columns: have %s want %s" %
               (len(row), len(header))),
              file=sys.stderr)
        print(("** start: Row is %s" % (row,)), file=sys.stderr)
        assert False
      if cleanp and can_pos != None and sci_pos != None:
        if clean(row, can_pos, sci_pos):
          clean_count += 1

      if landmark_pos != None: 
        l = row[landmark_pos]
        if l != None and l != '':
          e = int(l)
          # enum landmark: %i[no_landmark minimal abbreviated extended full]
          if   e == 1: row[landmark_pos] = 'minimal'
          elif e == 2: row[landmark_pos] = 'abbreviated'
          elif e == 3: row[landmark_pos] = 'extended'
          elif e == 4: row[landmark_pos] = 'full'
          else: row[landmark_pos] = None
      if source_pos != None and row[source_pos]:
        sources = row[source_pos].split(',')
        if len(sources) > 1 and ':' in sources[0] and ':' in sources[1]:
          row[source_pos] = sources[0]
      # Assign taxon ids to nodes (usually synonyms) that don't have them
      if taxon_id_pos != None and not row[taxon_id_pos]:
        row[taxon_id_pos] = "start.py.%s" % (count,)
      writer.writerow(row)
      count += 1
  print("start: %s rows, %s columns, %s cleaned" %
        (count, len(header), clean_count),
        file=sys.stderr)
  if trimmed > 0:
    # Ignoring extra values is appropriate behavior for DH 0.9.  But
    # elsewhere we might want ragged input to be treated as an error.
    print("start: trimmed extra values from %s rows" % (trimmed,),
          file=sys.stderr)
    
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

def clean(row, can_pos, sci_pos):
  c = row[can_pos]
  s = row[sci_pos]
  if s == None or s == '':
    # if is_scientific(c): row[sci_pos] = c
    return False
  if is_scientific(s):
    # if not c: row[can_pos] = s
    return False
  # s is nonnull and not 'scientific'
  if c == None or c == '':
    # swap
    row[sci_pos] = None
    row[can_pos] = s
    # print("start: c := s", file=sys.stderr) - frequent in DH 1.1
    return True
  if is_scientific(c):
    # should gnparse!
    if c.startswith(s):
      row[sci_pos] = c
      row[can_pos] = s
      # never in DH 1.1
      print("start: swap s/c", file=sys.stderr)
      return True
    return False
  if c == s:
    row[sci_pos] = None
    # print("start: flush s", file=sys.stderr) - happens all the time in 1.1
    return True
  return False

sci_re = re.compile(", [0-9]{4}\\b")

def is_scientific(name):
  sci_re.search(name)

def windex(header, fieldname):
  if fieldname in header:
    return header.index(fieldname)
  else:
    return None

def csv_parameters(path):
  if ".csv" in path:
    return (",", '"', csv.QUOTE_MINIMAL)
  else:
    return ("\t", "\a", csv.QUOTE_NONE)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Normalize to csv and check that number of columns is same in every row.
    CSV rows are written to standard output.
    """)
  parser.add_argument('--input', default=None,
                      help='name of input file.  TSV assumed unless name contains ".csv"')

  parser.add_argument('--clean', dest='clean', action='store_true',
                      help='clean up scientificName and canonicalName a little bit')
  parser.add_argument('--no-clean', dest='clean', action='store_false')
  parser.set_defaults(clean=True)

  args=parser.parse_args()
  start_csv(args.input, sys.stdout, args.clean)
