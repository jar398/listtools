#!/usr/bin/env python3

# Apply a delta

import sys, argparse, csv
from util import windex, correspondence, apply_correspondence

mode_column = "mode"
new_pk_column = "new_pk"

# Read old state from inport
# Write new state to output

def apply_delta(inport, deltaport, pk_col, outport):
  reader1 = csv.reader(inport)
  header1 = next(reader1)
  old_pk_pos1 = windex(header1, pk_col)
  assert old_pk_pos1 != None

  reader2 = csv.reader(deltaport)
  header2 = next(reader2)
  mode_pos = windex(header2, "mode")
  old_pk_pos2 = windex(header2, pk_col)
  new_pk_pos = windex(header2, "new_pk")
  assert mode_pos != None
  assert old_pk_pos2 != None
  assert new_pk_pos != None

  header2[old_pk_pos2] = "old_pk"
  header2[new_pk_pos] = pk_col

  # row2 = row from delta.  pk column has original pk, new_pk has new pk.
  def convert_row(row2):
    row3 = row2 + []
    del row3[old_pk_pos2]    # old_pk
    del row3[mode_pos]
    return row3

  header3 = convert_row(header2)

  writer = csv.writer(outport)
  writer.writerow(header3)

  def write_row(row2):
    writer.writerow(convert_row(row2))

  # Turn a file 1 row into a delta row
  corr_13 = correspondence(header1, header3)

  # Cf. diff.py
  row1 = None
  row2 = None
  added = 0
  removed = 0
  changed = 0
  continued = 0
  count = 0

  while True:
    if row1 == None:
      try:
        row1 = next(reader1)
        if len(row1) != len(header1):
          print("** Row %s of stdin is ragged" % (count,), file=sys.stderr)
          assert False
        pk1 = row1[old_pk_pos1]
      except StopIteration:
        row1 = False
    if row2 == None:
      try:
        row2 = next(reader2)
        if len(row2) != len(header2):
          print("** Row %s of %s is ragged" % (count,), file=sys.stderr)
          assert False
        pk2 = row2[old_pk_pos2]
      except StopIteration:
        row2 = False
    assert row1 != None and row2 != None    # should be obvious
    if row1 == False and row2 == False:
      break
    if count % 500000 == 0:
      print("# apply: %s" % count, file=sys.stderr)

    if row1 and (not row2 or pk1 < pk2):
      # CARRY OVER.
      writer.writerow(apply_correspondence(corr_13, row1))
      row1 = None
      continued += 1
    elif row2 and (not row1 or pk1 > pk2):
      if row2[mode_pos] == "remove":
        # key2 is key for row in file 1 to remove... ???
        removed += 1
      elif row2[mode_pos] == "add":
        write_row(row2)
        added += 1
      else:
        print("Invalid mode %s for %s" % (row2[mode_pos], pk2),
              file=sys.stderr)
        assert False
      row2 = None
    else:
      # row1 updated -> row2
      if row2[mode_pos] == "update":
        write_row(row2)
        changed += 1
      elif row2[mode_pos] == "remove":
        # key2 is key for row in file 1 to remove
        removed += 1
      else:
        print("Invalid mode: %s for %a" % (row2[mode_pos], pk2),
              file=sys.stderr)
        assert False
      row1 = row2 = None
    count += 1

  print("# diff: Added:     %s" % added, file=sys.stderr)
  print("# diff: Removed:   %s" % removed, file=sys.stderr)
  print("# diff: Changed:   %s" % changed, file=sys.stderr)
  print("# diff: Continued: %s" % continued, file=sys.stderr)



if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Standard input is file to be updated.
    Standard output is updated version.
    """)
  parser.add_argument('--delta',
                      help='name of file specifying delta')
  parser.add_argument('--pk',
                      help='name of column containing primary key')
  # List of fields stored in database or graphdb should be an arg.
  args=parser.parse_args()
  with open(args.delta, "r") as inport2:
    apply_delta(sys.stdin, inport2, args.pk, sys.stdout)
