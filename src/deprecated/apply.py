#!/usr/bin/env python3

# Apply a delta

import sys, argparse, csv
from util import windex, correspondence, apply_correspondence

mode_column = "mode"
new_pk_column = "new_pk"

# Read old checklist from inport (sorted by old primary key pk)
# Read delta from deltaport (sorted by old primary key pk)
# Write new state to output

def apply_delta(inport, deltaport, pk_col, outport):
  reader1 = csv.reader(inport)
  header1 = next(reader1)
  old_pk_pos1 = windex(header1, pk_col)
  assert old_pk_pos1 != None

  reader2 = csv.reader(deltaport)
  header2 = next(reader2)
  # Every delta has mode and new_pk columns, as well as a primary key
  # column (typically taxonID).
  mode_pos = windex(header2, "mode")
  new_pk_pos2 = windex(header2, "new_pk")
  remark_pos = windex(header2, "remark")
  old_pk_pos2 = windex(header2, pk_col)
  assert mode_pos != None
  assert old_pk_pos2 != None
  assert new_pk_pos2 != None

  header3 = header2 + []
  lose = sorted([mode_pos, new_pk_pos2, remark_pos])
  lose.reverse()
  for pos in lose:
    del header3[pos]
  pk_pos3 = windex(header3, pk_col)

  # Turn a file 1 row into a delta row
  corr_13 = correspondence(header1, header3)
  # This may be too clever... column new_pk in delta has to go to
  # column taxonID in the new B fie
  hacked_header3 = header3 + []
  hacked_header3[pk_pos3] = "new_pk"    # to match delta
  corr_23 = correspondence(header2, hacked_header3)

  writer = csv.writer(outport)
  writer.writerow(header3)

  def write_row(row2):
    row3 = apply_correspondence(corr_23, row2)
    writer.writerow(row3)

  # Cf. diff.py
  row1 = None
  row2 = None
  added = 0
  removed = 0
  changed = 0
  continued = 0
  count1 = count2 = 0

  while True:
    if row1 == None:
      try:
        row1 = next(reader1)
        if count1 % 500000 == 0:
          print("# apply: row %s" % count1, file=sys.stderr)
        count1 += 1
        if len(row1) != len(header1):
          print("** Row %s of stdin is ragged" % (count1,), file=sys.stderr)
          assert False
        pk1 = row1[old_pk_pos1]    # primary key in A
      except StopIteration:
        row1 = False
        pk1 = None
    if row2 == None:
      try:
        row2 = next(reader2)
        if count2 % 500000 == 0:
          print("# apply: delta %s" % count2, file=sys.stderr)
        count2 += 1
        if len(row2) != len(header2):
          print("** Row %s of %s is ragged" % (count2,), file=sys.stderr)
          assert False
        pk2 = row2[old_pk_pos2]    # primary key in A
      except StopIteration:
        row2 = False
        pk2 = None
    assert row1 != None and row2 != None    # should be obvious
    if row1 == False and row2 == False:
      break

    if row1 and (not row2 or pk1 < pk2):
      # CARRY OVER.  Primary key is unchanged.
      row3 = apply_correspondence(corr_13, row1)
      writer.writerow(row3)
      row1 = None
      continued += 1
    elif row2 and row2[mode_pos] == "add":
      # insert new row and continue around loop
      write_row(row2)
      row2 = None
      added += 1
    elif row2 and (not row1 or pk1 > pk2):
      print("Invalid mode '%s' for %s < %s (need to sort?)" % (row2[mode_pos], pk2, pk1),
            file=sys.stderr)
      assert False
    else:
      assert row1 and row2
      assert pk1 == pk2
      # row1 updated -> row2
      if row2[mode_pos] == "update":
        write_row(row2)
        row1 = row2 = None
        changed += 1
      elif row2[mode_pos] == "remove":
        # key2 is key for row in file 1 to remove
        row1 = row2 = None
        removed += 1
      else:
        print("Invalid mode %s for %s = %s" % (row2[mode_pos], pk2, pk1),
              file=sys.stderr)
        assert False

  print("-- apply: added %s, removed %s, updated %s, continued %s" %
        (added, removed, changed, continued),
        file=sys.stderr)

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
