#!/usr/bin/env python3

# comet = ☄  erase = ⌫  recycle = ♲


"""
Here's what we're trying to do, as an n^2 method:

  For (record 1, record 2) pair, compute match score.
  Consider all candidate matches (x1, x2) where 
    either x1 = record 1 or x2 = record 2.
  Designate (record 1, record 2) a match if it has the highest score
    among all of these candidates.
  (I.e. if record 1 = unique best match to record 2 AND v.v.)

With some indexing, we can do it in approximately linear time.
"""

# -----

import sys, io, argparse, csv
from functools import reduce
from util import ingest_csv, windex, MISSING, \
                 correspondence, precolumn, apply_correspondence
from match_records import compute_sum, ingest_sum

def prepare_delta(inport1, inport2, pk_col_arg, indexed, managed_arg,
              sum_file, outport):
  global INDEX_BY, pk_pos1, pk_pos2, pk_col
  pk_col = pk_col_arg
  INDEX_BY = indexed.split(",")    # kludge

  a_table = ingest_csv(csv.reader(inport1), pk_col)
  b_table = ingest_csv(csv.reader(inport2), pk_col)

  if sum_file:
    with open(sum_file, "r") as cofile:
      cop = ingest_sum(csv.reader(cofile), pk_col)
  else:
    cop = compute_sum(a_table, b_table, pk_col, INDEX_BY)

  (header1, _) = a_table
  (header2, _) = b_table
  managed = [man
             for man in managed_arg.split(",")
             if man in header1 or man in header2]

  write_delta(cop, a_table, b_table, managed, outport)

# Write sum in the form of an EOL-ish delta

def write_delta(cop, a_table, b_table, managed, outport):
  global pk_col
  global pk_pos1, pk_pos2

  (header1, all_rows1) = a_table
  (header2, all_rows2) = b_table
  pk_pos1 = windex(header1, pk_col)
  pk_pos2 = windex(header2, pk_col)

  header3 = ['mode', 'new_pk', 'remark'] + managed
  mode_pos3 = 0
  new_pk_pos3 = 1
  remark_pos3 = 2
  # new_pk_pos3 will be the primary key within the delta itself
  old_pk_pos3 = 3 + windex(managed, pk_col)

  corr_13 = correspondence(header1, header3)
  corr_23 = correspondence(header2, header3)

  managed_positions2 = [windex(header2, name) for name in managed]

  writer = csv.writer(outport)

  def write_row(mode, row, old, new, co, remark):
    if new: assert co == new
    row[mode_pos3] = mode
    row[old_pk_pos3] = old or MISSING
    row[new_pk_pos3] = co
    row[remark_pos3] = remark
    writer.writerow(row)

  writer.writerow(header3)

  corr_12 = correspondence(header1, header2)
  add_count = 0
  carry_count = 0
  update_count = 0
  remove_count = 0
  stats = [[0, 0, 0, ([], [], [])] for col in header2]

  for (key3, (key1, key2, remark)) in cop.items():
    if key1 != None and key2 != None:
      # Carry or difference in field values
      row1 = all_rows1[key1]
      row2 = all_rows2[key2]
      if analyze_changes(row1, row2, managed_positions2, corr_12, stats):
        write_row("update", apply_correspondence(corr_23, row2),
                  key1, key2, key3, remark)
        update_count += 1
      else:
        # Carry - no output
        carry_count += 1
    elif key1 != None:
      row1 = all_rows1[key1]
      write_row("remove", apply_correspondence(corr_13, row1),
                key1, key2, key3, remark)
      remove_count += 1
    else:
      row2 = all_rows2[key2]
      write_row("add", apply_correspondence(corr_23, row2),
                key1, key2, key3, remark)
      add_count += 1

  # All done with the delta.  Reporting now.

  print("-- delta: %s carries, %s additions, %s removals, %s updates" %
        (carry_count, add_count, remove_count, update_count,),
        file=sys.stderr)
  for j in range(0, len(header2)):
    (q, c, d, (qs, cs, ds)) = stats[j]   # (set, modified, erased)
    if q > 0:
      print("--   %s: %s set" % (header2[j], q),
            file=sys.stderr)
      for row2 in qs:
        print("--     %s = %s: %s" %
              (row2[pk_pos2], row2[windex(header2, "canonicalName")], row2[j]),
              file=sys.stderr)
    if c > 0:
      print("--   %s: %s modified" % (header2[j], c),
            file=sys.stderr)
      for (row1, row2) in cs:
        print("--     %s = %s: %s -> %s" %
              (row2[pk_pos2],
               row2[windex(header2, "canonicalName")],
               row1[precolumn(corr_12, j)],
               row2[j]),
              file=sys.stderr)
    if d > 0:
      print("--   %s: %s cleared" % (header2[j], d),
            file=sys.stderr)
      for row1 in ds:
        print("--     %s = %s: %s" %
              (row1[pk_pos1],
               row1[windex(header1, "canonicalName")],
               row1[precolumn(corr_12, j)]),
              file=sys.stderr)

# for readability
SAMPLES = 3

# managed_positions2 = positions in header2 of the change-managed columns.
# Side effect: increment counters per column of added, updated, removed rows

def analyze_changes(row1, row2, managed_positions2, corr_12, stats):
  for pos2 in managed_positions2:
    if pos2 != None:
      value2 = row2[pos2]
      pos1 = precolumn(corr_12, pos2)
      value1 = MISSING
      if pos1 != None: value1 = row1[pos1]
      ss = stats[pos2]
      if value1 != value2:
        (q, c, d, (qs, cs, ds)) = ss
        if value1 == MISSING:
          ss[0] += 1            # 0 = field set
          if q < SAMPLES: qs.append(row2)
        elif value2 == MISSING: 
          ss[2] += 1            # 2 = field erased
          if d < SAMPLES: ds.append(row1)
        else:
          ss[1] += 1            # 1 = field modified
          if c < SAMPLES: cs.append((row1, row2))
        return True
  return False

# Test
def test1():
  inport1 = io.StringIO(u"taxonID,bar\n1,2")
  inport2 = io.StringIO(u"taxonID,bar\n1,3")
  prepare_delta(inport1, inport2, sys.stdout)

def test():
  inport1 = io.StringIO(u"taxonID,bar\n1,dog\n2,cat")
  inport2 = io.StringIO(u"taxonID,bar\n91,cat\n93,pig")
  prepare_delta(inport1, inport2, sys.stdout)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Standard input is file containing initial state.
    Standard output is the final state, consisting of specified state 
    file annotated with ids for initial state records, via matching.
    """)
  parser.add_argument('--target',
                      help='name of file specifying target state')
  parser.add_argument('--pk',
                      default="taxonID",
                      help='name of column containing primary key')
  # Order is important
  indexed="taxonID,scientificName,canonicalName"
  parser.add_argument('--index',
                      default=indexed,
                      help='names of columns to match on')
  managed=indexed+",taxonRank,taxonomicStatus,nomenclaturalStatus,source,datasetID"
  parser.add_argument('--manage',
                      default=managed,
                      help='names of columns under version control')
  parser.add_argument('--matches',
                      default=None,
                      help='file containing A/B matches')
  # List of fields stored in database or graphdb should be an arg.
  args=parser.parse_args()
  with open(args.target, "r") as inport2:
    prepare_delta(sys.stdin, inport2, args.pk, args.index, args.manage,
                  args.matches, sys.stdout)
