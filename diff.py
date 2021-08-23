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
from util import read_csv, windex, MISSING, \
                 correspondence, precolumn, apply_correspondence

def matchings(inport1, inport2, pk_col, indexed, managed, outport):
  global INDEX_BY, pk_pos1, pk_pos2
  INDEX_BY = indexed.split(",")    # kludge

  (header1, all_rows1) = read_csv(inport1, pk_col)
  (header2, all_rows2) = read_csv(inport2, pk_col)

  pk_pos1 = windex(header1, pk_col)
  pk_pos2 = windex(header2, pk_col)
  assert pk_pos1 != None 
  assert pk_pos2 != None
  # Positions in header2 of the managed columns ("fields of interest")
  foi_positions = [windex(header2, name) for name in managed.split(",")]

  rows2_by_property = index_rows_by_property(all_rows2, header2)
  (best_rows_in_file1, best_rows_in_file2) = \
    find_best_matches(header1, header2, all_rows1, all_rows2,
                      pk_col, rows2_by_property)

  writer = csv.writer(outport)

  # Generate a delta row from a 2nd-input row.
  def convert_row(row2, mode, key1):
    delta_row = [mode, row2[pk_pos2]] + row2
    delta_row[pk_pos2 + 2] = key1
    return delta_row

  # Primary key is key1; key2 goes to "new_pk" column
  def write_row(mode, key1, row2):
    writer.writerow(convert_row(row2, mode, key1))

  # key1 goes into the primary key column, while
  # key2 goes into the "new_pk" column.  Other columns are those from 
  # the second file.  The output is sorted and processed according to
  # the old file's primary key.  A change to the primary key comes by
  # setting a new record's primary key to the value in the 'new_pk'
  # column.
  modified_header2 = header2 + []
  modified_header2[pk_pos2] = "new_pk"
  write_row("mode", pk_col, modified_header2)

  def flush_row(row1):
    # We don't really need the values: could say [MISSING for x in header2]
    fake = apply_correspondence(corr_12, row1)
    fake[pk_pos2] = MISSING
    write_row("remove", key1, fake)

  # Now emit the matches.  
  # Process 1st file for changes and deletions
  corr_12 = correspondence(header1, header2)
  carry_count = 0
  update_count = 0
  remove_count = 0
  stats = [[0, 0, 0, ([], [], [])] for col in header2]
  seen = {}
  for (key1, row1) in all_rows1.items():
    best_rows2 = best_rows_in_file2.get(key1)
    if best_rows2:
      (score, rows2) = best_rows2
      (matchp, mode) = check_match([row1], rows2, score,
                                   best_rows_in_file1)
      if matchp:
        row2 = rows2[0]
        seen[row2[pk_pos2]] = True
        if analyze_changes(row1, row2, foi_positions, corr_12, stats):
          write_row("update", key1, row2)
          update_count += 1
        else:
          # write_row("carry", key1, row2)
          carry_count += 1
      else:
        # Delete row1.
        # No need to report; check_match has already done that.
        flush_row(row1)
        remove_count += 1
    else:
      # Unmatched; remove
      flush_row(row1)
      remove_count += 1

  # Find additions in 2nd file
  add_count = 0
  for (key2, row2) in all_rows2.items():
    if not key2 in seen:
      if key2 in all_rows1:
        print("Collision: %s" % key2)
      write_row("add", key2, row2)
      add_count += 1
  seen = None

  print("%s carries, %s additions, %s removals, %s updates" %
        (carry_count, add_count, remove_count, update_count,),
        file=sys.stderr)
  for j in range(0, len(header2)):
    (a, c, d, (qs, cs, ds)) = stats[j]
    if a > 0:
      x = [row2[pk_pos2] for row2 in qs]
      print("  %s: %s set %s" % (header2[j], a, x),
            file=sys.stderr)
    if c > 0:
      x = [row1[pk_pos1] for row1 in cs]
      print("  %s: %s modified %s" % (header2[j], c, x),
            file=sys.stderr)
    if d > 0:
      x = [row1[pk_pos1] for row1 in ds]
      print("  %s: %s cleared %s" % (header2[j], d, x),
            file=sys.stderr)

# for readability
SAMPLES = 3

# foi_positions = positions in header2 of the change-managed columns.
# Side effect: increment counters per column of added, updated, removed rows

def analyze_changes(row1, row2, foi_positions, corr_12, stats):
  for pos2 in foi_positions:
    if pos2 != None:
      value2 = row2[pos2]
      pos1 = precolumn(corr_12, pos2)
      value1 = MISSING
      if pos1 != None: value1 = row1[pos1]
      ss = stats[pos2]
      if value1 != value2:
        (a, c, d, (qs, cs, ds)) = ss
        if value1 == MISSING:
          ss[0] += 1
          if a < SAMPLES: qs.append(row2)
        elif value2 == MISSING: 
          ss[2] += 1
          if d < SAMPLES: ds.append(row1)
        else:
          ss[1] += 1
          if c < SAMPLES: cs.append(row1)
        return True
  return False

WAD_SIZE = 4

def check_match(rows1, rows2, score, best_rows_in_file1):
  global pk_pos1, pk_pos2
  keys1 = [row1[pk_pos1] for row1 in rows1]
  keys2 = [row2[pk_pos2] for row2 in rows2]
  if len(rows1) > 1:
    if len(rows1) < WAD_SIZE:
      # keys1 = [row1[pk_pos1] for row1 in rows1]
      print("Tie: multiple old %s -> new %s (score %s)" %
            (keys1, keys2[0], score),
            file=sys.stderr)
    return (False, "contentious")
  elif len(rows2) > 1:
    if len(rows2) < WAD_SIZE:
      # does not occur in 0.9/1.1
      print("Tie: old %s -> multiple new %s (score %s)" %
            (keys1[0], keys2, score),
            file=sys.stderr)
    return (False, "ambiguous")
  elif len(rows2) == 0:
    return (False, "unevaluated")
  elif score < 100:
    print("Old %s match to new %s is too weak to use (score %s)" % (keys1[0], keys2[0], score),
          file=sys.stderr)
    return (False, "weak")
  else:
    key1 = keys1[0]
    key2 = keys2[0]
    best_rows1 = best_rows_in_file1.get(key2)
    if best_rows1:
      (score3, rows3) = best_rows1
      pk_pos1 = 0    # FIX ME FIX ME
      if len(rows3) > 1:
        if len(rows3) < WAD_SIZE:
          keys3 = [row3[pk_pos1] for row3 in rows3]
          print("Old %s (score %s) colliding at %s" % (keys3, score, key2,),
                file=sys.stderr)
        return (False, "collision")
      elif score < score3:
        # Don't report; this happens way too often
        return (False, "inferior")
      else:
        assert score == score3
        row3 = rows3[0]
        key3 = row3[pk_pos1]
        if key1 != key3:
          print("%s = %s != %s, score %s, back %s" % (key1, key2, key3, score, score3),
                file=sys.stderr)
          assert key1 == key3
        return (True, "match")
    else:
      return (False, "unconsidered")

# Positions in header2 of columns to be indexed (properties)

def indexed_positions(header, index_by):
  return [windex(header, col)
          for col in index_by
          if windex(header, col) != None]

# One weight for each column in file A

def get_weights(header_b, header_a, index_by):
  weights = [(1 if col in header_b else 0) for col in header_a]

  # Censor these
  mode_pos = windex(header_b, "mode")
  if mode_pos != None: weights[mode_pos] = 0

  w = 100
  loser = index_by + []
  loser.reverse()               # I hate this
  for col in loser:
    pos = windex(header_a, col)
    if pos != None:
      weights[pos] = w
    w = w + w
  return weights

def find_best_matches(header1, header2, all_rows1, all_rows2,
                      pk_col, rows2_by_property):
  global pk_pos1, pk_pos2
  assert len(all_rows2) > 0
  corr_12 = correspondence(header1, header2)
  print("Correspondence: %s" % (corr_12,), file=sys.stderr)
  positions = indexed_positions(header1, INDEX_BY)
  print("Indexed: %s" % positions, file=sys.stderr)
  weights = get_weights(header1, header2, INDEX_BY)    # parallel to header2
  print("Weights: %s" % weights, file=sys.stderr)
  no_info = (-1, [])

  best_rows_in_file2 = {}    # key1 -> (score, rows2)
  best_rows_in_file1 = {}    # key2 -> (score, rows1)
  prop_count = 0
  for (key1, row1) in all_rows1.items():
    # The following check is also enforced by start.py... flush them here?
    best2_so_far = no_info
    best_rows_so_far2 = no_info

    for prop in row_properties(row1, header1, positions):
      if prop_count % 500000 == 0:
        print(prop_count, file=sys.stderr)
      prop_count += 1
      for row2 in rows2_by_property.get(prop, []):
        key2 = row2[pk_pos2]
        score = compute_score(row1, row2, corr_12, weights)
        best_rows_so_far1 = best_rows_in_file1.get(key2, no_info)

        # Update best file2 match for row1
        (score2_so_far, rows2) = best_rows_so_far2
        if score > score2_so_far:
          best_rows_so_far2 = (score, [row2])
        elif score == score2_so_far and not row2 in rows2:
          if len(rows2) < WAD_SIZE: rows2.append(row2)

        # Update best file1 match for row2
        (score1_so_far, rows1) = best_rows_so_far1
        if score > score1_so_far:
          best_rows_in_file1[key2] = (score, [row1])
        elif score == score1_so_far and not row1 in rows1:
          if len(rows1) < WAD_SIZE: rows1.append(row1)

    if best_rows_so_far2 != no_info:
      best_rows_in_file2[key1] = best_rows_so_far2

  print("%s properties" % prop_count, file=sys.stderr)
  if len(all_rows1) > 0 and len(all_rows2) > 0:
    assert len(best_rows_in_file1) > 0
    assert len(best_rows_in_file2) > 0
  return (best_rows_in_file1, best_rows_in_file2)

def compute_score(row1, row2, corr_12, weights):
  s = 0
  for j in range(0, len(row2)):
    w = weights[j]
    if w != 0:
      v2 = row2[j]
      i = precolumn(corr_12, j)
      if i != None:
        v1 = row1[i]
      else:
        v1 = MISSING
      if v1 == MISSING or v2 == MISSING:
        d = 1
      elif v1 == v2:
        d = 100
      else:
        d = 0
      s += w * d
  return s

LIMIT=100

def index_rows_by_property(all_rows, header):
  positions = indexed_positions(header, INDEX_BY)
  rows_by_property = {}
  entry_count = 0
  for (key, row) in all_rows.items():
    for property in row_properties(row, header, positions):
      rows = rows_by_property.get(property)
      if rows != None:
        if len(rows) <= LIMIT:
          if len(rows) == LIMIT:
            print("%s+ rows with property %s" % (LIMIT, property,),
                  file=sys.stderr)
          rows.append(row)
          entry_count += 1
      else:
        rows_by_property[property] = [row]
        entry_count += 1
  print("%s properties" % (len(rows_by_property),),
        file=sys.stderr)
  return rows_by_property

# Future: exclude really ephemeral properties like taxonID

def row_properties(row, header, positions):
  return [(header[i], row[i])
          for i in range(0, len(header))
          if (i in positions and
              row[i] != MISSING)]

# Test
def test1():
  inport1 = io.StringIO(u"taxonID,bar\n1,2")
  inport2 = io.StringIO(u"taxonID,bar\n1,3")
  matchings(inport1, inport2, sys.stdout)

def test():
  inport1 = io.StringIO(u"taxonID,bar\n1,dog\n2,cat")
  inport2 = io.StringIO(u"taxonID,bar\n91,cat\n93,pig")
  matchings(inport1, inport2, sys.stdout)


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
  indexed="taxonID,EOLid,scientificName,canonicalName"
  parser.add_argument('--index',
                      default=indexed,
                      help='names of columns to match on')
  managed=indexed+",taxonRank,taxonomicStatus,datasetID"
  parser.add_argument('--manage',
                      default=managed,
                      help='names of columns under version control')
  # List of fields stored in database or graphdb should be an arg.
  args=parser.parse_args()
  with open(args.target, "r") as inport2:
    matchings(sys.stdin, inport2, args.pk, args.index, args.manage, sys.stdout)
