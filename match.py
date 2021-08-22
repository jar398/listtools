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

# CONFIGURATION:

INDEX_BY = \
  ["EOLid", "source", "scientificName", "canonicalName", "taxonID"]

# Fields sufficient to say it's a continuation not an update...
#  kludge for testing purposes, really all columns should matter
INTERESTING_FIELDS = \
  ["canonicalName", "scientificName"]

# -----

import sys, io, argparse, csv
from functools import reduce
from util import read_csv, windex, MISSING, \
                 correspondence, precolumn, apply_correspondence

def matchings(inport1, inport2, pk_col, indexed, managed, outport):
  global INDEX_BY
  INDEX_BY = indexed.split(",")    # kludge

  (header1, all_rows1) = read_csv(inport1, pk_col)
  (header2, all_rows2) = read_csv(inport2, pk_col)

  pk_pos1 = windex(header1, pk_col)
  pk_pos2 = windex(header2, pk_col)
  assert pk_pos1 != None 
  assert pk_pos2 != None
  foi_positions = [windex(header2, name) for name in managed.split(",")]

  rows2_by_property = get_rows_by_property(all_rows2, header2)
  (best_in_file1, best_in_file2) = \
    find_best_matches(header1, header2, all_rows1, all_rows2,
                      pk_col, rows2_by_property)

  writer = csv.writer(outport)

  # Generate a delta row from a 2nd-input row.
  def convert_row(row2, mode, key1):
    row3 = [mode, row2[pk_pos2]] + row2
    row3[pk_pos2 + 2] = key1
    return row3

  # Primary key is key1; key2 goes to "new_pk" column
  def write_row(row2, mode, key1):
    writer.writerow(convert_row(row2, mode, key1))

  # key1 goes into the primary key column, while
  # key2 goes into the "new_pk" column.  Other columns are those from 
  # the second file.  The output is sorted and processed according to
  # the old file's primary key.  A change to the primary key comes by
  # setting a new record's primary key to the value in the 'new_pk'
  # column.
  modified_header2 = header2 + []
  modified_header2[pk_pos2] = "new_pk"
  write_row(modified_header2, "mode", pk_col)

  def flush_row(row1):
    # We don't really need the values: could say [MISSING for x in header2]
    fake = apply_correspondence(corr_12, row1)
    fake[pk_pos2] = MISSING
    write_row(fake, "remove", key1)

  # Now emit the matches.  
  # Process 1st file for changes and deletions
  corr_12 = correspondence(header1, header2)
  carry_count = 0
  update_count = 0
  remove_count = 0
  for (key1, row1) in all_rows1.items():
    best2 = best_in_file2.get(key1)
    if best2:
      (score, keys2) = best2
      (matchp, mode) = check_match([key1], keys2, score, best_in_file1)
      if matchp:
        key2 = keys2[0]
        row2 = all_rows2[key2]
        if unchanged(row1, row2, foi_positions, corr_12):
          # write_row(row2, "carry", key1)
          carry_count += 1
        else:
          write_row(row2, "update", key1)
          update_count += 1
      else:
        # Delete row1.
        # No need to report; check_match has already done that.
        flush_row(row1)
        remove_count += 1
    else:
      # Unmatched; remove
      flush_row(row1)
      remove_count += 1

  print("%s carries, %s updates, %s removals" %
        (carry_count, update_count, remove_count,),
        file=sys.stderr)

  # Find additions in 2nd file
  add_count = 0
  for (key2, row2) in all_rows2.items():
    best1 = best_in_file1.get(key2)
    if best1:
      (score, keys1) = best1
      (matchp, mode) = check_match(keys1, [key2], score, best_in_file1)
      if matchp:
        pass                    # covered in previous loop
      else:
        # Condition has already been reported.  No report here.
        write_row(row2, "add", key1)    # ?
        add_count += 1
    else:
      # ?  need a key1-like key for sorting. use key2, tweaked in
      # order to avoid possible collision.
      write_row(row2, "add", key2 + '+')
      add_count += 1

  print("%s addition" % (add_count,), file=sys.stderr)

def unchanged(row1, row2, foi_positions, corr_12):
  for pos2 in foi_positions:
    if pos2 != None:
      pos1 = precolumn(corr_12, pos2)
      value1 = MISSING
      if pos1 != None: value1 = row1[pos1]
      if value1 != row2[pos2]:
        return False
  return True

WAD_SIZE = 4

def check_match(keys1, keys2, score, best_in_file1):
  if len(keys1) > 1:
    if len(keys1) < WAD_SIZE:
      print("Tie: multiple old %s -> new %s score %s" %
            (keys1, keys2[0], score),
            file=sys.stderr)
    return (False, "contentious")
  elif len(keys2) > 1:
    if len(keys2) < WAD_SIZE:
      # does not occur in 0.9/1.1
      print("Tie: old %s -> multiple new %s score %s" %
            (keys1[0], keys2, score),
            file=sys.stderr)
    return (False, "ambiguous")
  elif len(keys2) == 0:
    return (False, "unevaluated")
  elif score < 100:
    print("Old %s match to new %s is too weak to use (score %s)" % (keys1[0], keys2[0], score),
          file=sys.stderr)
    return (False, "weak")
  else:
    key1 = keys1[0]
    key2 = keys2[0]
    best1 = best_in_file1.get(key2)
    if best1:
      (score3, keys3) = best1
      if len(keys3) > 1:
        if len(keys3) < WAD_SIZE:
          print("Old %s tied (score %s) with other old(s) for new %s" % (key1, score, key2,),
                file=sys.stderr)
        return (False, "tie")
      else:
        key3 = keys3[0]
        if score != score3:
          print("Old %s (score %s) defeated by old %s (score %s) for new %s" % (key1, score, key3, score3, key2,),
                file=sys.stderr)
          return (False, "defeated")
        else:
          if key1 != key3:
            print("%s = %s != %s, score %s, back %s" % (key1, key2, key3, score, score3),
                  file=sys.stderr)
            assert key1 == key3
          return (True, "match")
    else:
      return (False, "unconsidered")


def indexed_positions(header):
  return [windex(header, col)
          for col in INDEX_BY
          if windex(header, col) != None]

# One weight for each column in file A

def get_weights(header_a, header_b):
  weights = [(1 if col in header_b else 0) for col in header_a]

  # Censor these
  mode_pos = windex(header_b, "mode")
  if mode_pos != None: weights[mode_pos] = 0

  w = 100
  loser = INDEX_BY + []
  loser.reverse()               # I hate this
  for col in loser:
    pos = windex(header_a, col)
    if pos != None:
      weights[pos] = w
    w = w + w
  return weights

def find_best_matches(header1, header2, all_rows1, all_rows2,
                      pk_col, rows2_by_property):
  assert len(all_rows2) > 0
  corr_12 = correspondence(header1, header2)
  print("Correspondence: %s" % (corr_12,), file=sys.stderr)
  positions = indexed_positions(header1)
  print("Indexed: %s" % positions, file=sys.stderr)
  weights = get_weights(header2, header1)
  print("Weights: %s" % weights, file=sys.stderr)
  pk_pos = windex(header1, pk_col)
  if pk_pos == None:
    print("No %s column in %s" % (pk_col, header1,),
          file=sys.stderr)
  assert pk_pos != None
  no_info = (-1, [])

  best_in_file2 = {}    # key1 -> (score, key2)
  best_in_file1 = {}    # key2 -> (score, key1)
  prop_count = 0
  for (key1, row1) in all_rows1.items():
    # The following check is also enforced by start.py... flush them here?
    best2_so_far = no_info

    for prop in row_properties(row1, header1, positions):
      if prop_count % 100000 == 0:
        print(prop_count, file=sys.stderr)
      prop_count += 1
      for key2 in rows2_by_property.get(prop, []):
        row2 = all_rows2[key2]
        score = compute_score(row1, row2, corr_12, weights)
        best1_so_far = best_in_file1.get(key2, no_info)

        # Update best file2 match for row1
        (score2_so_far, keys2) = best2_so_far
        if score > score2_so_far:
          best2_so_far = (score, [key2])
        elif score == score2_so_far and not key2 in keys2:
          if len(keys2) < WAD_SIZE: keys2.append(key2)

        # Update best file1 match for row2
        (score1_so_far, keys1) = best1_so_far
        if score > score1_so_far:
          best_in_file1[key2] = (score, [key1])
        elif score == score1_so_far and not key1 in keys1:
          if len(keys1) < WAD_SIZE: keys1.append(key1)

    if best2_so_far != no_info:
      best_in_file2[key1] = best2_so_far

  print("%s properties" % prop_count, file=sys.stderr)
  if len(all_rows1) > 0:
    assert len(best_in_file2) > 0
  return (best_in_file1, best_in_file2)

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

def get_rows_by_property(all_rows, header):
  positions = indexed_positions(header)
  by_property = {}
  entry_count = 0
  for (key, row) in all_rows.items():
    for property in row_properties(row, header, positions):
      keys = by_property.get(property)
      if keys != None:
        if len(keys) <= LIMIT:
          if len(keys) == LIMIT:
            print("%s+ rows with property %s" % (LIMIT, property,),
                  file=sys.stderr)
          keys.append(key)
          entry_count += 1
      else:
        by_property[property] = [key]
        entry_count += 1
  print("%s properties" % (len(by_property),),
        file=sys.stderr)
  return by_property

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
  indexed="EOLid,source,scientificName,canonicalName,taxonID"
  parser.add_argument('--index',
                      default=indexed,
                      help='names of columns to match on')
  managed=indexed+"taxonRank,taxonomicStatus,datasetID"
  parser.add_argument('--manage',
                      default=managed,
                      help='names of columns under version control')
  # List of fields stored in database or graphdb should be an arg.
  args=parser.parse_args()
  with open(args.target, "r") as inport2:
    matchings(sys.stdin, inport2, args.pk, args.index, args.manage, sys.stdout)
