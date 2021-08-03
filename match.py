#!/usr/bin/env python3

# comet = ☄  erase = ⌫  recycle = ♲


"""
Here's what we're trying to do, as an n^2 method:

  For (record 1, record 2) pair, compute match score.
  Consider all matches (x1, x2) where 
    either x1 = record 1 or x2 = record 2.
  Designate (record 1, record 2) a match if it has the highest score
    among all of these (x1, x2) pairs.
  (I.e. if record 1 = unique best match to record 2 AND v.v.)

With some indexing, we can do it in approximately linear time.
"""

# CONFIGURATION:

# Most important column first
INDEX_BY = \
  ["EOLid", "source", "scientificName", "canonicalName"]

MANAGED_FIELDS = \
  ["canonicalName", "scientificName"]
"""
  ["EOLid", "source", "scientificName", "canonicalName"] + \
  ["taxonRank", "taxonomicStatus", "landmark_status"]
"""

# -----

import sys, io, argparse, csv
from functools import reduce
from util import read_csv, windex, map_row, MISSING

AMBIGUOUS = '¿'

def matchings(inport1, inport2, pk_col, outport):
  (header1, all_rows1) = read_csv(inport1, pk_col)
  (header2, all_rows2) = read_csv(inport2, pk_col)

  pk_pos1 = windex(header1, pk_col)
  pk_pos2 = windex(header2, pk_col)
  mode_pos = windex(header2, "mode")
  previous_pos = windex(header2, "previous_pk")
  foi_positions = [windex(header2, name) for name in MANAGED_FIELDS]

  rows2_by_property = get_rows_by_property(all_rows2, header2)
  (best_in_file1, best_in_file2) = \
    find_best_matches(header1, header2, all_rows1, all_rows2,
                      pk_col, rows2_by_property)

  def item_sort_key(item):
    (key1, (score, key2)) = item
    return (score, key1, key2)
  writer = csv.writer(outport)

  def write_row(row2, mode, previous):
    if previous_pos == None:
      row2 = [previous] + row2
    else:
      row2[previous_pos] = previous
    if mode_pos == None:
      row2 = [mode] + row2
    else:
      row2[mode_pos] = mode
    writer.writerow(row2)

  # Now emit the matches.
  write_row(header2, "mode", "previous_pk")

  # Note deletions of old rows
  correspondence = [windex(header1, column) for column in header2]
  remove_count = 0
  for (key1, row1) in sorted(all_rows1.items(),
                             key=lambda item: item[0]):
    best2 = best_in_file2.get(key1)
    if best2:
      (score, key2) = best2
      (matchp, mode) = check_match(key1, key2, score, best_in_file1)
      if not matchp:
        # Delete row1.  Kept rows are taken care of next
        # We don't need the values: [MISSING for x in header2]
        row2 = map_row(correspondence, row1)
        row2[pk_pos2] = MISSING
        write_row(row2, "remove", key1)
        remove_count += 1
  print("%s removals" % (remove_count,), file=sys.stderr)

  # Carry over everything from file2
  carry_count = 0
  update_count = 0
  add_count = 0
  for (key2, row2) in sorted(all_rows2.items(),
                             key=lambda item: item[0]):
    best1 = best_in_file1.get(key2)
    if best1:
      (score, key1) = best1
      if key1 != AMBIGUOUS:
        (matchp, mode) = check_match(key1, key2, score, best_in_file1)
        if matchp:
          if unchanged(all_rows1[key1], all_rows2[key2], foi_positions, correspondence):
            write_row(row2, "carry", key1)
            carry_count += 1
          else:
            write_row(row2, "update", key1)
            update_count += 1
        else:
          print("Missed match %s to %s because %s" % (key1, key2, mode), file=sys.stderr)
          write_row(row2, "add", key1)
          add_count += 1
    else:
      # Wholly new
      write_row(row2, "add", MISSING)
      add_count += 1
  print("%s carries, %s updates, %s additions" % (carry_count, update_count, add_count,), file=sys.stderr)

def unchanged(row1, row2, foi_positions, correspondence):
  for pos2 in foi_positions:
    if pos2 != None:
      pos1 = correspondence[pos2]
      value1 = MISSING
      if pos1 != None: value1 = row1[pos1]
      if value1 != row2[pos2]:
        return False
  return True

def check_match(key1, key2, score, best_in_file1):
  assert key1 != AMBIGUOUS
  if key2 == AMBIGUOUS:
    # does not occur in 0.9/1.1
    print("Tie: id %s -> multiple rows" % (key1,),
          file=sys.stderr)
    return (False, "ambiguous")
  elif key2 == None:
    return (False, "unevaluated")
  elif score < 100:
    print("Old %s match to new %s is too weak to use (score %s)" % (key1, key2, score),
          file=sys.stderr)
    return (False, "weak")
  else:
    best1 = best_in_file1.get(key2)
    if best1:
      (score3, key3) = best1
      if score != score3:
        print("Old %s (score %s) defeated by old %s (score %s) for new %s" % (key1, score, key3, score3, key2,),
              file=sys.stderr)
        return (False, "defeated")
      elif key3 == AMBIGUOUS:
        print("Old %s tied (score %s) with other old(s) for new %s" % (key1, score, key2,),
              file=sys.stderr)
        return (False, "tie")
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

def get_weights(header, header_b):
  weights = [(1 if col in header_b else 0) for col in header]

  # Censor these
  mode_pos = windex(header_b, "mode")
  previous_pos = windex(header_b, "previous_pk")
  if mode_pos != None: weights[mode_pos] = 0
  if previous_pos != None: weights[previous_pos] = 0

  w = 100
  loser = INDEX_BY + []
  loser.reverse()               # I hate this
  for col in loser:
    pos = windex(header, col)
    if pos != None:
      weights[pos] = w
    w = w + w
  return weights

def find_best_matches(header1, header2, all_rows1, all_rows2,
                      pk_col, rows2_by_property):
  assert len(all_rows2) > 0
  correspondence = [windex(header1, column) for column in header2]
  print("Correspondence: %s" % correspondence, file=sys.stderr)
  positions = indexed_positions(header1)
  print("Indexed: %s" % positions, file=sys.stderr)
  weights = get_weights(header2, header1)
  print("Weights: %s" % weights, file=sys.stderr)
  pk_pos = windex(header1, pk_col)
  if pk_pos == None:
    print("No %s column in %s" % (pk_col, header1,),
          file=sys.stderr)
  assert pk_pos != None
  no_info = (-1, None)

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
        score = compute_score(row1, row2, correspondence, weights)
        best1_so_far = best_in_file1.get(key2, no_info)

        # Update best file2 match for row1
        (score2_so_far, key2_so_far) = best2_so_far
        if score > score2_so_far:
          best2_so_far = (score, key2)
        elif score == score2_so_far and key2 != key2_so_far:
          best2_so_far = (score, AMBIGUOUS)

        # Update best file1 match for row2
        (score1_so_far, key1_so_far) = best1_so_far
        if score > score1_so_far:
          best_in_file1[key2] = (score, key1)
        elif score == score1_so_far and key1 != key1_so_far:
          best_in_file1[key2] = (score, AMBIGUOUS)

    if best2_so_far != no_info:
      best_in_file2[key1] = best2_so_far

  print("%s properties" % prop_count, file=sys.stderr)
  if len(all_rows1) > 0:
    assert len(best_in_file2) > 0
  return (best_in_file1, best_in_file2)

def compute_score(row1, row2, correspondence, weights):
  s = 0
  for j in range(0, len(row2)):
    w = weights[j]
    if w != 0:
      if row2[j] != MISSING:
        i = correspondence[j]
        if i != None:
          if row1[i] != MISSING:
            s += w
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
            print("%s rows with property %s" % (LIMIT, property,),
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
  parser.add_argument('--state',
                      help='name of file to be annotated with matches to initial')
  parser.add_argument('--pk', default=None,
                      help='name of column containing primary key')
  # List of fields stored in database or graphdb should be an arg.
  args=parser.parse_args()
  with open(args.state, "r") as inport2:
    matchings(sys.stdin, inport2, args.pk, sys.stdout)
