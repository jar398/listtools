#!/usr/bin/env python3

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

# Most important first
INDEX_BY = ["EOLid", "source", "scientificName", "canonicalName"]

# Idea: distinguish between 'carry' and 'change' ??
GRAPHDB_COLUMNS = INDEX_BY +\
  ["taxonRank", "taxonomicStatus", "landmark"]

# -----

import sys, io, argparse, csv
from functools import reduce

MISSING = ''
AMBIGUOUS = '¿'

def matchings(inport1, inport2, outport):
  reader1 = csv.reader(inport1)
  header1 = next(reader1)
  reader2 = csv.reader(inport2)
  header2 = next(reader2)
  all_rows2 = read_rows(reader2, header2)
  rows2_by_property = get_rows_by_property(all_rows2, header2)
  (best_in_file1, best_in_file2, frontier) = \
    find_best_matches(reader1, header1, header2, all_rows2, rows2_by_property)
  emit_matches(best_in_file1, best_in_file2, header2, all_rows2, outport)
  emit_additions(best_in_file1, all_rows2, frontier, outport)

def emit_additions(best_in_file1, all_rows2, frontier, outport):
  writer = csv.writer(outport)
  for key2 in sorted(all_rows2.keys()):
    if not (key2 in best_in_file1):
      writer.writerow([frontier, "add"] + all_rows2[key2])
      frontier += 1

def emit_matches(best_in_file1, best_in_file2, header2, all_rows2, outport):
  def item_sort_key(item):
    (key1, (score, key2)) = item
    return (score, bleh(key1), bleh(key2))
  writer = csv.writer(outport)
  # Now emit the matches.
  writer.writerow(["aligned_id", "status"] + header2)
  blank = [MISSING for x in header2]
  for (key1, (score, key2)) in sorted(best_in_file2.items(),
                                    key=item_sort_key):
    mtype = match_type(key1, key2, score, best_in_file1)
    if score < 0: score = MISSING
    writer.writerow([key1, mtype] + all_rows2.get(key2, blank))
    # perhaps an explanation for the score and reason, as well??

def bleh(id):
  if id == None:
    return (-1, None)           # sorts before 0
  if (isinstance(id, int)):
    return (id, None)
  if (id.isdigit() or
      (id[0] == '-' and id[1:].isdigit())):
    return (int(id), id)
  else:
    return (1000000000, id)

def match_type(key1, key2, score, best_in_file1):
  assert key1 != AMBIGUOUS
  if key2 == AMBIGUOUS:
    # does not occur in 0.9/1.1
    return "ambiguous"
  else:
    best1 = best_in_file1.get(key2)
    if best1:
      (score3, id3) = best1
      if score != score3:
        return "remove"  # was "weak"
      elif id3 == AMBIGUOUS:
        return "remove"  # was "coambiguous"
      else:
        if key1 != id3:
          print("%s = %s != %s, score %s, back %s" % (key1, key2, id3, score, score3),
                file=sys.stderr)
          assert key1 == id3
        return "carry"     # the ones that matter
    else:
      return "remove"    # need better word

def indexed_positions(header):
  return [windex(header, col)
          for col in INDEX_BY
          if windex(header, col) != None]

def get_weights(header, header2):
  weights = [(1 if col in header2 else 0) for col in header]
  w = 100
  loser = INDEX_BY + []
  loser.reverse()               # I hate this
  for col in loser:
    pos = windex(header, col)
    if pos != None:
      weights[pos] = w
    w = w + w
  return weights

def find_best_matches(reader1, header1, header2, all_rows2, rows2_by_property):
  correspondence = [windex(header2, column) for column in header1]
  print("Correspondence: %s" % correspondence, file=sys.stderr)
  positions = indexed_positions(header1)
  print("Indexed: %s" % positions, file=sys.stderr)
  weights = get_weights(header1, header2)
  print("Weights: %s" % weights, file=sys.stderr)
  aligned_id_pos = windex(header1, "aligned_id")
  assert aligned_id_pos != None
  no_info = (-1, None)

  best_in_file2 = {}    # key1 -> (score, [record2, ...])  aligned id
  best_in_file1 = {}    # key2 -> (score, [record1, ...])  row number
  all_rows1 = {}
  count = 0
  frontier = 0
  for row1 in reader1:
    key1 = int(row1[aligned_id_pos])
    if key1 >= frontier: frontier = key1 + 1

    # The following check is also enforced by start.py... flush them here?
    assert not key1 in all_rows1; all_rows1[key1] = row1
    best2_so_far = no_info

    for prop in row_properties(row1, header1, positions):
      if count % 10000 == 0:
        print(count, file=sys.stderr)
      count += 1
      for rownum2 in rows2_by_property.get(prop, []):
        row2 = all_rows2[rownum2]
        score = compute_score(row1, row2, correspondence, weights)
        key2 = rownum2
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

    best_in_file2[key1] = best2_so_far

  return (best_in_file1, best_in_file2, frontier)

def compute_score(row1, row2, correspondence, weights):
  s = 0
  for i in range(0, len(row1)):
    w = weights[i]
    if w != 0:
      if row1[i] != MISSING:
        j = correspondence[i]
        if j != None:
          if row2[j] != MISSING:
            s += w
  return s

LIMIT=100

def get_rows_by_property(all_rows, header):
  positions = indexed_positions(header)
  by_property = {}
  entry_count = 0
  for rownum in all_rows:
    row = all_rows[rownum]
    for property in row_properties(row, header, positions):
      rownums = by_property.get(property)
      if rownums != None:
        if len(rownums) <= LIMIT:
          if len(rownums) == LIMIT:
            print("%s rows with property %s" % (LIMIT, property,),
                  file=sys.stderr)
          rownums.append(rownum)
          entry_count += 1
      else:
        by_property[property] = [rownum]
        entry_count += 1
  print("%s properties" % (len(by_property),),
        file=sys.stderr)
  return by_property

def read_rows(reader, header):
  all_rows = {}
  rownum = 1
  for row in reader:
    all_rows[rownum] = row
    rownum += 1
  print("%s rows" % (rownum,),
        file=sys.stderr)
  return all_rows

# Future: exclude really ephemeral properties like taxonID

def row_properties(row, header, positions):
  return [(header[i], row[i])
          for i in range(0, len(header))
          if (i in positions and
              row[i] != MISSING)]

def windex(header, fieldname):
  if fieldname in header:
    return header.index(fieldname)
  else:
    return None

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
    TBD.  Output is state with ids from input via alignment.
    """)
  parser.add_argument('--state',
                      help='name of file containing unaliged CSV input')
  args=parser.parse_args()
  with open(args.state, "r") as infile:
    matchings(sys.stdin, infile, sys.stdout)
