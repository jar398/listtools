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

import sys, io, argparse, csv
from functools import reduce

MISSING = ''
AMBIGUOUS = '¿'

# Most important first
INDEX_BY = ["EOLid", "source", "scientificName", "canonicalName"]

# Idea: distinguish between 'carry' and 'change' ??
GRAPHDB_COLUMNS = INDEX_BY +\
  ["taxonRank", "taxonomicStatus", "landmark"]

def matchings(inport1, inport2, outport):
  reader1 = csv.reader(inport1)
  header1 = next(reader1)
  reader2 = csv.reader(inport2)
  header2 = next(reader2)
  all_rows2 = read_rows(reader2, header2)
  rows2_by_property = get_records_by_property(all_rows2, header2)
  (best_in_file1, best_in_file2) = \
    find_best_matches(reader1, header1, header2, rows2_by_property)
  emit_matches(best_in_file1, best_in_file2, outport)
  emit_additions(best_in_file1, all_rows2, outport)

def emit_additions(best_in_file1, all_rows2, outport):
  writer = csv.writer(outport)
  for id2 in all_rows2.keys():
    if not (id2 in best_in_file1):
      writer.writerow([MISSING, id2, MISSING, "add"])

def emit_matches(best_in_file1, best_in_file2, outport):
  def item_sort_key(item):
    (id1, (score, id2)) = item
    return (score, bleh(id1), bleh(id2))
  writer = csv.writer(outport)
  # Now emit the matches.
  writer.writerow(["id1", "id2", "score", "note"])
  for (id1, (score, id2)) in sorted(best_in_file2.items(),
                                    key=item_sort_key):
    mtype = match_type(id1, id2, score, best_in_file1)
    if score < 0: score = MISSING
    writer.writerow([id1, id2, score, mtype])
    # perhaps an explanation for the score and reason, as well??

def bleh(id):
  if id == None:
    return (-1, None)           # sorts before 0
  if (id.isdigit() or
      (id[0] == '-' and id[1:].isdigit())):
    return (int(id), id)
  else:
    return (1000000000, id)

def match_type(id1, id2, score, best_in_file1):
  assert id1 != AMBIGUOUS
  if id2 == AMBIGUOUS:
    # does not occur in 0.9/1.1
    return "ambiguous"
  else:
    best1 = best_in_file1.get(id2)
    if best1:
      (score3, id3) = best1
      if score != score3:
        return "remove"  # was "weak"
      elif id3 == AMBIGUOUS:
        return "remove"  # was "coambiguous"
      else:
        if id1 != id3:
          print("%s = %s != %s, score %s, back %s" % (id1, id2, id3, score, score3),
                file=sys.stderr)
          assert id1 == id3
        return "carry"     # the ones that matter
    else:
      return "remove"    # need better word

def id_position(header):
  pos = windex(header, "id")
  # Kludge for when we're too lazy to run the file through start.sh
  if pos == None:
    pos = windex(header, "taxonID")
    if pos == None:
      pos = windex(header, "EOLid")
      if pos == None:
        print(header, file=sys.stderr)
        assert pos
  print("id column is %s" % header[pos], file=sys.stderr)
  return pos

def indexed_positions(header):
  return [windex(header, col)
          for col in INDEX_BY
          if windex(header, col) != None]

def get_weights(header, header2):
  weights = [(1 if col in header2 else 0) for col in header]
  weights[id_position(header)] = 0
  w = 100
  loser = INDEX_BY + []
  loser.reverse()               # I hate this
  for col in loser:
    pos = windex(header, col)
    if pos != None:
      weights[pos] = w
    w = w + w
  return weights

def find_best_matches(reader1, header1, header2, rows2_by_property):
  correspondence = [windex(header2, column) for column in header1]
  print("Correspondence: %s" % correspondence, file=sys.stderr)
  positions = indexed_positions(header1)
  print("Indexed: %s" % positions, file=sys.stderr)
  weights = get_weights(header1, header2)
  print("Weights: %s" % weights, file=sys.stderr)
  id1_pos = id_position(header1)
  id2_pos = id_position(header2)
  no_info = (-1, None)

  best_in_file2 = {}    # key1 -> (score, [record2, ...])
  best_in_file1 = {}    # key2 -> (score, [record1, ...])
  all_rows1 = {}
  count = 0
  for row1 in reader1:
    key1 = fix_id(row1, id1_pos)

    # The following check is also enforced by start.py... flush them here?
    assert not key1 in all_rows1; all_rows1[key1] = row1
    best2_so_far = no_info

    for prop in row_properties(row1, header1, positions):
      if count % 10000 == 0:
        print(count, file=sys.stderr)
      count += 1
      for row2 in rows2_by_property.get(prop, []):
        score = compute_score(row1, row2, correspondence, weights)
        key2 = row2[id2_pos]
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

  return (best_in_file1, best_in_file2)

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

def fix_id(row, id_pos):        # For better sort order
  id = row[id_pos]
  assert id != MISSING
  return id

LIMIT=100

def get_records_by_property(all_rows, header):
  id_pos = id_position(header)
  positions = indexed_positions(header)
  by_property = {}
  entry_count = 0
  for key in all_rows:
    row = all_rows[key]
    for property in row_properties(row, header, positions):
      records = by_property.get(property)
      if records != None:
        if len(records) <= LIMIT:
          if len(records) == LIMIT:
            print("%s records with property %s" % (LIMIT, property,),
                  file=sys.stderr)
          records.append(row)
          entry_count += 1
      else:
        by_property[property] = [row]
        entry_count += 1
  print("%s properties" % (len(by_property),),
        file=sys.stderr)
  return by_property

def read_rows(reader, header):
  id_pos = windex(header, "id")
  all_rows = {}
  for row in reader:
    key = fix_id(row, id_pos)
    assert key
    assert key != MISSING
    assert not key in all_rows
    all_rows[key] = row
  print("%s rows" % (len(all_rows),),
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
    TBD
    """)
  parser.add_argument('--file2',
                      help='name of file containing second CSV input')
  args=parser.parse_args()
  with open(args.file2, "r") as infile:
    matchings(sys.stdin, infile, sys.stdout)
