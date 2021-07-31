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


import sys, io, csv, argparse

MISSING = ''
AMBIGUOUS = '¿'

INDEX_BY = ["EOLid", "source", "scientificName", "canonicalName"]

def matchings(inport1, inport2, outport):
  reader1 = csv.reader(inport1)
  header1 = next(reader1)
  reader2 = csv.reader(inport2)
  header2 = next(reader2)
  rows2_by_property = get_records_by_property(reader2, header2)
  (best_in_file1, best_in_file2) = \
    find_best_matches(reader1, header1, header2, rows2_by_property)
  emit_matches(best_in_file1, best_in_file2, outport)

def emit_matches(best_in_file1, best_in_file2, outport):
  writer = csv.writer(outport)
  # Now emit the matches.
  writer.writerow(["key1", "key2", "score"])
  for (key1, (score, key2)) in best_in_file2.items():
    assert isinstance(key1, str)
    assert isinstance(key2, str)
    assert isinstance(score, int)
    assert key1 != AMBIGUOUS
    if key2 != AMBIGUOUS:
      rebest = best_in_file1.get(key2)
      if rebest:
        (score3, key3) = rebest
        assert isinstance(key3, str)
        assert isinstance(score3, int)
        if score == score3 and key3 != AMBIGUOUS:
          if key1 != key3:
            print("%s = %s != %s, score %s, back %s" % (key1, key2, key3, score, score3),
                  file=sys.stderr)
          assert key1 == key3
          writer.writerow([key1, key2, score])

def primary_key_position(header):
  key = windex(header, "taxonID")
  if key == None: key = windex(header, "EOLid")
  if key == None:
    print(header, file=sys.stderr)
    assert key
  return key

def indexed_positions(header):
  return [windex(header, col)
          for col in INDEX_BY
          if windex(header, col) != None]

def find_best_matches(reader1, header1, header2, rows2_by_property):
  correspondence = [windex(header2, column) for column in header1]
  print("Correspondence: %s" % correspondence, file=sys.stderr)
  positions = indexed_positions(header1)
  print("Indexed: %s" % positions, file=sys.stderr)
  pk1_pos = primary_key_position(header1)
  pk2_pos = primary_key_position(header2)
  no_info = (-1, None)

  best_in_file2 = {}    # key1 -> (score, [record2, ...])
  best_in_file1 = {}    # key2 -> (score, [record1, ...])
  seen_keys = {}
  count = 0
  for row1 in reader1:
    key1 = row1[pk1_pos]
    assert key1 != MISSING
    assert not key1 in seen_keys; seen_keys[key1] = True
    best2_so_far = no_info

    for prop in row_properties(row1, header1, positions):
      if count % 10000 == 0:
        print(count, file=sys.stderr)
      count += 1
      for row2 in rows2_by_property.get(prop, []):
        score = compute_score(row1, row2, correspondence)
        key2 = row2[pk2_pos]
        row2_primary_key = row2[pk2_pos]
        best1_so_far = best_in_file1.get(row2_primary_key, no_info)

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

    (_, key2) = best2_so_far
    if key2 != None and key2 != AMBIGUOUS:
      best_in_file2[key1] = best2_so_far

  return (best_in_file1, best_in_file2)

def compute_score(row1, row2, correspondence):
  s = 0
  for i in range(0, len(row1)):
    if row1[i] != MISSING:
      j = correspondence[i]
      if j != None:
        if row2[j] != MISSING:
          s += 1
  return s

LIMIT=100

def get_records_by_property(reader, header):
  by_property = {}
  seen_keys = {}
  pos = primary_key_position(header)
  positions = indexed_positions(header)
  entry_count = 0
  for row in reader:
    key = row[pos]
    assert key != MISSING
    assert not key in seen_keys
    seen_keys[key] = True
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
  print("%s keys, %s properties" % (len(seen_keys), len(by_property)),
        file=sys.stderr)
  return by_property

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
