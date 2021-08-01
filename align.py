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
INDEX_BY = ["EOLid", "source", "scientificName", "canonicalName"]

GRAPHDB_COLUMNS = INDEX_BY +\
  ["taxonRank", "taxonomicStatus", "landmark"]

# -----

import sys, io, argparse, csv
from functools import reduce

MISSING = ''
AMBIGUOUS = '¿'

def matchings(inport1, inport2, pk_col, outport):
  reader1 = csv.reader(inport1)
  header1 = next(reader1)
  reader2 = csv.reader(inport2)
  header2 = next(reader2)
  pk_pos1 = windex(header1, pk_col)
  pk_pos2 = windex(header2, pk_col)

  all_rows2 = read_rows(reader2, pk_pos2)
  rows2_by_property = get_rows_by_property(all_rows2, header2)
  (best_in_file1, best_in_file2, all_rows1) = \
    find_best_matches(reader1, header1, header2, all_rows2, pk_col, rows2_by_property)
  inj = injection(all_rows1, all_rows2, best_in_file1, best_in_file2,
                  windex(header1, pk_col), windex(header2, pk_col))

  def item_sort_key(item):
    (key1, (score, key2)) = item
    return (score, key1, key2)
  writer = csv.writer(outport)
  # Now emit the matches.
  writer.writerow(["status"] + header2)

  accepted_pos = windex(header2, "acceptedNameUsageID") if pk_col == "taxonID" else None
  parent_pos   = windex(header2, "parentNameUsageID")   if pk_col == "taxonID" else None
  print("Accepted in %s, parent in %s" % (accepted_pos, parent_pos), file=sys.stderr)

  # Preserve anything in file2 that was also in file1
  # List things in file1 not preserved for file2
  for (key1, (score, key2)) in sorted(best_in_file2.items(),
                                      key=item_sort_key):
    if score < 0: score = MISSING
    row2 = all_rows2.get(key2)
    if row2 == None:
      mtype = "remove"
      row2 = [MISSING for x in header2]
      row2[pk_pos2] = key1
    else:
      mtype = "carry"
      row2[pk_pos2] = inj[key2]
      if accepted_pos != None:
        row2[accepted_pos] = inj.get(row2[accepted_pos])
      if parent_pos != None:
        row2[parent_pos] = inj.get(row2[parent_pos])
    writer.writerow([mtype] + row2)
    # perhaps an explanation for the score and reason, as well??

  # Add new things from file2 that weren't in file1
  for key2 in sorted(all_rows2.keys()):
    if not (key2 in best_in_file1):
      row2 = all_rows2[key2]
      row2[pk_pos2] = inj[key2]
      # !!!!! also need to map parent and accepted.
      writer.writerow(["add"] + row2)

# Map from rows2 key (row number) to rows1 key (primary key)

def injection(rows1, rows2, best_in_file1, best_in_file2, pk_pos1, pk_pos2):
  def item_sort_key(item):
    (key1, (score, key2)) = item
    return (score, key1, key2)
  inj = {}
  # file1 is keyed by pk, file2 is keyed by row number...
  for (key1, (score, key2)) in sorted(best_in_file2.items(),
                                      key=item_sort_key):
    if check_match(key1, key2, score, best_in_file1):
      inj[key2] = key1
  for (key2, row2) in sorted(rows2.items()):
    if not key2 in inj:
      row2 = rows2[key2]
      pk2 = row2[pk_pos2]
      if pk2 == MISSING:
        print("No pk in column %s: %s" % (pk_pos2, row2,), file=sys.stderr)
      inj[key2] = choose_pk(pk2, rows1)
  return inj

def choose_pk(pk2, rows1):
  assert pk2 != MISSING
  if pk2 in rows1:
    # Already mapped, i.e. collision
    foo = pk2.rsplit(1)
    if len(foo) == 2 and foo[1].isdigit():
      stem = foo[0]
      spin = int(foo[1])
    else:
      stem = pk2
      spin = 1
    while True:
      candidate = "%s^%s" % (stem, spin)
      print("Considering %s" % (candidate,), file=sys.stderr)
      if candidate in rows1:
        spin += 1
      else:
        break
    return candidate
  else:
    # Use pk2 for primary key; it doesn't conflict
    return pk2

def check_match(key1, key2, score, best_in_file1):
  assert key1 != AMBIGUOUS
  if key2 == AMBIGUOUS:
    # does not occur in 0.9/1.1
    print("Tie: id %s -> multiple rows" % (key1,),
          file=sys.stderr)
    return False
  else:
    best1 = best_in_file1.get(key2)
    if best1:
      (score3, key3) = best1
      if score != score3:
        print("Id %s defeated by id %s for row %s" % (key1, key3, key2,),
              file=sys.stderr)
        return False  # was "weak"
      elif key3 == AMBIGUOUS:
        print("Id %s tied with other id(s) for row %s" % (key1, key2,),
              file=sys.stderr)
        return False  # was "coambiguous"
      elif score < 100:
        print("Id %s match to row %s is too weak to use (score %s)" % (key1, key2, score),
              file=sys.stderr)
        return False
      else:
        if key1 != key3:
          print("%s = %s != %s, score %s, back %s" % (key1, key2, key3, score, score3),
                file=sys.stderr)
          assert key1 == key3
        # maybe "change" as well?
        return True
    else:
      return False    # need better word


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

def find_best_matches(reader1, header1, header2, all_rows2, pk_col, rows2_by_property):
  correspondence = [windex(header2, column) for column in header1]
  print("Correspondence: %s" % correspondence, file=sys.stderr)
  positions = indexed_positions(header1)
  print("Indexed: %s" % positions, file=sys.stderr)
  weights = get_weights(header1, header2)
  print("Weights: %s" % weights, file=sys.stderr)
  pk_pos = windex(header1, pk_col)
  if pk_pos == None:
    print("No %s column in %s" % (pk_col, header1,),
          file=sys.stderr)
  assert pk_pos != None
  no_info = (-1, None)

  best_in_file2 = {}    # key1 -> (score, [record2, ...])  aligned id
  best_in_file1 = {}    # key2 -> (score, [record1, ...])  row number
  all_rows1 = {}
  count = 0
  for row1 in reader1:
    key1 = row1[pk_pos]
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

  return (best_in_file1, best_in_file2, all_rows1)

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

def read_rows(reader, pk_pos):
  all_rows = {}
  for row in reader:
    pk = row[pk_pos]
    assert pk != MISSING
    all_rows[pk] = row
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
    TBD.  Output is state with ids from input via alignment.
    """)
  parser.add_argument('--state',
                      help='name of file containing unaliged CSV input')
  parser.add_argument('--pk', default=None,
                      help='name of column containing primary key')
  args=parser.parse_args()
  with open(args.state, "r") as infile:
    matchings(sys.stdin, infile, args.pk, sys.stdout)
