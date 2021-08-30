#!/usr/bin/env python3

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
                 correspondence, precolumn

def main(inport1, inport2, pk_col, indexed, outport):
  index_by = indexed.split(",")    # kludge
  (header1, all_rows1) = read_csv(inport1, pk_col)
  (header2, all_rows2) = read_csv(inport2, pk_col)

  cop = compute_coproduct(header1, header2, all_rows1, all_rows2,
                          pk_col, index_by)

  write_coproduct(cop, pk_col, outport)

def write_coproduct(cop, pk_col, outport):
  writer = csv.writer(outport)
  writer.writerow([pk_col, pk_col + "_A", pk_col + "_B", "remark"])
  for (key3, (key1, key2, remark)) in cop.items():
    writer.writerow([key3, key1, key2, remark])

def compute_coproduct(header1, header2, all_rows1, all_rows2, pk_col_arg, index_by):
  global INDEX_BY, pk_col, pk_pos1, pk_pos2
  pk_col = pk_col_arg
  INDEX_BY = index_by    # kludge

  pk_pos1 = windex(header1, pk_col)
  pk_pos2 = windex(header2, pk_col)
  assert pk_pos1 != None 
  assert pk_pos2 != None

  seen2 = {}                    # injection B -> A+B
  cop = {}

  (best_rows_in_file2, best_rows_in_file1) = \
    find_best_matches(header1, header2, all_rows1, all_rows2, pk_col)
  print("# %s A rows with B match(es), %s B rows with A match(es)" %
        (len(best_rows_in_file2), len(best_rows_in_file1)),
        file=sys.stderr)

  def connect(key1, key2, remark):
    # Choose an id in the coproduct
    if key2:
      key3 = key2
      seen2[key2] = key3
    elif key1 in all_rows2:     # id collision
      row1 = all_rows1[key1]
      key3 = "%s$%s" % (key1, stable_hash(row1))
      print("-- id %s is used differently in the two inputs.\n" + \
            "--   %s will replace it for the A input" % (key1, key3),
            file=sys.stderr)
    else:
      key3 = key1
    # Establish correspondences
    cop[key3] = (key1, key2, remark)
    if not remark.startswith("."):
      print("# %s" % (remark,), file=sys.stderr)

  def find_match(key1, best_rows_in_file2, best_rows_in_file1):
    key2 = None
    (row2, remark) = check_match(key1, best_rows_in_file2)
    if row2 != None:
      candidate = row2[pk_pos2]
      (back1, remark2) = check_match(candidate, best_rows_in_file1)
      if back1 != None:
        if back1[pk_pos1] == key1:
          # Mutual match!
          key2 = candidate
          remark = MISSING
        else:
          remark = (".round trip failed: %s -> %s -> %s" %
                    (key1, row2[pk_pos2], back1[pk_pos1]))
      else:
        # Probably an ambiguity {row1, row1'} <-> row2
        remark = ".2 " + remark2
    else:
      # No match in B, or an ambiguity row1 <-> {row2, row2'}
      remark = ".1 " + remark
    return (key2, remark)

  # -----

  for (key1, row1) in all_rows1.items():
    (key2, remark) = find_match(key1, best_rows_in_file2, best_rows_in_file1)
    # Store the correspondence
    connect(key1, key2, remark)

  for (key2, row2) in all_rows2.items():
    if not key2 in seen2:
      (key1, remark) = find_match(key2, best_rows_in_file1, best_rows_in_file2)
      if key1:
        # shouldn't happen
        remark = "surprising match: %s <-> %s" % (key2, key1)
      # Addition - why did it fail?
      # Unique match means outcompeted, probably?
      connect(None, key2, remark)

  # Print stats on outcome
  aonly = bonly = matched = 0
  for (key3, (key1, key2, remark)) in cop.items():
    if key1 != None and key2 != None:
      matched += 1
    elif key1 == None:
      bonly += 1
    else:
      aonly += 1
  print("-- matching records: %s matched, %s in A unmatched, %s in B unmatched" %
        (matched, aonly, bonly),
        file=sys.stderr)

  return cop

WAD_SIZE = 4

def check_match(key1, best_rows_in_file2):
  global pk_pos1, pk_pos2
  best2 = best_rows_in_file2.get(key1)
  if best2:
    (score2, rows2) = best2
    if len(rows2) == 1:
      return (rows2[0], ".unique")    # unique match
    elif len(rows2) < WAD_SIZE:
      keys2 = [row2[pk_pos2] for row2 in rows2]
      return (None, ("ambiguous: %s -> %s (score %s)" %
                     (key1, keys2, score2)))
    else:
      return (None, ".highly ambiguous")
  else:
    return (None, ".no matches")

def find_best_matches(header1, header2, all_rows1, all_rows2, pk_col):
  global pk_pos1, pk_pos2
  assert len(all_rows2) > 0
  corr_12 = correspondence(header1, header2)
  print("# correspondence: %s" % (corr_12,), file=sys.stderr)
  positions = indexed_positions(header1, INDEX_BY)
  print("# indexed: %s" % positions, file=sys.stderr)
  weights = get_weights(header1, header2, INDEX_BY)    # parallel to header2
  print("# weights: %s" % weights, file=sys.stderr)
  rows2_by_property = index_rows_by_property(all_rows2, header2)
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
        print("# %s values" % prop_count, file=sys.stderr)
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

  print("# delta: indexed %s values" % prop_count, file=sys.stderr)
  if len(all_rows1) > 0 and len(all_rows2) > 0:
    assert len(best_rows_in_file1) > 0
    assert len(best_rows_in_file2) > 0
  return (best_rows_in_file2, best_rows_in_file1)

# Positions in header2 of columns to be indexed (properties)

def indexed_positions(header, index_by):
  return [windex(header, col)
          for col in index_by
          if windex(header, col) != None]

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
            print("# %s+ rows with property %s" % (LIMIT, property,),
                  file=sys.stderr)
          rows.append(row)
          entry_count += 1
      else:
        rows_by_property[property] = [row]
        entry_count += 1
  print("# %s rows indexed by values" % (len(rows_by_property),),
        file=sys.stderr)
  return rows_by_property

# Future: exclude really ephemeral properties like taxonID

def row_properties(row, header, positions):
  return [(header[i], row[i])
          for i in range(0, len(header))
          if (i in positions and
              row[i] != MISSING)]

# Inverse of write_coproduct

def read_coproduct(cofile, pk_col):
  reader = csv.reader(cofile)
  header = next(reader)
  pk_pos = windex(header, pk_col)
  pk_a_pos = windex(header, pk_col + "_A")
  pk_b_pos = windex(header, pk_col + "_B")
  remark_pos = windex(header, "remark")
  assert pk_pos != None and pk_a_pos != None and pk_b_pos != None
  assert remark_pos != None
  cop = {}
  for row in reader:
    a = row[pk_a_pos]
    if a == MISSING: a = None
    a = row[pk_b_pos]
    if b == MISSING: b = None
    cop[row[pk_pos]] = (a, b, row[remark+pos])
  return cop

# Test
def test1():
  inport1 = io.StringIO(u"taxonID,bar\n1,2")
  inport2 = io.StringIO(u"taxonID,bar\n1,3")
  main(inport1, inport2, "taxonID", "taxonID,bar", sys.stdout)

def test():
  inport1 = io.StringIO(u"taxonID,bar\n1,dog\n2,cat")
  inport2 = io.StringIO(u"taxonID,bar\n91,cat\n93,pig")
  main(inport1, inport2, "taxonID", "taxonID,bar", sys.stdout)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Standard input is the 'left' or 'A' usage records.
    Standard output is the record-level coproduct A+B of the two inputs,
    with rows of the form [id-in-A+B, id-in-A, id-in-B, remark].
    """)
  parser.add_argument('--target',
                      help="name of file containing 'right' or 'B' usages")
  parser.add_argument('--pk',
                      default="taxonID",
                      help='name of column containing primary key')
  # Order is important
  indexed="taxonID,scientificName,canonicalName"
  parser.add_argument('--index',
                      default=indexed,
                      help='names of columns to match on')
  args=parser.parse_args()
  with open(args.target, "r") as inport2:
    main(sys.stdin, inport2, args.pk, args.index, sys.stdout)