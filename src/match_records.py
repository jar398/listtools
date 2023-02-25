#!/usr/bin/env python3

"""
This code was written prior to property.py and align.py, and would
probably be more pleasant to read if it used property.py.

Here's what we're trying to do, as an n^2 method:

  For (record 1, record 2) pair, compute match score.
  Consider all candidate matches (x1, x2) where 
    either x1 = record 1 or x2 = record 2.
  Designate (record 1, record 2) a match if it has the highest score
    among all of these candidates.
  (I.e. if record 1 = unique best match to record 2 AND v.v.)

With some indexing, we can do it in approximately linear time.

A row is represented according to what the CSV reader returns (i.e. a
list of strings).

pk_col is the name of the primary key column (typically "taxonID"),
while pk_pos is the position in the row list of the primary key.

"""

# Largest number of alternative matches to provide
WAD_SIZE = 4

# -----

import sys, io, argparse, csv
from functools import reduce
from util import ingest_csv, windex, MISSING, \
                 correspondence, precolumn, stable_hash, log
from rcc5 import *

MUTUAL = EQ       # relationship to use for mutual match
LOSER = NOINFO    # relationship to use for runners-up

# 20 : 8 : 20 : 8

index_by_limit = 20    # max number of reverse entries to keep
silly = 16

unweights = None

default_index_by = ["scientificName", "tipe",
                    "canonicalName", "canonicalStem",
                    "managed_id"]

# Row iterables -> row iterable

def match_records(a_reader, b_reader, pk_col="taxonID", index_by=default_index_by):
  a_table = ingest_csv(a_reader, pk_col)
  b_table = ingest_csv(b_reader, pk_col)
  cop = compute_sum(a_table, b_table, pk_col, index_by)
  return generate_sum(cop, pk_col)

# Sum -> row iterable, suitable for writing to output CSV file.
# B is priority, so treat A matches as annotations on it
# Part of this module's API.
# coproduct arg is a generator, result is a generator.

def generate_sum(coproduct, pk_col):
  yield ["match_id", "relationship", pk_col, "basis_of_match"]
  aonly = bonly = matched = ambig = 0
  for (key1, ship, key2, remark) in coproduct:

    # Print stats on outcome
    if ship == MUTUAL:
      matched += 1
    elif key1 == None:
      bonly += 1
    else:
      aonly += 1
    if "ambiguous" in remark:
      ambig += 1

    if ship != NOINFO or remark:
      ship_sym = rcc5_symbol(ship)
      yield (key1, ship_sym, key2, remark)

  print("-- match_records: %s matches, %s in A unmatched, %s in B unmatched" %
        (matched, aonly, bonly),
        file=sys.stderr)
  print("-- match_records: %s ambiguous" %
        (ambig),
        file=sys.stderr)

# table * table -> sum (cf. delta.py ?)

def compute_sum(a_table, b_table, pk_col_arg, index_by):
  global INDEX_BY, pk_col, pk_pos1, pk_pos2
  pk_col = pk_col_arg
  assert len(index_by) < index_by_limit
  INDEX_BY = index_by    # kludge

  (header1, all_rows1) = a_table
  (header2, all_rows2) = b_table
  prepare_rows(header1, all_rows1)
  prepare_rows(header2, all_rows2)

  pk_pos1 = windex(header1, pk_col)
  pk_pos2 = windex(header2, pk_col)
  assert pk_pos1 != None 
  assert pk_pos2 != None

  (best_rows_in_file2, best_rows_in_file1) = \
    find_best_matches(header1, header2, all_rows1, all_rows2, pk_col)

  # Returns (relationship, string, comment)
  def find_match(key1, best_rows_in_file2, best_rows_in_file1,
                 pk_pos1, pk_pos2, key1_in_A):
    (row2, rows2, remark2) = check_match(key1, best_rows_in_file2, pk_pos2, key1_in_A)
    # remark = MISSING means no matches
    if row2 != None:            # len(rows2) == 1
      candidate = row2[pk_pos2]
      (back1, rows1, remark1) = check_match(candidate, best_rows_in_file1, pk_pos1, not key1_in_A)
      assert remark1
      if back1 != None:         # len(rows1) == 1
        if back1[pk_pos1] == key1:
          # Mutual match!
          key2 = candidate
          # remark2 and remark1 give fields
          if remark2 == remark1:
            remark = remark2
          else:
            # Shouldn't happen I think ??
            remark = (key2, "%s|%s" % (remark2, remark1))
          answer = (MUTUAL, key2, remark)
        else:
          b1 = back1[pk_pos1]
          # wish there was a name here so we can see what's going on
          # should be sensitive to key1_in_A
          if key1_in_A:
            answer = (LOSER, row2[pk_pos2],
                      ("match not mutual: %s -> %s" %
                       (row2[pk_pos2], b1)))
          else:
            answer = (LOSER, row2[pk_pos2],
                      ("match not mutual: %s <- %s" %
                       (b1, row2[pk_pos2])))
      else:
        # Probably an ambiguity {row1, row1'} <-> row2
        answer = (NOINFO, TOP, remark1)
    else:
      # No match in B, or an ambiguity row1 <-> {row2, row2'}
      answer = (NOINFO, TOP, remark2)
    return answer

  TOP = MISSING

  # -----

  seen2 = {}                    # injection B -> A+B
  punt = 0

  for (key1, row1) in all_rows1.items():
    (rel2, key2, remark) = find_match(key1, best_rows_in_file2,
                                      best_rows_in_file1, pk_pos1, pk_pos2,
                                      True)
    # Store the correspondence.  key2 may be None.
    yield (key1, rel2, key2, remark)
    if key2 != None:
      seen2[key2] = True

  copunt = 0

  for (key2, row2) in all_rows2.items():
    if not key2 in seen2:
      (rel1, key1, remark) = find_match(key2, best_rows_in_file1,
                                        best_rows_in_file2, pk_pos2, pk_pos1,
                                        False)
      if key1:
        # shouldn't happen
        remark = "! surprising match: %s <-> %s" % (key2, key1)
      # Addition - why did it fail?
      # Unique match means outcompeted, probably?
      rel2 = reverse_relationship(rel1) # <= to >=
      yield (key1, rel2, key2, remark)

  if punt + copunt > 0:
    print("-- Senior synonyms: %s in A, %s in B" % (punt, copunt),
          file=sys.stderr)

# Makes sure the match is mutual and unique.
# Returns a row.

def check_match(key1, best_rows_in_file2, pk_pos2, key1_in_A):
  assert len(unweights) > 0
  best2 = best_rows_in_file2.get(key1)
  if best2:
    (score2, rows2) = best2
    reason = match_reason(score2, unweights)
    if reason != None:
      if len(rows2) == 1:
        assert rows2[0]
        return (rows2[0], rows2, reason)
      elif len(rows2) < WAD_SIZE:
        keys2 = [row2[pk_pos2] for row2 in rows2]
        alts = ';'.join(keys2)
        if key1_in_A:
          return (None, rows2,
                        ("ambiguous (%s): -> %s" %
                         (reason, alts)))
        else:
          return (None, rows2,
                        ("coambiguous (%s): <- %s" %
                         (reason, alts)))
      else:
        return (None, (), "too many matches (%s)" % reason)
    else:
      return (None, (), "weak matches only (%x)" % score2)
  else:
    return (None, (), MISSING)

# For each row in each input (A/B), compute the set of all best
# (i.e. highest scoring) matches in the opposite input.

def find_best_matches(header1, header2, all_rows1, all_rows2, pk_col):
  global pk_pos1, pk_pos2, unweights # Passed in
  assert len(all_rows2) > 0
  corr_12 = correspondence(header1, header2)
  positions = indexed_positions(header1, INDEX_BY)
  (a_weights, unweights) = get_weights(header1, header2, INDEX_BY)
  if False:
    print("# correspondence: %s" % (corr_12,), file=sys.stderr)
    print("# indexed: %s" % positions, file=sys.stderr)
    print("# A weights: %s" % a_weights, file=sys.stderr)
    print("# unweights: %s" %
          ", ".join(["%s: %s" % (x, y) for (x, y) in unweights]),
          file=sys.stderr)
  rows2_by_property = index_rows_by_property(all_rows2, header2)
  no_info = (-1, [])

  best_rows_in_file2 = {}    # key1 -> (score, rows2)
  best_rows_in_file1 = {}    # key2 -> (score, rows1)
  prop_count = 0
  for (key1, row1) in all_rows1.items():
    assert len(row1) == len(header1)
    # The following check is also enforced by start.py... flush them here?
    best2_so_far = no_info
    best_rows_so_far2 = no_info

    for prop in row_properties(row1, header1, positions):
      if prop_count > 0 and prop_count % 500000 == 0:
        print("# %s property values..." % prop_count, file=sys.stderr)
      prop_count += 1
      for row2 in rows2_by_property.get(prop, []):
        assert len(row2) == len(header2)
        key2 = row2[pk_pos2]
        score = compute_score(row1, row2, corr_12, a_weights)
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

  print("# match_records: indexed %s property values" % prop_count, file=sys.stderr)
  if len(all_rows1) > 0 and len(all_rows2) > 0:
    assert len(best_rows_in_file1) > 0
    assert len(best_rows_in_file2) > 0
  return (best_rows_in_file2, best_rows_in_file1)

# Positions in header2 of columns to be indexed (properties)

def indexed_positions(header, index_by):
  p = []
  for col in index_by:
    w = windex(header, col)
    if w != None:
      p.append(w)
    else:
      print("* No %s column present" % col, file=sys.stderr)
  return p

# Higher score = more similar

def compute_score(row1, row2, corr_12, a_weights):
  score = 0
  for j in range(0, len(row2)):
    i = precolumn(corr_12, j)
    if i != None:
      v1 = row1[i]
      v2 = row2[j]
      if v1 != MISSING and v2 != MISSING and v1 == v2:
        weight = a_weights[i]
        if weight != 0:
          score += weight
        else:
          # give a tiny bit of weight for unweighted columns
          # this could give random results...
          score += 1
  return score

def match_reason(score, unweights):
  for (weight, col) in unweights:
    if (score & weight) > 0: return col
  print("# score = %x unweight len = %x" % (score, len(list(unweights))),
        file=sys.stderr)
  return None

# Initialize weights used in score computation.  One for each column
# in file A.  Returns (weights, unweights) where
# weights = powers of 2 or zero, indexed by A-columns
# unweights = (weight, name) for each indexed column

def get_weights(header_b, header_a, index_by):
  a_weights = [0 for col in header_a]
  unweights = []
  for i in range(0, len(index_by)):
    col = index_by[i]
    pos = windex(header_a, col)   # Position in an A row
    if pos != None:
      weight = 1 << (silly + len(index_by) - i)
      a_weights[pos] = weight
    unweights.append((weight, col))
  return (a_weights, unweights)

LIMIT=10

def index_rows_by_property(all_rows, header):
  positions = indexed_positions(header, INDEX_BY)
  rows_by_property = {}
  entry_count = 0
  for (key, row) in all_rows.items():
    assert len(row) == len(header)
    for property in row_properties(row, header, positions):
      rows = rows_by_property.get(property)
      if rows != None:
        if len(rows) <= LIMIT:
          if len(rows) == LIMIT:
            log("# %s+ rows with property %s" % (LIMIT, property,))
            log("#   e.g. %s" % (rows[0],))
          rows.append(row)
          entry_count += 1
      else:
        rows_by_property[property] = [row]
        entry_count += 1
  # print("# match_records: %s rows indexed by values" % (len(rows_by_property),),
  #       file=sys.stderr)
  return rows_by_property

# Future: exclude really ephemeral properties like taxonID

def row_properties(row, header, positions):
  return [(header[i], row[i])
          for i in range(0, len(header))
          if (i in positions and
              row[i] != MISSING)]

# Had to do this for the Norway/Sweden comparison.
# Want various ways of expressing 'accepted' to be = in the score computation.
# Turns out this matters...

def prepare_rows(header, table):
  tstat_col = windex(header, "taxonomicStatus")
  n = 0
  if tstat_col != None:
    for (id, row) in table.items():
      tstat = row[tstat_col]
      if (tstat == MISSING or
          tstat.startswith("accepted") or
          tstat.startswith("valid")):
        row[tstat_col] = "accepted"

# Test
def test1():
  inport1 = io.StringIO(u"taxonID,bar\n1,2")
  inport2 = io.StringIO(u"taxonID,bar\n1,3")
  main(inport1, inport2, "taxonID", ["taxonID", "bar"], sys.stdout)

def test():
  inport1 = io.StringIO(u"taxonID,bar\n1,dog\n2,cat")
  inport2 = io.StringIO(u"taxonID,bar\n91,cat\n93,pig")
  main(inport1, inport2, "taxonID", ["taxonID", "bar"], sys.stdout)

# -----------------------------------------------------------------------------
# Command line interface

def main(inport1, inport2, pk_col, index_by, outport):
  gen = match_records(csv.reader(inport1), csv.reader(inport2),
                      pk_col=pk_col, index_by=index_by)
  writer = csv.writer(outport)
  for row in gen: writer.writerow(row)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    Standard input is the 'left' or 'A' usage records.
    Standard output is the record-level sum A+B of the two inputs,
    with rows of the form [id-in-A+B, id-in-A, id-in-B, remark].
    """)
  parser.add_argument('--A',
                      help="name of file containing 'left' or 'A' usages")
  parser.add_argument('--B',
                      help="name of file containing 'right' or 'B' usages")
  parser.add_argument('--pk',
                      default="taxonID",
                      help='name of column containing primary key')
  # Order is important
  parser.add_argument('--index',
                      default=None,
                      help='names of columns to match on')
  args=parser.parse_args()
  with open(args.A, "r") as inport1:
    with open(args.B, "r") as inport2:
      main(inport1, inport2, args.pk,
           args.index.split(",") if args.index else None,
           sys.stdout)

