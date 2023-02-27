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

# 20 : 8 : 20 : 8

index_by_limit = 20    # max number of reverse entries to keep
silly = 16

unweights = None

default_index_by = ["scientificName",
                    "canonicalName", "canonicalStem", "tipe",
                    "managed_id"]
#default_discriminators = ["namePublishedInYear", "nomenclaturalStatus"]
default_discriminators = []

TOP = MISSING

# Row iterables -> row iterable

def match_records(a_reader, b_reader, pk_col="taxonID",
                  index_by=None,
                  discriminators=None):
  if index_by == None: index_by = default_index_by
  if discriminators == None: discriminators = default_discriminators
  assert index_by != None, index_by
  assert discriminators != None, discriminators

  a_table = ingest_csv(a_reader, pk_col)
  b_table = ingest_csv(b_reader, pk_col)
  # match_checklists returns a generator
  cop = match_checklists(a_table, b_table, pk_col, index_by, discriminators)
  yield from generate_match_report(cop, pk_col)

# Sum -> row iterable, suitable for writing to output CSV file.
# B is priority, so treat A matches as annotations on it
# Part of this module's API.
# coproduct arg is a generator, result is a generator.

def generate_match_report(coproduct, pk_col):
  yield [pk_col, "relationship", "match_id", "direction", "kind", "basis_of_match"]
  kind_counters = {}
  for (key1, ship, key2, direction, kind, basis_of_match) in coproduct:

    dk = (direction, kind)
    counter = kind_counters.get(dk)
    if not counter:
      counter = [0]
      kind_counters[dk] = counter
    counter[0] += 1

    if kind != "no match":
      ship_sym = rcc5_symbol(ship)
      yield (key1, ship_sym, key2, direction, kind, basis_of_match)

  for (dk, counter) in kind_counters.items():
    (direction, kind) = dk
    log("-- match: %s %s: %s" % (direction, kind, counter[0]))

# table * table -> sum (cf. delta.py ?)

def match_checklists(a_table, b_table, pk_col_arg, index_by, discriminators):
  global header1, header2, pk_col, pk_pos1, pk_pos2
  pk_col = pk_col_arg
  (header1, all_rows1) = a_table   # dict all_rows1: key -> row
  (header2, all_rows2) = b_table
  pk_pos1 = windex(header1, pk_col)
  assert pk_pos1 != None 
  pk_pos2 = windex(header2, pk_col)
  assert pk_pos2 != None

  assert len(index_by) + len(discriminators) < index_by_limit

  prepare_rows(header1, all_rows1)
  prepare_rows(header2, all_rows2)

  (best_results2, best_results1) = \
    find_best_matches(header1, header2, all_rows1, all_rows2, 
                      index_by, discriminators)

  # -----

  log("-- match: %s rows1, %s with matches" % (len(all_rows1), len(best_results2)))

  for (key1, row1) in all_rows1.items():
    (rows2, kind1, rel1, basis1) = \
      analyze_matches(row1, best_results2, best_results1, True)
    if len(rows2) >= WAD_SIZE:
      rows2 = (); kind1 = "too many matches"
    # Store the correspondence.  key2 may be None.
    key2 = '|'.join(map(lambda row2: row2[pk_pos2], rows2))
    if kind1 == "mutual":
      yield (key1, rel1, key2, "A<->B", kind1, basis1)
    else:
      yield (key1, rel1, key2, "A->B", kind1, basis1)

  log("# %s rows1, %s with matches" % (len(all_rows2), len(best_results1)))

  for (key2, row2) in all_rows2.items():
    (rows1, kind2, rel2, basis2) = \
      analyze_matches(row2, best_results1, best_results2, False)
    if len(rows1) >= WAD_SIZE:
      rows1 = (); kind2 = "too many matches"
    key1 = '|'.join(map(lambda row1: row1[pk_pos1], rows1))
    if kind2 == "mutual":
      pass
    else:
      yield (key1, rel2, key2, "B->A", kind2, basis2)

no_result = (-1, [])

# Returns (relationship, string, comment)
def analyze_matches(row_a, best_results_b, best_results_a, key_a_in_A):
  pk_pos_a = pk_pos1 if key_a_in_A else pk_pos2
  pk_pos_b = pk_pos2 if key_a_in_A else pk_pos1
  key_a = row_a[pk_pos_a]
  (score, rows_b) = best_results_b.get(key_a, no_result)
  rows_a = ()
  if True:
    def reciprocal(row_b):
      key_b = row_b[pk_pos_b]
      (_, rows_a) = best_results_a.get(key_b, no_result)
      return row_a in rows_a
    reciprocal_rows_b = [row_b for row_b in rows_b if reciprocal(row_b)]
    if len(reciprocal_rows_b) >= 1:
      rows_b = reciprocal_rows_b
      if len(reciprocal_rows_b) == 1:
        kind = "mutual"
      else:
        kind = "ambiguous"
    elif len(rows_b) > 0:
      kind = "no mutual matches"
    else: # len(rows_b) == 0:
      kind = "no match"
  else:
    if len(rows_b) == 1:
      row_b = rows_b[0]
      assert isinstance(row_b[0], str)
      pk_pos_b = pk_pos2 if key_a_in_A else pk_pos1
      key_b = row_b[pk_pos_b]
      (_, rows_a) = best_results_a.get(key_b, no_result)
      key_a = row_a[pk_pos_a]
      assert len(rows_a) > 0, (key_a, key_b)
      pk_pos_a = pk_pos1 if key_a_in_A else pk_pos2
      keys_a = map(lambda row:row[pk_pos_a], rows_a)
      if key_a in keys_a:
        if len(rows_a) == 1:
          kind = "mutual"
        else:
          kind = "in competition"
      else:
        kind = "lost competition"
    elif len(rows_b) > 1:
      kind = "ambiguous"
    else:
      kind = "no match"
  rel = EQ if kind == "mutual" else NOINFO
  return (rows_b, kind, rel, match_basis(score, unweights))

# For each row in each input (A/B), compute the set of all best
# (i.e. highest scoring) matches in the opposite input.
# Returns (key1 -> (score, [row2, ...]), key2 -> (score, [row1, ...]))

def find_best_matches(header1, header2, all_rows1, all_rows2, 
                      index_by, discriminators):
  global pk_pos1, pk_pos2, unweights # Passed in
  assert len(all_rows2) > 0
  corr_12 = correspondence(header1, header2)
  positions = get_positions(header1, index_by)
  (a_weights, unweights) = get_weights(header1, header2, index_by, discriminators)
  if False:
    log("# correspondence: %s" % (corr_12,))
    log("# indexed: %s" % positions)
    log("# A weights: %s" % a_weights)
    log("# unweights: %s" %
          ", ".join(["%s: %s" % (x, y) for (x, y) in unweights]))
  rows2_by_property = index_rows_by_property(all_rows2, header2, index_by)

  best_results2 = {}    # key1 -> (score, [row2, ...])
  best_results1 = {}    # key2 -> (score, [row1, ...])
  prop_count = 0
  for (key1, row1) in all_rows1.items():
    assert isinstance(row1[0], str)
    assert len(row1) == len(header1)

    # best_results2: key1 -> (score, [row2...])
    result2_so_far = best_results2.get(key1, None) or no_result

    # Compare matching rows in rows2 to row1.
    for prop in row_properties(row1, header1, positions):
      if prop_count > 0 and prop_count % 500000 == 0:
        log("# %s property values..." % prop_count)
      prop_count += 1

      for row2 in rows2_by_property.get(prop, ()):
        assert isinstance(row2[0], str)
        assert len(row2) == len(header2)
        key2 = row2[pk_pos2]
        # key2 -> (score, [row1 ...])
        result1_so_far = best_results1.get(key2, None) or no_result

        score = compute_score(row1, row2, corr_12, a_weights)
        if False:
          log("# compare %s to %s (%s) because %s" % (key1, key2, score, prop))

        # Does row2 improve on existing result for row1?
        (score2_so_far, rows2) = result2_so_far # keyed by row1
        if score > score2_so_far:
          result2_so_far = (score, [row2])
        elif score == score2_so_far and not row2 in rows2:
          if len(rows2) < WAD_SIZE:
            rows2.append(row2)

        # Does row1 improve on existing result for row2?
        (score1_so_far, rows1) = result1_so_far # keyed by row2
        if score > score1_so_far:
          result1_so_far = (score, [row1])
        elif score == score1_so_far and not row1 in rows1:
          if len(rows1) < WAD_SIZE:
            rows1.append(row1)

        # key2 -> (score, [row1...])
        if len(result1_so_far[1]) > 0:
          best_results1[key2] = result1_so_far

    # key1 -> (score, [row2...])
    if len(result2_so_far[1]) > 0:
      best_results2[key1] = result2_so_far

  #log("# match_records: indexed %s property values" % prop_count)
  if len(all_rows1) > 0 and len(all_rows2) > 0:
    assert len(best_results1) > 0
    assert len(best_results2) > 0
  return (best_results2,           # Keyed by key1
          best_results1)           # Keyed by key2

# Higher score = more similar

def compute_score(row1, row2, corr_12, a_weights):
  assert len(a_weights) == len(row1), (a_weights, row1)
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
          score += 2
  return score

def match_basis(score, unweights):
  for (weight, col) in unweights:
    if (score & weight) > 0: return col
  log("# score = %x unweight len = %x" % (score, len(list(unweights))))
  return None

# Initialize weights used in score computation.  One for each column
# in file A.  Returns (weights, unweights) where
# weights = powers of 2 or zero, indexed by A-columns
# unweights = (weight, name) for each indexed column

def get_weights(header_a, header_b, index_by, discriminators):
  all = index_by + discriminators
  a_weights = [0 for col in header_a]
  unweights = []
  for i in range(0, len(all)):  # 0 = more important
    col = all[i]
    pos_a = windex(header_a, col)   # Position in an A row
    pos_b = windex(header_b, col)   # Position in a B row
    if pos_a != None and pos_b != None:
      weight = 1 << (silly + len(all) - i)
      a_weights[pos_a] = weight
    if i < len(index_by):
      unweights.append((weight, col))
  return (a_weights, unweights)

LIMIT=10

def index_rows_by_property(all_rows, header, index_by):
  positions = get_positions(header, index_by)
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
  # log("# match_records: %s rows indexed by values" % (len(rows_by_property),))
  return rows_by_property

# Positions in header2 of columns to be indexed (properties), as a list

def get_positions(header, index_by):
  p = []
  for col in index_by:
    w = windex(header, col)
    if w != None:
      p.append(w)
    else:
      log("* No %s column present" % col)
  return p

# Generate all the (propertyname, value) pairs selected from row at
# the given positions.
# A property is a pair (propertyname, value).
# Future: exclude really ephemeral properties like taxonID

def row_properties(row, header, positions):
  for pos in positions:
    if row[pos] != MISSING and row[pos] != "NA":
      yield (header[pos], row[pos])

# Had to do this for the Norway/Sweden comparison...
# Want various ways of expressing 'accepted' to be = in the score computation.
# Turns out this matters...

def prepare_rows(header, dict): # dict: key -> row
  tstat_col = windex(header, "taxonomicStatus")
  nstat_col = windex(header, "nomenclaturalStatus")
  for (id, row) in dict.items():
    if tstat_col != None:
      tstat = row[tstat_col]
      if (tstat == MISSING or
          tstat.startswith("accepted") or
          tstat.startswith("valid")):
        row[tstat_col] = "accepted"
    if nstat_col != None:
      if row[nstat_col] == MISSING:
        row[tstat_col] = "available"

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

