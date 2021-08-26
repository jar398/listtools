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

def matchings(inport1, inport2, pk_col_arg, indexed, managed_arg, outport):
  global INDEX_BY, pk_pos1, pk_pos2, pk_col
  pk_col = pk_col_arg
  INDEX_BY = indexed.split(",")    # kludge

  (header1, all_rows1) = read_csv(inport1, pk_col)
  (header2, all_rows2) = read_csv(inport2, pk_col)

  managed = [man
             for man in managed_arg.split(",")
             if man in header1 or man in header2]

  pk_pos1 = windex(header1, pk_col)
  pk_pos2 = windex(header2, pk_col)
  assert pk_pos1 != None 
  assert pk_pos2 != None

  cop = compute_coproduct(header1, header2,
                          all_rows1, all_rows2)

  write_delta(cop, header1, header2, all_rows1, all_rows2, managed, outport)

def compute_coproduct(header1, header2, all_rows1, all_rows2):
  inl = {}
  inr = {}
  out = {}

  (best_rows_in_file2, best_rows_in_file1) = \
    find_best_matches(header1, header2, all_rows1, all_rows2, pk_col)
  print("# %s A rows with B match(es), %s B rows with A match(es)" %
        (len(best_rows_in_file2), len(best_rows_in_file1)),
        file=sys.stderr)

  def connect(key1, key2, remark):
    # Choose an id in the coproduct
    if key2:
      key3 = key2
    elif key1 in all_rows2:     # id collision
      row1 = all_rows1[key1]
      key3 = "%s$%s" % (key1, stable_hash(row1))
      print("-- id %s is used differently in the two inputs.\n" + \
            "--   %s will replace it for first input" % (key1, key3),
            file=sys.stderr)
    else:
      key3 = key1
    # Establish correspondences
    if key1:
      inl[key1] = key3
    if key2:
      inr[key2] = key3
    out[key3] = (key1, key2, remark)
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
          remark = ".mutual"
        else:
          remark = (".round trip failed: %s -> %s -> %s" %
                    (key1, row2[pk_pos2], back1[pk_pos1]))
      else:
        # Probably an ambiguity {row1, row1'} <-> row2
        remark = ".2 " + remark2
    else:
      # Probably an ambiguity row1 <-> {row2, row2'}
      remark = ".1 " + remark
    return (key2, remark)

  # -----

  for (key1, row1) in all_rows1.items():
    (key2, remark) = find_match(key1, best_rows_in_file2, best_rows_in_file1)
    # Store the correspondence
    connect(key1, key2, remark)

  for (key2, row2) in all_rows2.items():
    if not key2 in inr:
      (key1, remark) = find_match(key2, best_rows_in_file1, best_rows_in_file2)
      if key1:
        # shouldn't happen
        remark = "surprising match: %s <-> %s" % (key2, key1)
      # Addition - why did it fail?
      # Unique match means outcompeted, probably?
      connect(None, key2, remark)

  # Print stats on outcome
  aonly = bonly = matched = 0
  for (key3, (key1, key2, remark)) in out.items():
    if key1 != None and key2 != None:
      matched += 1
    elif key1 == None:
      bonly += 1
    else:
      aonly += 1
  print("# coproduct: %s A only, %s A and B, %s B only" %
        (aonly, matched, bonly),
        file=sys.stderr)

  return (inl, inr, out)

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

# Write coproduct in the form of an EOL-ish delta

def write_delta(cop, header1, header2, all_rows1, all_rows2, managed, outport):
  global pk_col

  header3 = ['mode', 'new_pk', 'remark'] + managed
  mode_pos3 = 0
  new_pk_pos3 = 1
  remark_pos3 = 2
  # new_pk_pos3 will be the primary key within the delta itself
  old_pk_pos3 = 3 + windex(managed, pk_col)

  corr_13 = correspondence(header1, header3)
  corr_23 = correspondence(header2, header3)

  managed_positions2 = [windex(header2, name) for name in managed]

  (inl, inr, out) = cop

  writer = csv.writer(outport)

  def write_row(mode, row, old, new, co, remark):
    if new: assert co == new
    row[mode_pos3] = mode
    row[old_pk_pos3] = old or MISSING
    row[new_pk_pos3] = co
    row[remark_pos3] = remark
    writer.writerow(row)

  writer.writerow(header3)

  corr_12 = correspondence(header1, header2)
  add_count = 0
  carry_count = 0
  update_count = 0
  remove_count = 0
  stats = [[0, 0, 0, ([], [], [])] for col in header2]

  for (key3, (key1, key2, remark)) in out.items():
    if key1 != None and key2 != None:
      # Carry or difference in field values
      row1 = all_rows1[key1]
      row2 = all_rows2[key2]
      if analyze_changes(row1, row2, managed_positions2, corr_12, stats):
        write_row("update", apply_correspondence(corr_23, row2),
                  key1, key2, key3, remark)
        update_count += 1
      else:
        # Carry - no output
        carry_count += 1
    elif key1 != None:
      row1 = all_rows1[key1]
      write_row("remove", apply_correspondence(corr_13, row1),
                key1, key2, key3, remark)
      remove_count += 1
    else:
      row2 = all_rows2[key2]
      write_row("add", apply_correspondence(corr_23, row2),
                key1, key2, key3, remark)
      add_count += 1

  print("-- delta: %s carries, %s additions, %s removals, %s updates" %
        (carry_count, add_count, remove_count, update_count,),
        file=sys.stderr)
  for j in range(0, len(header2)):
    (a, c, d, (qs, cs, ds)) = stats[j]
    if a > 0:
      x = [row2[pk_pos2] for row2 in qs]
      print("--   %s: %s set %s" % (header2[j], a, x),
            file=sys.stderr)
    if c > 0:
      x = [row1[pk_pos1] for row1 in cs]
      print("--   %s: %s modified %s" % (header2[j], c, x),
            file=sys.stderr)
    if d > 0:
      x = [row1[pk_pos1] for row1 in ds]
      print("--   %s: %s cleared %s" % (header2[j], d, x),
            file=sys.stderr)

# for readability
SAMPLES = 3

# managed_positions2 = positions in header2 of the change-managed columns.
# Side effect: increment counters per column of added, updated, removed rows

def analyze_changes(row1, row2, managed_positions2, corr_12, stats):
  for pos2 in managed_positions2:
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
  indexed="taxonID,scientificName,canonicalName"
  parser.add_argument('--index',
                      default=indexed,
                      help='names of columns to match on')
  managed=indexed+",taxonRank,taxonomicStatus,nomenclaturalStatus,source,datasetID"
  parser.add_argument('--manage',
                      default=managed,
                      help='names of columns under version control')
  # List of fields stored in database or graphdb should be an arg.
  args=parser.parse_args()
  with open(args.target, "r") as inport2:
    matchings(sys.stdin, inport2, args.pk, args.index, args.manage, sys.stdout)
