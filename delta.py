#!/usr/bin/env python3

# Prepare CQL queries that can accomplish the patch.
# This really ought to be done in Ruby.


import sys, os, csv, argparse, pathlib
from util import windex

# Later: do it only for records provided by a particular resource (e.g. traits)

def write_cql_files(dest, pk_col):
  node_type = "Usage"
  prepare_add_cql(node_type, dest)
  prepare_remove_cql(node_type, pk_col, dest)
  prepare_update_cql(node_type, pk_col, dest)

def get_header(fname):
  # Just get the header
  with open(fname, "r") as infile:
    reader = csv.reader(infile)
    header = next(reader)
  return header

def prepare_add_cql(node_type, dest):
  # Now prepare the CQL file
  csv_name = os.path.join(dest, "add.csv")
  cql_name = os.path.join(dest, "add.cql")
  props = get_header(csv_name)
  pk_col = windex(props, "previous_pk")
  assert pk_col != None
  del props[pk_col]
  settings = ",\n  ".join(["%s: row.%s" % zz for zz in zip(props, props)])
  csv_url = "file:///%s" % (csv_name,)
  with open(cql_name, "w") as outfile:
    print("""LOAD CSV WITH HEADERS FROM '%s'
             AS row
             CREATE (%s {%s})
             RETURN COUNT(row)""" % (csv_url, node_type, settings),
          file=outfile)
  print("Wrote %s" % cql_name, file=sys.stderr)

def prepare_update_cql(node_type, pk_col, dest):
  csv_name = os.path.join(dest, "update.csv")
  cql_name = os.path.join(dest, "update.cql")

  # match previous_pk
  # update, including changing pk

  props = get_header(csv_name)
  del props[windex(props, "previous_pk")]
  settings = ",\n  ".join(["node.%s = row.%s" % zz for zz in zip(props, props)])

  csv_url = "file:///%s" % (csv_name,)
  with open(cql_name, "w") as outfile:
    print("""LOAD CSV WITH HEADERS FROM '%s'
             AS row
             MATCH (node:%s) 
             WHERE node.%s = row.previous_pk
             SET %s
             RETURN COUNT(row)""" % (csv_url, node_type, pk_col, settings),
          file=outfile)
  print("Wrote %s" % cql_name, file=sys.stderr)
  pass

def prepare_remove_cql(node_type, pk_col, dest):
  # Now prepare the CQL file
  csv_name = os.path.join(dest, "remove.csv")
  cql_name = os.path.join(dest, "remove.cql")
  assert pk_col in get_header(csv_name)
  settings = "%s: row.previous_pk" % (pk_col,)
  csv_url = "file:///%s" % (csv_name,)
  with open(cql_name, "w") as outfile:
    print("""LOAD CSV WITH HEADERS FROM '%s'
             AS row
             DETACH DELETE (%s {%s})
             RETURN COUNT(row)""" % (csv_url, node_type, settings),
          file=outfile)
  print("Wrote %s" % cql_name, file=sys.stderr)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    TBD
    """)
  parser.add_argument('--delta',
                      help='name of directory holding add.csv and the other csv files')
  parser.add_argument('--pk', default="taxonID",
                      help='name of column containing primary key')
  args=parser.parse_args()
  write_cql_files(args.delta, args.pk)
