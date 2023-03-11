#!/usr/bin/env python3

# Sources and sinks of CSV rows

import sys, csv, newick

def open(specifier, mode='r'):
  if not specifier: specifier = '-'
  return Rows(specifier, mode)

class Rows:
  def __init__(self, x, mode): 
    self.x = x

    if x.startswith('('):
      self.file = None
      assert not 'w' in mode
    else:
      if x == '-':
        self.close_required = False
        self.file = sys.stdout if 'w' in mode else sys.stdin
      else:
        self.close_required = True
        self.file = open(x, mode)
    self.open = True

  def __enter__(self):
    if self.file and self.close_required:
      self.file.__enter__()
    return self
      
  def __exit__(self, exc_type, exc_val, exc_tb):
    if self.file and self.close_required:
      self.file.__exit__(exc_type, exc_val, exc_tb)
    self.open = False

  def rows(self):             # returns a row generator
    if self.file:
      yield from csv.reader(self.file)
    else:
      yield from newick.parse_newick(self.x)

  def write_rows(self, row_generator):
    assert self.file
    writer = csv.writer(self.file)
    for row in row_generator:
      assert self.open
      writer.writerow(row)

# Test

if __name__ == '__main__':
  with open('(a,b*)c') as r:
    with open('-', 'w') as w:
      w.write_rows(r.rows())
