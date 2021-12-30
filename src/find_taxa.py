#!/usr/bin/env python3

# Get the name of the DwC taxa file
# Argument is directory containing meta.xml

# Ideally we read meta.xml to find out the name

import sys, os

def find(dir):
  for dir in (dir, os.path.dirname(dir)):
    for f in ("taxon", "taxa", "Taxon", "Taxa"):
      for x in ("csv", "tsv", "txt"):
        name = os.path.join(dir, "%s.%s" % (f, x))
        if os.path.isfile(name):
          return name

name = find(sys.argv[1])
if name:
  print(name)
else:
  print("Nope", file=sys.stderr)
  sys.exit(1)

