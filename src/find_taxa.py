#!/usr/bin/env python3

# Get the name of the DwC taxa file in a GBIF or CoL DwCA archive.
# GBIF moved the taxa file from dump/ to dump/backbone/ in 2022.
# Argument is path to directory where dumped.

# Ideally we should find meta.xml and read it to find out the name,
# but that would require me learning how to parse and crawl XML in
# python.

import sys, os

def find(dir):
  for dir in (dir, os.path.join(dir, "backbone")):
    for f in ("taxon", "taxa", "Taxon", "Taxa"):
      for x in ("csv", "tsv", "txt"):
        name = os.path.join(dir, "%s.%s" % (f, x))
        if os.path.isfile(name):
          return name

dir = sys.argv[1]
name = find(dir)
if name:
  print(name)
else:
  print("Didn't find a .taxon file or similar in %s" % dir, file=sys.stderr)
  sys.exit(1)

