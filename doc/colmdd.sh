#!/bin/bash

# Generate the alignment for the COL/MDD comparison in a planned
# article by Rees, Franz, and Sterner, 2026.  Inputs are the mammal
# subtree from Catalogue of Life 2024 and version 2.0 of the Mammal
# Diversity Database.  Outputs are a set of 'exemplars' used for
# occurrence codetermination and a species-level alignment report.
# By default the outputs are deposited in directory `./rfs26`.

# Here are steps you need to do in preparation to running this script.
# Review these sections before proceeding.
#
#  0. Get dependencies (regexp, gnparser) (see INSTALL.md)
#  1. Review the parameters, below, and override if necessary (especially `cache`)
#  2. Optional: to download checklists from their original sources, see 
#     instructions below
#
# The script will do its one-time setup if necessary, and then each
# time write (or overwrite) temporary and output files.

# Internet file locations are subject to change.  If you have trouble
# obtaining the original sources, you might be able to find copies in
# the 'artifacts' repository in the 'jar389' account on github.

# There are many other ways one might do this.  You could use
# different inputs, customize the pipeline, write your own 'makefile',
# etc.

# As of 2/2026 precomputed outputs are found here:
#   https://github.com/jar398/artifacts/tree/main/rfs26

set -e

# ---------------------------------------------------------------------------
# 1. PARAMETERS

# The defauls are my personal choices.  Please review and override for
# your local setup and preferences.

# The shell syntax 'var=path foo.sh' sets `var` while running foo.sh.
# If you override using this syntax or exported shell variables, you 
# can avoid having to modify this file.

# $work is where you want to put the temporary files and
# final alignment reports.
work="${work:-./rfs26}"

# $cache is a directory to contain resources you load from the
# Internet, mainly git clones.  To avoid circular nesting $tools ->
# $cache -> $tools you should put `cache` outside of `tools`, e.g. as
# its parent directory.
cache="${cache:-$work/g}"

# Git URL prefix.  Repositories are cloned as needed.
origin="${origin:-git@github.com:jar398}"

# tools is to be the clone of the listtools repository.
tools="${tools:-$cache/listtools}"

# mdd is where you have put or want to put the MDD to Darwin Core
# mapper.
mdd="${mdd:-$cache/MDD-DwC-mapping}"

# $artifacts is a cache of other artifacts from the internet, either one you
# create manually or a clone of the jar398 artifacts repo.
artifacts="${artifacts:-$cache/artifacts}"

# Location of 'gnparser' utility
gnparser=gnparser

# -----------------------------------------------------------------------------
# 2. TO DOWNLOAD SOURCE CHECKLISTS MANUALLY FROM ORIGINAL SOURCES:
#
#   If you choose to skip these steps, because they don't work or you
#   don't have the patience, the sources will be obtained from an
#   'artifacts' cache on github.
#
#   1. Choose a download directory for the sources, here written <path>.
#    mkdir -p <path>/col24-mammals (you choose <path>, must match parameter $artifacts below)
#    mkdir -p <path>/mdd2.0
#
#   2. Get COL 2024 Mammals:
#    Go to checklistbank.org.  Log in with your gbif.org user id.  Go to
#    'Downloads'.  Find COL24 (? the UI changes sometimes).  Prepare
#    with these parameters: DwCA output, Mammalia root taxon.  Download
#    the .zip file to <path>/col24-mammals/.
#
#   3. Get MDD v2.0:
#    Go to zenodo.org.  Search for "Mammal Diversity Database".  Select "Mammal
#    Diversity Database".  Select "Version v2.0".  Under "Files" select
#    "MDD_v2.0_6759species.csv".  Select "Download".  
#    Find MDD_v2.0_6759species.csv where your browser placed it on your 
#    computer (maybe in 'Downloads').  Move this csv file to <path>/mdd2.0/.

# -----------------------------------------------------------------------------
# SCRIPT

[ which gnparser ] || echo "No gnparser, see ../INSTALL.md"

# A place to work
mkdir -p $work
cd $work

# One-time setup:

# Configure parameters above (variables cache, origin, tools, mdd, artifacts, work).

# Create a local place for the git clones, if necessary.
mkdir -p $cache

# Clone repos.  Substitute for these URLs for your own
# github authentication method.
test -e $tools || \
   (cd $cache && git clone $origin/listtools.git)

test -e $mdd || \
   (cd $cache && git clone $origin/MDD-DwC-mapping.git)

# If you create artifacts "manually" they won't be retrieved from github
test -e $artifacts || \
   (cd $cache && git clone $origin/artifacts.git)

# -----------------------------------------------------------------------------

# Finally, run the tools:

# Source files for the various list tools:
run=$tools/src

# Prepare CoL for alignment.  You may have to adjust the URL.
unzip -u $artifacts/col24-mammals/2b1541bc-12e7-4200-833b-7ae02e1d5f35.zip -d col24-mammals
$run/clean.py --input `$run/find_taxa.py col24-mammals` >col24-mammals-clean.csv
$run/extract_names.py < col24-mammals-clean.csv \
	| $gnparser -s \
	| $run/use_gnparse.py --source col24-mammals-clean.csv \
        > col24-mammals.csv

# Prepare MDD v2.0 for alignment
$M/src/explore_data.py \
   --input $artifacts/mdd2.0/MDD_v2.0_6759species.csv \
   --output mdd2.0-dwc.csv
$run/clean.py --input mdd2.0-dwc.csv > mdd2.0-clean.csv
$run/extract_names.py < mdd2.0-clean.csv \
	| $gnparser -s \
	| $run/use_gnparse.py --source mdd2.0-clean.csv > mdd2.0.csv

time ($run/exemplar.py --A col24-mammals.csv --B mdd2.0.csv > exemplars.csv) \
  2>&1 | tee exemplars.dribble 

time ($run/align.py --A col24-mammals.csv --B mdd2.0.csv \
                  --exemplars exemplars.csv > report.csv) \
  2>&1 | tee report.dribble 

if false; then
  mkdir -p ~/g/artifacts/rfs26
  cp -p run/report.* ~/g/artifacts/rfs26/
fi

echo
echo "Script ran to ompletion."
echo "Exemplar definitions are in" $run/exemplars.csv
echo "Species alignment report is in" $run/report.csv
