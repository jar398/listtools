#!/bin/bash

# Generate the alignment for the COL/MDD comparison in a planned
# article by Rees, Franz, and Sterner, 2026.

# Here are steps you need to do in preparation to running this script.
# Review these sections before proceeding.
#
#  0. Get dependencies (regexp, gnparser) (see INSTALL.md)
#  1. Review the parameters, below, and override if necessary (especially Q)
#  2. Optional: to download checklists from their original sources, see 
#     instructions below
#
# The script will do its one-time setup if necessary, and then each
# time write temporary and output files to ./rfs or $W.

# File locations are subject to change.  In particular, if you have
# trouble obtaining the checklist sources, you might be able to find
# copies in the 'artifacts' repository in the 'jar389' account on
# github.

# There are many other ways to do this.  You can customize the
# pipeline, use different inputs, use 'make', etc.

set -e

# ---------------------------------------------------------------------------
# 1. PARAMETERS

# The setting of Q is my personal choice; you will likely want to
# override Q, and perhaps other parameter settings, either by changing
# this file or by exporting it from your shell before running this
# script.
# (I use short names to make this file easier on the eyes.)
#
# The shell syntax 'Q=path foo.sh' which sets Q while running foo.sh
# can be useful.

# From stackexchange:    VAR1="${VAR1:-default value}"

# $W (for 'work') is where you want to put the temporary files and
# final reports.
W="${W:-./rfs26}"

# $Q is a directory to contain resources you load from the Internet.
# On my setup this directory is the parent of a listtools repository
# clone, which you probably already have.  To avoid nesting
# $L -> $Q -> $L you might want to put Q outside of L, e.g. as
# $L's parent directory.
Q="${Q:-$W/g}"

# Git URL prefix.  Repositories are cloned as needed.
G="${G:-git@github.com:jar398}"

# L is the clone of the repository where you put or want to put
# listtools.
L="${L:-$Q/listtools}"

# M is where you have put or want to put the MDD to DW mapper.
M="${M:-$Q/MDD-DwC-mapping}"

# $A is a cache of other artifacts from the internet, either one you
# create manually or a clone of the jar398 artifacts repo.
A="${A:-$Q/artifacts}"

# Location of 'gnparser' utility
gnparser=gnparser

# -----------------------------------------------------------------------------
# 2. TO DOWNLOAD SOURCE CHECKLISTS MANUALLY:
#
#   If you choose to skip these steps, because they don't work or you
#   don't have the patience, the sources will be obtained from an
#   'artifacts' cache on github.
#
#   1. Choose a download directory for the sources, here written <path>.
#    mkdir -p <path>/col24-mammals (you choose <path>, must match parameter $A below)
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

# A place to work
mkdir -p $W
cd $W

# One-time setup:

# Configure parameters below (variables Q, L, M, A, W).

# Create a local place for the git clones, if necessary.
mkdir -p $Q

# Clone repos.  Substitute for these URLs for your own
# github authentication method.
test -e $L || \
   (cd $Q && git clone $G/listtools.git)

test -e $M || \
   (cd $Q && git clone $G/MDD-DwC-mapping.git)

test -e $A || \
   (cd $Q && git clone $G/artifacts.git)

# -----------------------------------------------------------------------------

# Finally, run the tools:

# Source files for the various list tools:
T=$L/src

# Prepare CoL for alignment.  You may have to adjust the URL.
unzip -u $A/col24-mammals/2b1541bc-12e7-4200-833b-7ae02e1d5f35.zip -d col24-mammals
$T/clean.py --input `$T/find_taxa.py col24-mammals` >col24-mammals-clean.csv
$T/extract_names.py < col24-mammals-clean.csv \
	| $gnparser -s \
	| $T/use_gnparse.py --source col24-mammals-clean.csv \
        > col24-mammals.csv

# Prepare MDD v2.0 for alignment
$M/src/explore_data.py \
   --input $A/mdd2.0/MDD_v2.0_6759species.csv \
   --output mdd2.0-dwc.csv
$T/clean.py --input mdd2.0-dwc.csv > mdd2.0-clean.csv
$T/extract_names.py < mdd2.0-clean.csv \
	| $gnparser -s \
	| $T/use_gnparse.py --source mdd2.0-clean.csv > mdd2.0.csv

time ($T/exemplar.py --A col24-mammals.csv --B mdd2.0.csv > exemplars.csv) \
  2>&1 | tee exemplars.dribble 

time ($T/align.py --A col24-mammals.csv --B mdd2.0.csv \
                  --exemplars exemplars.csv > report.csv) \
  2>&1 | tee report.dribble 

if false; then
  mkdir -p ~/g/artifacts/rfs26
  cp -p run/report.* ~/g/artifacts/rfs26/
fi
