#!/bin/bash

# Generate the alignment for the COL/MDD comparison in the mss.

# First, perform the one-time manual setup ("installation").  Then run
# this file as an ordinary shell script using 'bash' (other shells
# might work too).  Output goes to ./rfs26.

# Personally, I compute the demo alignment with:
#   (W=run Q=~/g ~/g/listtols/doc/colmdd.sh)
# which puts the final alignment .csv file and temporary files in
# ./rfs26.

# There are many other ways to do this.  Clones and CSVs can go in
# other places, customize the pipeline, different inputs, etc.

# Shell variables.
# Short names just to make this file easier on the eyes.
(($Q)) || Q=~/g
(($L)) || L=$Q/listtools
(($M)) || M=$Q/MDD-DwC-mapping
(($A)) || A=$Q/artifacts
(($W)) || W=./rfs26
# Source files for the various list tools:
P=$L/src

# A place to work
mkdir -p $W
cd $W

# -----------------------------------------------------------------------------

# One-time setup:

# Create a local place to work and cd into it, if desired.  You can
# use a different location more compatible with your local setup.

# Get dependencies (regexp, gnparser) (see documentation ??)

# MANUAL STEP: Clone repos.  Substitute for these URLs for your own
# github authentication method.
#   (cd $Q && git clone git@github.com:jar398/listtools.git)
#   (cd $Q && git clone git@github.com:jar398/MDD-DwC-mapping.git)

# MANUAL STEP: Get source checklists.
#  GET SOURCE CHECKLISTS METHOD A:
#   (cd $Q && git clone git@github.com:jar398/artifacts.git)
#
#  GET SOURCE CHECKLISTS METHOD B:
#   Get COL 2024 Mammals:
#   mkdir -p $A
#   Go to checklistbank.org.  Log in with GBIF user id.  Go to
#   'Downloads'.  Find COL24 (? the UI changes sometimes).  Prepare
#   with these parameters: DwCA output, Mammalia root taxon.  Download
#   the .zip file to $A/col24/.
#
#   Get MDD v2.0:
#   Go to zenodo.org.  Search for "Mammal Diversity Database".  Select "Mammal
#   Diversity Database".  Select "Version v2.0".  Under "Files" select
#   "MDD_v2.0_6759species.csv".  Select "Download".  
#   Find MDD_v2.0_6759species.csv on your computer (e.g. in 'Downloads'). 
#   Move it to $A/mdd2.0/.

# -----------------------------------------------------------------------------

# Prepare CoL for alignment.  You may have to adjust the URL.
unzip -u $A/col24-mammals/2b1541bc-12e7-4200-833b-7ae02e1d5f35.zip -d col24-mammals
$P/clean.py --input `$P/find_taxa.py col24-mammals` >col24-mammals-clean.csv
$P/extract_names.py < col24-mammals-clean.csv \
	| gnparser -s \
	| $P/use_gnparse.py --source col24-mammals-clean.csv \
        > col24-mammals.csv

# Prepare MDD v2.0 for alignment
$M/src/explore_data.py \
   --input $A/mdd2.0/MDD_v2.0_6759species.csv \
   --output mdd2.0-dwc.csv
$P/clean.py --input mdd2.0-dwc.csv > mdd2.0-clean.csv
$P/extract_names.py < mdd2.0-clean.csv \
	| gnparser -s \
	| $P/use_gnparse.py --source mdd2.0-clean.csv > mdd2.0.csv

time ($P/align.py --A col24-mammals.csv --B mdd2.0.csv > report.csv) \
  2>&1 | tee report.dribble 

if False; do
  mkdir -p ~/g/artifacts/rfs26
  cp -p run/report.* ~/g/artifacts/rfs26/
done
