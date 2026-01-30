#!/bin/bash

set -v
set -e

# Generate the alignment for the COL/MDD comparison in the mss.

# First, perform the one-time manual setup ("installation").  Then run
# this file as an ordinary shell script using 'bash' (other shells
# might work too).  Output goes to ./rfb26.

# Personally, I compute the demo alignment with:
#   (W=run Q=~/g ~/g/listtols/doc/make-colmdd.sh)
# which puts the final alignment .csv file and temporary files in
# ./run.

# There are many other ways to do this.  Clones and CSVs can go in
# other places, customize the pipeline, different inputs, etc.

# Shell variables.
# Short names just to make this file easier on the eyes.
(($Q)) || L=~/g
(($L)) || L=$Q/listtools
(($M)) || M=$Q/MDD-DwC-mapping
(($A)) || A=$Q/artifacts
(($W)) || W=./run
mkdir -p $W

# Source files for the various list tools:
P=$L/src

# -----------------------------------------------------------------------------

# Create a local place to work and cd into it, if desired.  You can
# use a different location more compatible with your local setup.

# I've been setting up (once per computer) here:
mkdir -p $Q

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
unzip -u $A/col24-mammals/2b1541bc-12e7-4200-833b-7ae02e1d5f35.zip -d $W/col24-mammals
$P/clean.py --input `$P/find_taxa.py $W/col24-mammals` >$W/col24-mammals-clean.csv
$P/extract_names.py < $W/col24-mammals-clean.csv \
	| gnparser -s \
	| $P/use_gnparse.py --source $W/col24-mammals-clean.csv \
        > $W/col24-mammals.csv

# Prepare MDD v2.0 for alignment
$M/src/explore_data.py \
   --input $A/mdd2.0/MDD_v2.0_6759species.csv \
   --output $W/mdd2.0-dwc.csv
$P/clean.py --input $W/mdd2.0-dwc.csv > $W/mdd2.0-clean.csv
$P/extract_names.py < $W/mdd2.0-clean.csv \
	| gnparser -s \
	| $P/use_gnparse.py --source $W/mdd2.0-clean.csv > $W/mdd2.0.csv

time $P/plugin.py --A $W/col24-mammals.csv --B $W/mdd2.0.csv >report.csv
