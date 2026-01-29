# Generate the alignment for the COL/MDD comparison in the mss.

# I've been setting up (once per computer) with:
#   mkdir -p ~/g
#   cd ~/g
#   clone ... URL for listtools repo ...
#   clone ... URL for MDD-DwC-mapping repo ...
#   (cd MDD-DwC-mapping; make)
#   cd listtools; mkdir -p run

# Then I compute the alignment with:
#   (pushd run; L=.. ../doc/colmdd/colmdd.sh)
# which puts the alignment CSV file in run/.

# There are many other ways to do this.  Clones and CSVs can go in
# other places, customize the pipeline, different inputs, etc.

# Directory containing listtools (github checkout)
L=~/g/listtools

# Source files for execution
P=$L/src

# COL - there is a copy sequestered in the in/ directory.

# How to make col-mammals.csv ?  Get mammals from COL24 from checklistbank,
# then extract the taxa file.
# (Before that: get mammals zip file from checklistbank.)
# I think what we need is in listtools/in/col24.7.tsv ?
$P/clean.py --input $L/in/col24.7.tsv >col24-mammals-clean.csv
$P/extract_names.py < col24-mammals-clean.csv \
	| gnparser -s \
	| $P/use_gnparse.py --source col24-mammals-clean.csv > col24-mammals.csv

# MDD

# Start with listtools/in/MDD_v2.0_6759species.csv
# (Before that: Go to Zenodo to get MDD 2.0.  get zip file.
#  Convert MDD format to DwC.)

$L/../MDD-DwC-mapping/src/explore_data.py \
   --input $L/in/MDD_v2.0_6759species.csv \
   --output mdd2.0-dwc.csv
$P/clean.py --input mdd2.0-dwc.csv >mdd2.0-clean.csv
$P/extract_names.py < mdd2.0-clean.csv \
	| gnparser -s \
	| $P/use_gnparse.py --source mdd2.0-clean.csv > mdd2.0.csv

time $P/plugin.py --A col24-mammals.csv --B mdd2.0.csv >report.csv
