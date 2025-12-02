# Sample run of the alignment tool.

# Perform manual steps.  Then run this file as an ordinary shell
# script using bash or any other Unix shell.

# MANUAL STEP: Get listtools.
# Go to https://github.com/jar398/listtools.
# Clone the repo (under 'Code' / 'Clone' / HTTP or the protocol of your choice).
# Edit L below to point to your local clone of listtools.

L=~/g/listtools
P=$L/src

# MANUAL STEP: Get COL 2024.
# Go to checklistbank.org.  Log in.  Go to 'Downloads'.  Find COL24.
# Download with these parameters: DwCA output, Mammalia root taxon.
# Save .zip file to $L/in/.

# Prepare CoL for alignment
unzip -u $L/in/2b1541bc-12e7-4200-833b-7ae02e1d5f35.zip -d col24-mammals
$P/clean.py --input `$P/find_taxa.py col24-mammals` >col24-mammals-clean.csv
$P/extract_names.py < col24-mammals-clean.csv \
	| gnparser -s \
	| $P/use_gnparse.py --source col24-mammals-clean.csv > col24-mammals.csv

# MANUAL STEP: Get tool for converting MDD files to DwC.
# Go to https://github.com/jar398/MDD-DwC-mapping.
# Clone the repo (under 'Code' / 'Clone' / HTTP or the protocol of your choice).
# Edit M below to point to your local clone or MDD-DwC-mapping.

M=~/g/MDD-DwC-mapping

# MANUAL STEP: Get MDD 2.0.
# Go to zenodo.org.  Search for "Mammal Diversity Database".  Select "Mammal
# Diversity Database".  Select "Version v2.0".  Under "Files" select
# "MDD_v2.0_6759species.csv".  Select "Download".  
# Find MDD_v2.0_6759species.csv on your computer (e.g. in 'Downloads'). 
# Move it to $L/in/.
# Run conversion program MDD-DwC-mapping/ on it.

# Prepare MDD 2.0 for alignment
$M/src/explore_data.py --input $L/in/MDD_v2.0_6759species.csv \
  --output mdd2.0-raw.csv
$P/clean.py --input mdd2.0-raw.csv >mdd2.0-clean.csv
$P/extract_names.py < mdd2.0-clean.csv \
	| gnparser -s \
	| $P/use_gnparse.py --source mdd2.0-clean.csv > mdd2.0.csv

# Align
$P/exemplar.py --A col24-mammals.csv --B mdd2.0.csv >col24-mammals-mdd2.0-exemplars.csv
$P/plugin.py --A col24-mammals.csv --B mdd2.0.csv --exemplars col24-mammals-mdd2.0-exemplars.csv \
  > alignment-report.csv
