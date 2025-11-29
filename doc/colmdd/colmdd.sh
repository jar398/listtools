# Might require bash and/or gnu make

# Go to checklistbank.org.  Log in.  Go to 'Downloads'.  Find COL24.
# Download with these parameters: DwCA output, Mammalia root taxon.
# Save .zip file to $L/in/.

# Go to zenodo.org.  Search for "Mammal diversity".  Select "Mammal
# diversity database".  Select "Version v2.0".  Under "Files" select
# "MDD_v2.0_6759species.csv".  Click "Download".  
# Find MDD_v2.0_6759species.csv on your computer (e.g. in 'Downloads'). 
# Move it to $L/in/.
# Run conversion program MDD-DwC-mapping/ on it.

L=~/g/listtools
M=~/g/MDD-DwC-mapping
P=$L/src

# Prepare CoL for alignment
unzip $L/in/2b1541bc-12e7-4200-833b-7ae02e1d5f35.zip -d col24-mammals
$P/clean.py --input `$P/find_taxa.py col24-mammals` >col24-mammals-clean.csv
$P/extract_names.py < col24-mammals-clean.csv \
	| gnparser -s \
	| $P/use_gnparse.py --source col24-mammals-clean.csv > col24-mammals.csv

# Prepare MDD for alignment
$M/src/explore_data.py --input $L/in/MDD_v2.0_6759species.csv \
  --output mdd2.0-raw.csv
$P/clean.py --input mdd2.0-raw.csv >mdd2.0-clean.csv
$P/extract_names.py < mdd2.0-clean.csv \
	| gnparser -s \
	| $P/use_gnparse.py --source mdd2.0-clean.csv > mdd2.0.csv

# Align
$P/exemplar.py --A col24-mammals.csv --B mdd2.0.csv >exemplars.csv
$P/plugin.py --A col24-mammals.csv --B mdd2.0.csv --exemplars exemplars.csv \
  > alignment-report.csv
