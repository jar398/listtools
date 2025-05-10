# CoL 2024 mammals vs. MDD

W=work
A=$W/COL24.7
B=$W/MDD2.0
ABW=$W/$(notdir $A)_$(notdir $B)

include ../listtools.mk

all: $(ABW)_report.csv

# For each .zip file:
#   extract contents
#   get taxon file out of .zip file
# For each taxon .csv file:
#   clean up taxon file
#   extract_names, gnparser, use_gnparse

# COL (we use version 24.7 mammals) comes in a DwCA .zip file
# $(@F) = The file-within-directory part of the file name of the target. 
# $(<F) = The directory part and the file-within-directory part of the 
#   first prerequisite. 
A_zip?=in/50ed1b9d-9e7b-4aee-b611-0582ccb24d90.zip
$A.zip: $(A_zip)
	mkdir -p $W
	(cd $W && ln -sf ../$< $(@F))

# in work:
# link  ../in/50ed1b9d-9e7b-4aee-b611-0582ccb24d90.zip
# to	COL24.7.zip

$A.dump: $A.zip

$A.csv: $A.dump

# MDD does not have a zip file.  We need a DwC taxon file, and that
# has to be created using a conversion script.

$B-clean.csv: $B-raw.csv
	$P/clean.py --input $< >$@

$B-raw.csv: in/MDD_v2.0_6759species.csv
	../../MDD-DwC-mapping/src/explore_data.py \
	  --input $< \
	  --output $@
