# This is an example of the use of the list tools.

# Run with: 
#   make -f doc/diffpatch.makefile A=oldtable B=newtable
# The tables here would be oldtable.csv and newtable.csv

# Default example:
# Compare two versions of NCBI mammals

A=work/ncbi201505-mammals
B=work/ncbi202008-mammals
MAMMALIA=40674

# make A=work/ncbi201505-mammals B=work/ncbi202008-mammals

# time make A=work/ncbi201505 B=work/ncbi202008


# EOL examples:

# EOL DH 1.1 / 1.2 mammals only 
# A=work/dh11-mammals B=work/dh12-mammals MAMMALIA=EOL-000000627548

# EOL 0.9 / 1.1
# time make A=work/dh09 B=work/dh11

# EOL 0.9 / 1.1 mammals only
# time make A=work/dh09-mammals B=work/dh11-mammals MAMMALIA=EOL-000000627548

# EOL 1.1 / 1.2
# time make A=work/dh11 B=work/dh12

# Hierarchies
# DELTA_KEY=EOLid MANAGE=EOLid,parentEOLid,taxonID,landmark_status \
#   time make A=work/dh09-hier B=work/dh11-hier

# DELTA_KEY=EOLid MANAGE=EOLid,parentEOLid,taxonID,landmark_status \
#   time make A=work/dh11-hier B=work/dh12-hier

SHELL = /usr/bin/bash
RAKE = cd ../plotter && rake
P = src


# $< = first prerequisite

DELTA=work/delta-$(shell basename $A)-$(shell basename $B).csv
ROUND=work/round-$(shell basename $A)-$(shell basename $B).csv

all: $(ROUND)

# Columns that are managed by diff/patch
USAGE_KEY = taxonID
HIER_KEY = EOLid
DELTA_KEY ?= $(USAGE_KEY)

# Formerly: $P/project.py --keep $(MANAGE) <$< | ...
# and	    $P/sortcsv.py --key $(USAGE_KEY) <$< >$@.new

INDEX ?= taxonID,EOLid,scientificName,canonicalName
MANAGE ?= taxonID,scientificName,canonicalName,taxonRank,taxonomicStatus,nomenclaturalStatus,datasetID,source

$(DELTA): $A.csv $B.csv $P/delta.py
	@echo
	@echo "--- COMPUTING DELTA ---"
	set -o pipefail; \
	$P/delta.py --target $B.csv --pk $(DELTA_KEY) \
		    --index $(INDEX) --manage $(MANAGE) \
		    < $A.csv \
	| $P/sortcsv.py --key $(DELTA_KEY) \
	> $@.new
	@mv -f $@.new $@
	wc $@

# something:
# 	$P/scatter.py --dest $(basename $(DELTA)) < $(DELTA)

$(ROUND): $(DELTA) $A-narrow.csv $B-narrow.csv
	@echo
	@echo "--- APPLYING DELTA ---"
	set -o pipefail; \
	$P/apply.py --delta $< --pk $(DELTA_KEY) \
	    < $A-narrow.csv \
	| $P/sortcsv.py --key $(DELTA_KEY) \
	> $@.new
	@mv -f $@.new $@
	@echo "--- Comparing $@ to $B.csv ---"
	@wc $B-narrow.csv; wc $@

%-narrow.csv: %.csv
	$P/project.py --keep $(MANAGE) <$< >$@.new
	@mv -f $@.new $@

%-mammals.csv: %.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(MAMMALIA) < $< > $@.new
	@mv -f $@.new $@

# ----------------------------------------------------------------------
# EOL examples

inputs: dh work/dh09.csv work/dh11.csv
dh: work/dh09.csv work/dh11.csv work/dh12.csv
ASSEMBLY=prod

work/dh09.id:
	@mkdir -p work
	echo 1 >$@
work/dh11.id:
	@mkdir -p work
	echo 724 > $@

# about half a minute for DH 1.1

work/%.csv: work/%.id $P/start.py
	@mkdir -p work
	ID=$$(cat $<); \
	$P/start.py --input `$(RAKE) resource:taxa_path \
	       	          CONF=$(ASSEMBLY) \
		          ID=$$ID` \
		    --pk $(USAGE_KEY) \
	| $P/sortcsv.py --key $(USAGE_KEY) > $@.new
	@mv -f $@.new $@

DH12_LP="https://opendata.eol.org/dataset/tram-807-808-809-810-dh-v1-1/resource/02037fde-cc69-4f03-94b5-65591c6e7b3b"

work/dh12.csv: $P/start.py
	@mkdir -p work
	$P/start.py --input `$(RAKE) dwca:taxa_path OPENDATA=$(DH12_LP)` \
		    --pk taxonID \
	       > $@.new
	@mv -f $@.new $@

# in1=./deprecated/work/1-mam.csv
# in2=./deprecated/work/724-mam.csv

# EOL mammals root = page id 1642, usage id EOL-000000627548 (DH 1.1)
# MAMMALIA=EOL-000000627548

# Mammals root has a different usage id in DH 0.9
MAMMALIA09=-168130
work/dh09-mammals.csv: work/dh09.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(MAMMALIA09) < $< > $@.new
	@mv -f $@.new $@

work/dh11-mammals-hier.csv: work/dh11-mammals.csv work/dh11-map.csv $P/hierarchy.py
	$P/hierarchy.py --mapping work/dh11-map.csv \
		  < $< \
		  > $@.new
	@mv -f $@.new $@

# EOL dynamic hierarchy - usages mapped to pages

work/%-hier.csv: work/%.csv work/%-map.csv $P/hierarchy.py
	set -o pipefail; \
	$P/hierarchy.py --mapping $(<:.csv=-map.csv) \
			--keep landmark_status \
		  < $< \
	| $P/sortcsv.py --key $(HIER_KEY) > $@.new
	@mv -f $@.new $@

work/%-map.csv: work/%.id
	ID=$$(cat $<); \
	cp `$(RAKE) resource:map CONF=$(ASSEMBLY) ID=$$ID` $@.new
	@mv -f $@.new $@

work/dh12-map.csv: work/dh11-map.csv
	cp $< $@

# Deprecated ... ?

work/%-mapped.csv: work/%.csv work/%-map.csv $P/idmap.py
	$P/idmap.py --mapping $(<:.csv=-map.csv) \
		  < $< > $@.new
	@mv -f $@.new $@

# ----------------------------------------------------------------------
# ASU example

# N.b. NCBI taxonomy id for mammals is MAMMALIA=40674

foo: work/ncbi201505.csv

work/ncbi201505.url:
	echo ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2015-05-01.zip \
	  >$@
work/ncbi202008.url:
	echo ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2020-08-01.zip \
	  >$@

%/dump/names.dmp: %.url
	mkdir -p `dirname $@`
	wget -O `dirname $@`.zip $$(cat $<)
	unzip -d `dirname $@` `dirname $@`.zip
	touch `dirname $@`/*
.PRECIOUS: %/dump/names.dmp

# Convert NCBI taxdump to DwC form
%.csv: %/dump/names.dmp src/ncbi_to_dwc.py 
	$P/ncbi_to_dwc.py `dirname $<` \
	| $P/start.py --pk taxonID \
	| $P/sortcsv.py --key $(USAGE_KEY) > $@.new
	@mv -f $@.new $@
