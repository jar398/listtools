# This is an example of the use of the list tools.

# Run with: 
#   make -f doc/diffpatch.makefile A=oldtable B=newtable
# The tables here would be oldtable.csv and newtable.csv

# Default example:
# Compare two versions of NCBI mammals

A ?= work/ncbi201505-mammals
B ?= work/ncbi202008-mammals
MAMMALIA=40674

# make A=work/ncbi201505-mammals B=work/ncbi202008-mammals

# time make A=work/ncbi201505 B=work/ncbi202008


# EOL examples:

# EOL DH 1.1 / 1.2 mammals only 
# time make A=work/dh11-mammals B=work/dh12-mammals MAMMALIA=EOL-000000627548

# EOL 0.9 / 1.1 mammals only
# time make A=work/dh09-mammals B=work/dh11-mammals \
#   MAMMALIA=EOL-000000627548 INDEX=EOLid,scientificName,canonicalName

# EOL 0.9 / 1.1
# time make A=work/dh09 B=work/dh11

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

INDEX ?= taxonID,EOLid,scientificName,canonicalName,canonicalStem,altKey
MANAGE ?= taxonID,scientificName,canonicalName,taxonRank,taxonomicStatus,nomenclaturalStatus,datasetID,source

$(DELTA): $A.csv $B.csv $P/delta.py $P/match_records.py $P/property.py
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

$(ROUND): $(DELTA) $A-narrow.csv $B-narrow.csv $P/apply.py
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

%-narrow.csv: %.csv $P/project.py
	$P/project.py --keep $(MANAGE) <$< | \
	$P/sortcsv.py --key $(USAGE_KEY) >$@.new
	@mv -f $@.new $@

%-mammals.csv: %.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(MAMMALIA) < $< > $@.new
	@mv -f $@.new $@

%-mammals-gnp.csv: %-mammals.csv

RM=work/rm-$(shell basename $A)-$(shell basename $B).csv
ALIGNMENT=work/alignment-$(shell basename $A)-$(shell basename $B).csv
REPORT=work/report-$(shell basename $A)-$(shell basename $B).csv

report: $(REPORT)
alignment: $(ALIGNMENT)

$(RM): $A-gnp.csv $B-gnp.csv $P/match_records.py
	@echo
	@echo "--- COMPUTING RECORD MATCHES ---"
	$P/match_records.py --target $B-gnp.csv --pk $(DELTA_KEY) --index $(INDEX) \
		    < $A-gnp.csv > $(RM)

$(ALIGNMENT): $(RM) $P/align.py $P/property.py
	@echo
	@echo "--- COMPUTING ALIGNMENT ---"
	$P/align.py --target $B-gnp.csv --matches $(RM) \
		    < $A-gnp.csv > $(ALIGNMENT)

$(REPORT): $(ALIGNMENT) $P/report.py $P/property.py
	@echo
	@echo "--- PREPARING REPORT ---"
	$P/report.py --source $A-gnp.csv --alignment $(ALIGNMENT) \
		    < $B-gnp.csv > $(REPORT)

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
	$P/hierarchy.py --mapping $(basename $<)-map.csv) \
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
	$P/idmap.py --mapping $(basename $<)-map.csv) \
		  < $< > $@.new
	@mv -f $@.new $@

# ----------------------------------------------------------------------
# ASU example

# N.b. NCBI taxonomy id for mammals is MAMMALIA=40674

foo: work/ncbi201505.csv

work/ncbi201505.ncbi-url:
	echo ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2015-05-01.zip \
	  >$@
work/ncbi202008.ncbi-url:
	echo ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2020-08-01.zip \
	  >$@

%/dump/names.dmp: %.ncbi-url
	mkdir -p `dirname $@`
	wget -O `dirname $@`.zip $$(cat $<)
	unzip -d `dirname $@` `dirname $@`.zip
	touch `dirname $@`/*
.PRECIOUS: %/dump/names.dmp

# Convert NCBI taxdump to DwC form
%.csv: %/dump/names.dmp src/ncbi_to_dwc.py 
	$P/ncbi_to_dwc.py `dirname $<` \
	| $P/start.py --pk $(USAGE_KEY) \
	| $P/sortcsv.py --key $(USAGE_KEY) > $@.new
	@mv -f $@.new $@

# Adjoin altKey column.  To skip this change the rule to simply copy %.csv to %-gnp.csv above.
%-gnp.csv: %.csv $P/extract_names.py $P/use_parse.py
	$P/extract_names.py < $< \
	| gnparser -s \
	| $P/use_parse.py --source $< > $@.new
	@mv -f $@.new $@

# MDD

work/mdd1.7.csv: $(HOME)/Downloads/MDD_DwC_versions/MDD_v1.7_6567species_inDwC.csv
	cp -p $< $@
work/mdd1.7-gnp.csv: work/mdd1.7.csv
work/mdd1.6.csv: $(HOME)/Downloads/MDD_DwC_versions/MDD_v1.6_6557species_inDwC.csv
	cp -p $< $@
work/mdd1.6-gnp.csv: work/mdd1.6.csv
work/mdd1.2.csv: $(HOME)/Downloads/MDD_DwC_versions/MDD_v1.2_6485species_inDwC.csv
	cp -p $< $@
work/mdd1.2-gnp.csv: work/mdd1.2.csv
work/mdd1.1.csv: $(HOME)/Downloads/MDD_DwC_versions/MDD_v1.1_6526species_inDwC.csv
	cp -p $< $@
work/mdd1.1-gnp.csv: work/mdd1.1.csv


# GBIF, mammals = 359

work/gbif20190916.gbif-url:
	echo https://rs.gbif.org/datasets/backbone/2019-09-06/backbone.zip >$@

work/gbif20210303.gbif-url:
	echo https://rs.gbif.org/datasets/backbone/2021-03-03/backbone.zip >$@

%/dump/meta.xml: %.gbif-url
	mkdir -p `dirname $@`/dump
	wget -O `dirname $@`.zip $$(cat $<)
	unzip -d `dirname $@` `dirname $@`.zip
	touch `dirname $@`/*
.PRECIOUS: %/dump/meta.xml

# Ingest GBIF dump
%.csv: %/dump/meta.xml
	$P/start.py --pk $(USAGE_KEY) --input `dirname $<`/Taxon.tsv >$@.new
	@mv -f $@.new $@
.PRECIOUS: %.csv

fuu: work/gbif20210303-mammals.csv

fuuu:
	@echo $(dir a/b/c.d)
	@echo $(basename a/b/c.d)
