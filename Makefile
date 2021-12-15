# This is an example of the use of the list tools.

# Run with: 
#   make -f doc/diffpatch.makefile A=oldtable B=newtable
# The tables here would be oldtable.csv and newtable.csv

all:
	@echo "Please specify a target:"
	@echo "  diff      show new/removed/changed records"
	@echo "  round     demo diff/patch round trip"
	@echo "  report    report on taxon concept differences"
	@echo "Use make parameter syntax to specify the two checklists, e.g."
	@echo "  make A=work/ncbi201505-mammals B=work/ncbi202008-mammals report"
	@echo "(when the .csv file for work/foo is in work/foo.csv etc)"

# Default example:
# Compare two versions of NCBI mammals

A?=work/ncbi201505-mammals
B?=work/ncbi202008-mammals

# make A=work/ncbi201505-mammals B=work/ncbi202008-mammals demo
# time make A=work/ncbi201505 B=work/ncbi202008 demo
#  etc. etc.

# N.b. record_id has a specific meaning in this system.  The id in
#  some managed namespace, prefixed by something to identify the namespace.
# Example syntax:
#   NCBI:123, EOL:4567 (for page id), GBIF:8910, WORMS:2345   etc.

# ----- NCBI and GBIF examples:

# make A=work/ncbi201505-mammals B=work/ncbi202008-mammals report
# make A=work/ncbi201505 B=work/ncbi202008 report  # will take forever

# make A=work/gbif20190916-mammals B=work/gbif20210303-mammals report
# make A=work/ncbi202008-mammals B=work/gbif20210303-mammals report
# and so on.

# ----- BioKIC/ATCR examples:

# make A=work/mdd1.2-mammals B=work/mdd1.3 report
# make A=work/mdd1.2 B=work/mdd1.3 report
# make A=work/mdd1.6 B=work/mdd1.7 report
# make A=work/gbif20210303-mammals B=work/mdd1.0-mammals report

# ----- EOL examples:

# time make A=work/dh11-mammals B=work/dh12-mammals round
# time make A=work/dh09-mammals B=work/dh11-mammals round
# time make A=work/dh09 B=work/dh11 round
# time make A=work/dh11 B=work/dh12 round

# Hierarchies
# DELTA_KEY=EOLid MANAGE=EOLid,parentEOLid,taxonID,landmark_status \
#   time make A=work/dh09-hier B=work/dh11-hier report

# DELTA_KEY=EOLid MANAGE=EOLid,parentEOLid,taxonID,landmark_status \
#   time make A=work/dh11-hier B=work/dh12-hier report

# ----- Parameters

# $< = first prerequisite, $@ = target

SHELL?=/usr/bin/bash
RAKE?=cd ../plotter && rake
P?=src

INDEX?="scientificName,tipe,record_id,EOLid,canonicalName,canonicalStem"

# Primary key column
PRIMARY_KEY?=taxonID

# overridden to EOLid for EOL
DELTA_KEY?=$(PRIMARY_KEY)

# What is MANAGE_KEY for? I forget

# ----- General rules

DELTA=work/delta-$(shell basename $A)-$(shell basename $B).csv
ROUND=work/round-$(shell basename $A)-$(shell basename $B).csv

round: $(ROUND)
delta: $(DELTA)

$(RM): $A-gnp.csv $B-gnp.csv $P/match_records.py
	@echo
	@echo "--- COMPUTING RECORD MATCHES ---"
	$P/match_records.py --target $B-gnp.csv --pk $(DELTA_KEY) --index $(INDEX) \
		    < $A-gnp.csv > $@.new
	@mv -f $@.new $@

# Formerly: $P/project.py --keep $(MANAGE) <$< | ...
# and	    $P/sortcsv.py --key $(PRIMARY_KEY) <$< >$@.new

MANAGE?=taxonID,scientificName,canonicalName,taxonRank,taxonomicStatus,nomenclaturalStatus,datasetID,source

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

# Make the file smaller by eliminating unneeded columns

%-narrow.csv: %.csv $P/project.py
	$P/project.py --keep $(MANAGE) <$< | \
	$P/sortcsv.py --key $(PRIMARY_KEY) >$@.new
	@mv -f $@.new $@

%-mammals-gnp.csv: %-mammals.csv

# ----- NCBI-specific rules

MAMMALIA_NCBI=40674

ncbi%-mammals.csv: ncbi%.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(MAMMALIA_NCBI) < $< > $@.new
	@mv -f $@.new $@

gbif%-mammals.csv: gbif%.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(MAMMALIA_GBIF) < $< > $@.new
	@mv -f $@.new $@

dh1%-mammals.csv: dh1%.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(MAMMALIA_DH11) < $< > $@.new
	@mv -f $@.new $@

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
%.csv: %/dump/names.dmp src/ncbi_to_dwc.py $P/start.py
	$P/ncbi_to_dwc.py `dirname $<` \
	| $P/start.py --pk $(PRIMARY_KEY) \
	| $P/sortcsv.py --key $(PRIMARY_KEY) > $@.new
	@mv -f $@.new $@

# Adjoin tipe column.  To skip this change the rule to
# simply copy %.csv to %-gnp.csv above. 
%-gnp.csv: %.csv $P/extract_names.py $P/use_parse.py
	$P/extract_names.py < $< \
	| gnparser -s \
	| $P/use_parse.py --source $< > $@.new
	@mv -f $@.new $@

# ----- GBIF-specific rules

MAMMALIA_GBIF=359

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

# Ingest GBIF dump, convert TSV to CSV
%.csv: %/dump/meta.xml $P/start.py
	$P/start.py --pk $(PRIMARY_KEY) --input `dirname $<`/Taxon.tsv \
	  --managed taxonID --prefix GBIF: >$@.new
	@mv -f $@.new $@
.PRECIOUS: %.csv

work/gbif20210303-mammals-gnp.csv: work/gbif20210303-mammals.csv

# ----------------------------------------------------------------------
# ASU/BioKIC example

# MDD

work/mdd1.7.csv: $(HOME)/Downloads/MDD_DwC_versions/MDD_v1.7_6567species_inDwC.csv
	cp -p $< $@
work/mdd1.7-gnp.csv: work/mdd1.7.csv
work/mdd1.6.csv: $(HOME)/Downloads/MDD_DwC_versions/MDD_v1.6_6557species_inDwC.csv
	cp -p $< $@
work/mdd1.6-gnp.csv: work/mdd1.6.csv
work/mdd1.5.csv: $(HOME)/Downloads/MDD_DwC_versions/MDD_v1.5_6554species_inDwC.csv \
		 $P/start.py
	$P/start.py < $< --pk taxonID > $@.new
	@mv -f $@.new $@
work/mdd1.5-gnp.csv: work/mdd1.5.csv
work/mdd1.4.csv: $(HOME)/Downloads/MDD_DwC_versions/MDD_v1.4_6533species_inDwC.csv \
		 $P/start.py
	$P/start.py < $< --pk taxonID > $@.new
	@mv -f $@.new $@
work/mdd1.4-gnp.csv: work/mdd1.4.csv
work/mdd1.3.csv: $(HOME)/Downloads/MDD_DwC_versions/MDD_v1.3_6513species_inDwC.csv
	cp -p $< $@
work/mdd1.3-gnp.csv: work/mdd1.3.csv
work/mdd1.31.csv: $(HOME)/Downloads/MDD_DwC_versions/MDD_v1.31_6513species_inDwC.csv
	cp -p $< $@
work/mdd1.31-gnp.csv: work/mdd1.3.csv
work/mdd1.2.csv: $(HOME)/Downloads/MDD_DwC_versions/MDD_v1.2_6485species_inDwC.csv
	cp -p $< $@
work/mdd1.2-gnp.csv: work/mdd1.2.csv
work/mdd1.1.csv: $(HOME)/Downloads/MDD_DwC_versions/MDD_v1.1_6526species_inDwC.csv
	cp -p $< $@
work/mdd1.1-gnp.csv: work/mdd1.1.csv
work/mdd1.0.csv: $(HOME)/Downloads/MDD_DwC_versions/MDD_v1_6495species_JMamm_inDwC.csv \
		 $P/start.py
	$P/start.py < $< --pk taxonID > $@.new
	@mv -f $@.new $@
work/mdd1.0-gnp.csv: work/mdd1.0.csv

# ----- Record match rules, mainly for EOL purposes

MAMMALIA_DH11=EOL-000000627548
MAMMALIA_DH09=-168130

RM=work/rm-$(shell basename $A)-$(shell basename $B).csv
ALIGNMENT=work/alignment-$(shell basename $A)-$(shell basename $B).csv
REPORT=work/report-$(shell basename $A)-$(shell basename $B).csv

report: $(REPORT)
alignment: $(ALIGNMENT)

$(ALIGNMENT): $(RM) $P/align.py $P/property.py
	@echo
	@echo "--- COMPUTING ALIGNMENT ---"
	$P/align.py --target $B-gnp.csv --matches $(RM) \
		    < $A-gnp.csv > $@.new
	@mv -f $@.new $@

REPORT_OPTION ?=

$(REPORT): $(ALIGNMENT) $(RM) $P/report.py $P/property.py
	@echo
	@echo "--- PREPARING REPORT ---"
	$P/report.py --source $A-gnp.csv --alignment $(ALIGNMENT) \
		     --matches $(RM) \
		     $(REPORT_OPTIONS) \
		     < $B-gnp.csv > $@.new
	@mv -f $@.new $@

# ----------------------------------------------------------------------
# EOL examples

# EOL mammals root = page id 1642

HIER_KEY=EOLid

inputs: dh work/dh09.csv work/dh11.csv
dh: work/dh09.csv work/dh11.csv work/dh12.csv
ASSEMBLY=prod

work/dh09.eol-resource-id:
	@mkdir -p work
	echo 1 >$@
work/dh11.eol-resource-id:
	@mkdir -p work
	echo 724 > $@

# about half a minute for DH 1.1

work/%.csv: work/%.eol-resource-id $P/start.py
	@mkdir -p work
	ID=$$(cat $<); \
	$P/start.py --input `$(RAKE) resource:taxa_path \
	       	          CONF=$(ASSEMBLY) \
		          ID=$$ID` \
		    --pk $(PRIMARY_KEY) \
	| $P/sortcsv.py --key $(PRIMARY_KEY) > $@.new
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

# Mammals root has a different taxonID in DH 0.9
work/dh09-mammals.csv: work/dh09.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(MAMMALIA_DH09) < $< > $@.new
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

work/%-map.csv: work/%.eol-resource-id
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

