# $< = first prerequisite, $@ = target

# This has examples of the use of the list tools.
# Very much in flux!

all:
	@echo "Please specify a target:"
	@echo "  diff       show new/removed/changed records"
	@echo "  round      demo diff/patch round trip"
	@echo "  report     report on NCBI extensions differences 2015-2020"
	@echo "  mdd-report report on MDD extensions differences 1.6-1.7"
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

# N.b. managed_id has a specific meaning in this system; it's an id in
#  some managed namespace, prefixed by something to identify the namespace.
# Example syntax:
#   ncbi:123, eol:4567 (for page id), gbif:8910, worms:2345   etc.

# ----- 1. NCBI example:

ncbi-report:
	$(MAKE) A=work/ncbi201505-mammals B=work/ncbi202008-mammals report

# make A=work/ncbi201505 B=work/ncbi202008 report  # big, will take forever

# ----- 2. GBIF examples:

gbif-report:
	$(MAKE) A=work/gbif20190916-mammals B=work/gbif20210303-mammals report

# make A=work/ncbi202008-mammals B=work/gbif20210303-mammals report
# and so on.

# ----- 3. BioKIC/ATCR examples:

mdd-report:
	$(MAKE) A=work/mdd1.6 B=work/mdd1.7 report

# make A=work/mdd1.2-mammals B=work/mdd1.3 report
# make A=work/mdd1.2 B=work/mdd1.3 report
# make A=work/gbif20210303-mammals B=work/mdd1.0-mammals report

# ----- 4. EOL examples:

# Requires clone of 'plotter' repo:

eol-report:
	$(MAKE) A=work/dh11-mammals B=work/dh12-mammals report

# time make A=work/dh11-mammals B=work/dh12-mammals round
# time make A=work/dh09-mammals B=work/dh11-mammals round
# time make A=work/dh11 B=work/dh12 round
# time make A=work/dh09 B=work/dh11 round

# Hierarchies - columns are those that neo4j needs to know
# DELTA_KEY=EOLid MANAGE=EOLid,parentEOLid,taxonID,landmark_status \
#   time make A=work/dh09-hier B=work/dh11-hier report

# DELTA_KEY=EOLid MANAGE=EOLid,parentEOLid,taxonID,landmark_status \
#   time make A=work/dh11-hier B=work/dh12-hier report

# ----- General parameters

SHELL?=/usr/bin/bash
P?=src

INDEX?="scientificName,type,canonicalName,canonicalStem,managed_id"

# Primary key column
PRIMARY_KEY?=taxonID

# overridden to EOLid for EOL
DELTA_KEY?=$(PRIMARY_KEY)

# projection, to make the files smaller
KEEP?="taxonID,canonicalName,scientificName,type,canonicalStem,managed_id,parentNameUsageID,acceptedNameUsageID"

# not sure exactly what this means.  Add EOLid for EOL
MANAGED?=$(KEEP)

# This is for EOL.  Requires clone of 'plotter' repo
RAKE?=cd ../plotter && rake


# ----- General rules, top down

# MDD-style comparison report:

SUMMARY=work/$(shell basename $A)-$(shell basename $B)-summary.txt
REPORT=work/$(shell basename $A)-$(shell basename $B)-report.csv
MERGED=work/$(shell basename $A)-$(shell basename $B)-merged.csv
MATCHES=work/$(shell basename $A)-$(shell basename $B)-matches.csv
ROUND=work/$(shell basename $A)-$(shell basename $B)-round.csv
DELTA=work/$(shell basename $A)-$(shell basename $B)-delta.csv

report: $(REPORT)
REPORT_OPTIONS?=

$(REPORT): $(MERGED) $P/report.py $P/property.py
	@echo
	@echo "--- PREPARING REPORT ---"
	$P/report.py $(REPORT_OPTIONS) --summary $(SUMMARY) < $(MERGED) > $@.new
	@mv -f $@.new $@

# Merged checklist on which the report is based:

merged: $(MERGED)

$(MERGED): $(MATCHES) $P/merge.py $P/theory.py $P/workspace.py $P/checklist.py $P/property.py
	@echo
	@echo "--- MERGING ---"
	$P/merge.py --A $A.csv --B $B.csv --matches $(MATCHES) \
		    > $@.new
	@mv -f $@.new $@

# Record matches, required for merging:

matches: $(MATCHES)

$(MATCHES): $A.csv $B.csv $P/match_records.py
	@echo
	@echo "--- COMPUTING RECORD MATCHES ---"
	$P/match_records.py --A $A.csv --B $B.csv --pk $(DELTA_KEY) --index $(INDEX) \
		    	    > $@.new
	@mv -f $@.new $@

# Round trip, for test of record based diff/patch: (EOL demo. no
# hierarchy sensitivity)

round: $(ROUND)

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
	$P/project.py --keep $(KEEP) <$< | \
	$P/sortcsv.py --key $(PRIMARY_KEY) >$@.new
	@mv -f $@.new $@


# Delta, describing how to change current database state into new
# state:

delta: $(DELTA)

# Formerly: $P/project.py --keep $(KEEP) <$< | ...
# and	    $P/sortcsv.py --key $(PRIMARY_KEY) <$< >$@.new
# Files need to be sorted for delta

# Columns that use to decide whether a 'record has changed'
MANAGE?=taxonID,scientificName,canonicalName,taxonRank,taxonomicStatus,nomenclaturalStatus,datasetID,source

$(DELTA): $A.csv $B.csv $P/delta.py $P/match_records.py $P/property.py
	@echo
	@echo "--- COMPUTING DELTA ---"
	set -o pipefail; \
	$P/delta.py --A $A.csv --B $B.csv --pk $(DELTA_KEY) \
		    --index $(INDEX) --manage $(KEEP)
	| $P/sortcsv.py --key $(DELTA_KEY) \
	> $@.new
	@mv -f $@.new $@
	wc $@

# something:
# 	$P/scatter.py --dest $(basename $(DELTA)) < $(DELTA)

# ----- Mammals rules

# ncbi 40674, gbif 359, dh1 EOL-000000627548, dh0.9 -168130, page 1642

TAXON=Mammalia
taxon=mammals
TAXON_NCBI=$(TAXON)
TAXON_GBIF=$(TAXON)
TAXON_DH1=$(TAXON)
TAXON_DH09=$(TAXON)

ncbi%-$(taxon)-raw.csv: ncbi%-raw.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(TAXON_NCBI) < $< > $@.new
	@mv -f $@.new $@
.PRECIOUS: ncbi%-$(taxon)-raw.csv

gbif%-$(taxon)-raw.csv: gbif%-raw.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(TAXON_GBIF) < $< > $@.new
	@mv -f $@.new $@
.PRECIOUS: gbif%-$(taxon)-raw.csv

dh1%-$(taxon)-raw.csv: dh1%-raw.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(TAXON_DH1) < $< > $@.new
	@mv -f $@.new $@
.PRECIOUS: dh1%-$(taxon)-raw.csv

dh0%-$(taxon)-raw.csv: dh1%-raw.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(TAXON_DH09) < $< > $@.new
	@mv -f $@.new $@
.PRECIOUS: dh0%-$(taxon)-raw.csv


# ----- 1. NCBI-specific rules

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
ncbi%-raw.csv: work/ncbi%/dump/names.dmp src/ncbi_to_dwc.py $P/start.py
	$P/ncbi_to_dwc.py `dirname $<` \
	| $P/start.py --pk $(PRIMARY_KEY) \
	| $P/sortcsv.py --key $(PRIMARY_KEY) > $@.new
	@mv -f $@.new $@
.PRECIOUS: ncbi%-raw.csv

# Adjoin 'type' and 'year' columns.  To skip this change the rule to
# simply copy %-raw.csv to %.csv above. 
%.csv: %-raw.csv $P/extract_names.py $P/use_parse.py
	$P/extract_names.py < $< \
	| gnparser -s \
	| $P/use_parse.py --source $< > $@.new
	@mv -f $@.new $@

# ----- 2. GBIF-specific rules

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

# Ingest GBIF dump, convert TSV to CSV, add managed_id column
gbif%-raw.csv: gbif%/dump/meta.xml $P/start.py
	$P/start.py --pk $(PRIMARY_KEY) --input `dirname $<`/Taxon.tsv \
	  --managed gbif:taxonID >$@.new
	@mv -f $@.new $@
.PRECIOUS: gbif%-raw.csv


# ----- 3. ASU/BioKIC example

# MDD

MDDSOURCE?=$(HOME)/Downloads/MDD_DwC_versions
#MDDSOURCE?=$(HOME)/Downloads/MDD_DwC_versions.20211222

work/mdd1.0-source.csv: $(MDDSOURCE)/MDD_v1_6495species_JMamm_inDwC.csv
	cp -p $< $@
work/mdd1.1-source.csv: $(MDDSOURCE)/MDD_v1.1_6526species_inDwC.csv
	cp -p $< $@
work/mdd1.2-source.csv: $(MDDSOURCE)/MDD_v1.2_6485species_inDwC.csv
	cp -p $< $@
work/mdd1.3-source.csv: $(MDDSOURCE)/MDD_v1.3_6513species_inDwC.csv
	cp -p $< $@
work/mdd1.31-source.csv: $(MDDSOURCE)/MDD_v1.31_6513species_inDwC.csv
	cp -p $< $@
work/mdd1.4-source.csv: $(MDDSOURCE)/MDD_v1.4_6533species_inDwC.csv
	cp -p $< $@
work/mdd1.5-source.csv: $(MDDSOURCE)/MDD_v1.5_6554species_inDwC.csv
	cp -p $< $@
work/mdd1.6-source.csv: $(MDDSOURCE)/MDD_v1.6_6557species_inDwC.csv
	cp -p $< $@
work/mdd1.7-source.csv: $(MDDSOURCE)/MDD_v1.7_6567species_inDwC.csv
	cp -p $< $@

mdd%-raw.csv: work/mdd%-source.csv $P/start.py
	$P/start.py < $< --pk taxonID --managed mdd:taxonID > $@.new
	@mv -f $@.new $@
.PRECIOUS: mdd%-raw.csv

# ----- 4. EOL examples

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

# DH 1.2 hasn't yet been 'harvested' or 'published' in EOL, so we have to
# get it straight from opendata
DH12_LP="https://opendata.eol.org/dataset/tram-807-808-809-810-dh-v1-1/resource/02037fde-cc69-4f03-94b5-65591c6e7b3b"

work/dh12.taxafilename:
	@mkdir -p work
	echo `$(RAKE) dwca:taxa_path OPENDATA=$(DH12_LP)` >$@.new
	@mv -f $@.new $@

%.taxafilename: %.eol-resource-id
	@mkdir -p work
	@echo Resource id is $$(cat $<)
	ID=$$(cat $<); \
	  echo `$(RAKE) resource:taxa_path CONF=$(ASSEMBLY) ID=$$ID` >$@.new
	@mv -f $@.new $@
	@echo Taxa path is `cat $@`

# about half a minute for DH 1.1
# the managed_id can only be set if DH 0.9 has its records mapped to
# pages (see -mapped)

work/dh09-raw.csv: work/dh09.taxafilename $P/start.py
	@mkdir -p work
	$P/start.py --input `cat $<` \
                    --managed eol:EOLid \
		    --pk taxonID \
	       > $@.new
	@mv -f $@.new $@

# taxonID is managed in 1.1 and following, but not in 0.9

dh1%-raw.csv: dh1%.taxafilename $P/start.py
	@mkdir -p work
	cat $<
	$P/start.py --input `cat $<` \
                    --managed eolnode:taxonID \
		    --pk taxonID \
	| $P/sortcsv.py --key taxonID > $@.new
	@mv -f $@.new $@

# in1=./deprecated/work/1-mam.csv
# in2=./deprecated/work/724-mam.csv

work/dh11-$(taxon)-hier.csv: work/dh11-$(taxon).csv work/dh11-map-raw.csv $P/hierarchy.py
	$P/hierarchy.py --mapping work/dh11-map.csv \
		  < $< \
		  > $@.new
	@mv -f $@.new $@

# EOL dynamic hierarchy - usages mapped to pages

work/%-hier.csv: work/%-raw.csv work/%-map.csv $P/hierarchy.py
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

work/dh12-map.csv: work/dh11-map-raw.csv
	cp $< $@

# Deprecated ... ?

work/%-mapped.csv: work/%-raw.csv work/%-map.csv $P/idmap.py
	$P/idmap.py --mapping $(basename $<)-map.csv) \
		  < $< > $@.new
	@mv -f $@.new $@

# -----

tags:
	etags $P/*.py

test:
	$P/property.py
	$P/checklist.py
	$P/workspace.py
	$P/merge.py --test
	$(MAKE) A=work/test1 B=work/test2 report

test-report: 
	$(MAKE) A=work/test1 B=work/test2 report
	@echo 'A: $(TEST1)'
	@echo 'B: $(TEST2)'
	$P/newick.py < work/test1-test2-merged.csv

TEST1="((a,b)d,c)f1"
TEST2="(a,(b,c)e)f2"

work/test1-raw.csv:
	$P/newick.py --newick $(TEST1)  >$@

work/test2-raw.csv:
	$P/newick.py --newick $(TEST2) >$@
