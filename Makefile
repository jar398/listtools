# $< = first prerequisite, $@ = target

# This has examples of the use of the list tools.
# Very much in flux!

all:
	@echo "Please specify a target:"
	@echo "  demo       align and generate Euler/X form of alignment"
	@echo "  diff       show new/removed/changed records"
	@echo "  round      demo diff/patch round trip"
	@echo "  report     report on NCBI extensions differences 2015-2020"
	@echo "  mdd-report report on MDD extensions differences 1.6-1.7"
	@echo "Use make parameter syntax to specify the two checklists, e.g."
	@echo "  make A=work/ncbi201505-mammals B=work/ncbi202008-mammals report"

# Default example:
# Compare two versions of NCBI mammals

TAXON?=Mammalia
taxon?=mammals
A?=work/ncbi201505-$(taxon)
B?=work/ncbi202008-$(taxon)
ANAME?=A
BNAME?=B

# make A=ncbi201505-mammals B=ncbi202008-mammals demo
# time make A=work/ncbi201505 B=work/ncbi202008 demo
# time make A=work/ncbi201505 B=work/ncbi202008 demo
#  etc. etc.

# N.b. managed_id has a specific meaning in this system; it's an id in
#  some managed namespace, prefixed by something to identify the namespace.
# Example syntax:
#   ncbi:123, eol:4567 (for page id), gbif:8910, worms:2345   etc.

# ----- 1. NCBI example:

ncbi-report:
	$(MAKE) A=work/ncbi201505-$(taxon) B=work/ncbi202008-$(taxon) report

# make A=ncbi201505 B=ncbi202008 report  # big, will take forever

# ----- 2. GBIF examples:

gbif-report:
	$(MAKE) A=work/gbif20190916-$(taxon) B=work/gbif20210303-$(taxon) report

# make A=ncbi202008-mammals B=gbif20210303-mammals report
# and so on.

# ----- 3. BioKIC/ATCR examples:

mdd-demo:
	$(MAKE) A=work/mdd1.6 B=work/mdd1.7 ANAME=MDD1_6 BNAME=MDD1_7 demo

# make A=mdd1.2-mammals B=mdd1.3 report
# make A=mdd1.2 B=mdd1.3 report
# make A=gbif20210303-mammals B=mdd1.0-mammals report
# Prashant's request, see slack on 5/2/2022:
# make A=mdd1.7 B=gbif20220317-mammals report

# ----- 4. EOL examples:

# Requires clone of 'plotter' repo:

eol-report:
	$(MAKE) A=dh11-$(taxon) B=dh12-$(taxon) report

# time make A=dh11-mammals B=dh12-mammals round
# time make A=dh09-mammals B=dh11-mammals round
# time make A=dh11 B=dh12 round
# time make A=dh09 B=dh11 round

# Hierarchies - columns are those that neo4j needs to know
# DELTA_KEY=EOLid MANAGE=EOLid,parentEOLid,taxonID,landmark_status \
#   time make A=dh09-hier B=dh11-hier report

# DELTA_KEY=EOLid MANAGE=EOLid,parentEOLid,taxonID,landmark_status \
#   time make A=dh11-hier B=dh12-hier report

# ----- 5. CoL examples:

col-report:
	$(MAKE) A=col2019-$(taxon) B=col2021-$(taxon) report

# make A=col2021-mammals B=mdd1.7 report
# and so on.

# ----- 6. GBIF/MSW/MDD:

msw-demo:
	$(MAKE) A=msw3-$(taxon) B=mdd1.9-$(taxon) demo

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

SUMMARY=work/$(shell basename $A)-$(shell basename $B)-summary.csv
REPORT=work/$(shell basename $A)-$(shell basename $B)-report.csv
MERGED=work/$(shell basename $A)-$(shell basename $B)-merged.csv
ALIGNEDX=work/$(shell basename $A)-$(shell basename $B)-alignedx.csv
MATCHES=work/$(shell basename $A)-$(shell basename $B)-matches.csv
ROUND=work/$(shell basename $A)-$(shell basename $B)-round.csv
DELTA=work/$(shell basename $A)-$(shell basename $B)-delta.csv

DEMO=work/$(shell basename $A)-$(shell basename $B)-aligned.csv
EULERX=work/$(shell basename $A)-$(shell basename $B)-eulerx.txt
SHORT=work/$(shell basename $A)-$(shell basename $B)-short.txt
TIPWARDS=work/$(shell basename $A)-$(shell basename $B)-tipwards.csv

demo: $(DEMO)

$(DEMO): $P/demo.py $P/checklist.py $P/align.py $P/theory.py $P/workspace.py \
	   $A.csv $B.csv
	@echo
	@echo "--- PREPARING DEMO ---"
	$P/demo.py --A $A.csv --B $B.csv --Aname $(ANAME) --Bname $(BNAME) \
	  --eulerx $(EULERX).new --short $(SHORT).new > $@.new
	@mv -f $@.new $@
	@mv -f $(EULERX).new $(EULERX)
	@mv -f $(SHORT).new $(SHORT)

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

# Alignment of two checklists:

aligned: $(ALIGNED)

$(ALIGNED): $(MATCHES) $P/align.py $P/theory.py $P/workspace.py $P/checklist.py $P/property.py
	@echo
	@echo "--- ALIGNING ---"
	$P/align.py --A $A.csv --B $B.csv --matches $(MATCHES) \
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

# Any Darwin Core archive (DwCA):

%/dump/meta.xml: %.dwca-url
	mkdir -p `dirname $@`/dump
	wget -O `dirname $@`.zip $$(cat $<)
	unzip -d `dirname $@` `dirname $@`.zip
	touch `dirname $@`/*
.PRECIOUS: %/dump/meta.xml

# ----- Mammals rules

# ncbi 40674, gbif 359, dh1 EOL-000000627548, dh0.9 -168130, page 1642
# Most of these rules may be redundant now.  Cull at some point.

work/ncbi%-$(taxon)-raw.csv: in/ncbi%-raw.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(TAXON) < $< > $@.new
	@mv -f $@.new $@
.PRECIOUS: ncbi%-$(taxon)-raw.csv

work/gbif%-$(taxon)-raw.csv: in/gbif%-raw.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(TAXON) < $< > $@.new
	@mv -f $@.new $@
.PRECIOUS: gbif%-$(taxon)-raw.csv

work/dh1%-$(taxon)-raw.csv: in/dh1%-raw.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(TAXON) < $< > $@.new
	@mv -f $@.new $@
.PRECIOUS: dh1%-$(taxon)-raw.csv

work/dh0%-$(taxon)-raw.csv: in/dh0%-raw.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(TAXON) < $< > $@.new
	@mv -f $@.new $@
.PRECIOUS: dh0%-$(taxon)-raw.csv

work/col2021-mammals-raw.csv: in/col2021-raw.csv $P/subset.py
	$P/subset.py --hierarchy $< --root 6224G < $< > $@.new
	@mv -f $@.new $@
work/col2021-primates-raw.csv: in/col2021-raw.csv $P/subset.py
	$P/subset.py --hierarchy $< --root 3W7 < $< > $@.new
	@mv -f $@.new $@
work/%-$(taxon)-raw.csv: in/%-raw.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(TAXON) < $< > $@.new
	@mv -f $@.new $@
.PRECIOUS: work/%-$(taxon)-raw.csv

work/mdd1.7-$(taxon)-raw.csv: work/mdd1.7-raw.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(TAXON) < $< > $@.new
	@mv -f $@.new $@
work/mdd1.6-$(taxon)-raw.csv: work/mdd1.6-raw.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(TAXON) < $< > $@.new
	@mv -f $@.new $@
# Don't know why this doesn't work
work/mdd1.0-$(taxon)-raw.csv: work/mdd1.0-raw.csv $P/subset.py
	$P/subset.py --hierarchy $< --root $(TAXON) < $< > $@.new
	@mv -f $@.new $@

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

# Markus's use case from checklistbank
work/dyntaxa%-raw.csv: in/dyntaxa%.tsv $P/start.py
	$P/start.py --pk $(PRIMARY_KEY) --input $< \
	> $@.new
	@mv -f $@.new $@

work/arts%-raw.csv: in/arts%.tsv $P/start.py
	$P/start.py --pk $(PRIMARY_KEY) --input $< \
	> $@.new
	@mv -f $@.new $@

# Adjoin 'type' and 'year' columns.  To skip this change the rule to
# simply copy %-raw.csv to %.csv above. 
%.csv: %-raw.csv $P/extract_names.py $P/use_gnparse.py
	$P/extract_names.py < $< \
	| gnparser -s \
	| $P/use_gnparse.py --source $< > $@.new
	@mv -f $@.new $@

# ----- 2. GBIF-specific rules

work/gbif20190916.dwca-url:
	echo https://rs.gbif.org/datasets/backbone/2019-09-06/backbone.zip >$@

work/gbif20210303.dwca-url:
	echo https://rs.gbif.org/datasets/backbone/2021-03-03/backbone.zip >$@

# Ingest GBIF dump, convert TSV to CSV, add managed_id column
gbif%-raw.csv: gbif%/dump/meta.xml $P/start.py
	$P/start.py --pk $(PRIMARY_KEY) --input `src/find_taxa.py $<` \
	  --managed gbif:taxonID >$@.new
	@mv -f $@.new $@
.PRECIOUS: gbif%-raw.csv

# ----- 3. ASU/BioKIC example

# MDD

# Need to clone the pgasu/MDD-DwC-mapping repo and put the clone sister to this repo
# Get later versions at https://zenodo.org/record/7394529#.Y-z1dOLMI1I
MDDSOURCE?=../MDD-DwC-mapping/data
CONVERTMDD=mkdir -p work/mdd && python3 ../MDD-DwC-mapping/src/explore_data.py

work/mdd/mdd1.0.csv: $(MDDSOURCE)/MDD_v1_6495species_JMamm.csv
	$(CONVERTMDD) --input $< --output $@

work/mdd/mdd1.1.csv: $(MDDSOURCE)/MDD_v1.1_6526species.csv
	$(CONVERTMDD) --input $< --output $@
work/mdd/mdd1.2.csv: $(MDDSOURCE)/MDD_v1.2_6485species.csv
	$(CONVERTMDD) --input $< --output $@
work/mdd/mdd1.3.csv: $(MDDSOURCE)/MDD_v1.3_6513species.csv
	$(CONVERTMDD) --input $< --output $@
work/mdd/mdd1.31.csv: $(MDDSOURCE)/MDD_v1.31_6513species.csv
	$(CONVERTMDD) --input $< --output $@
work/mdd/mdd1.4.csv: $(MDDSOURCE)/MDD_v1.4_6533species.csv
	$(CONVERTMDD) --input $< --output $@
work/mdd/mdd1.5.csv: $(MDDSOURCE)/MDD_v1.5_6554species.csv
	$(CONVERTMDD) --input $< --output $@
work/mdd/mdd1.6.csv: $(MDDSOURCE)/MDD_v1.6_6557species.csv
	$(CONVERTMDD) --input $< --output $@
work/mdd/mdd1.7.csv: $(MDDSOURCE)/MDD_v1.7_6567species.csv
	$(CONVERTMDD) --input $< --output $@
work/mdd/mdd1.8.csv: $(MDDSOURCE)/MDD_v1.8_6591species.csv
	$(CONVERTMDD) --input $< --output $@
work/mdd/mdd1.9.csv: $(MDDSOURCE)/MDD_v1.9_6596species.csv
	$(CONVERTMDD) --input $< --output $@
work/mdd/mdd1.10.csv: $(MDDSOURCE)/MDD_v1.10_6615species.csv
	$(CONVERTMDD) --input $< --output $@

work/mdd%-raw.csv: work/mdd/mdd%.csv $P/start.py
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

# ----- 5. CoL

work/col2021.dwca-url:
	echo https://download.catalogueoflife.org/col/annual/2021_dwca.zip >$@

work/col2019.dwca-url:
	echo https://download.catalogueoflife.org/col/annual/2019_dwca.zip >$@

# Ingest GBIF dump, convert TSV to CSV, add managed_id column
# If CoL had record ids we could do --managed col:taxonID 
col%-raw.csv: col%/dump/meta.xml $P/start.py
	$P/start.py --pk $(PRIMARY_KEY) --input `src/find_taxa.py $<` \
	  >$@.new
	@mv -f $@.new $@
.PRECIOUS: col%-raw.csv

# ----- 6. ITIS

work/itis2022-mammals-raw.csv: work/itis2022-mammals/dump/meta.xml
	$P/start.py --pk $(PRIMARY_KEY) --input `dirname $<`/taxa_*.txt \
	  >$@.new
	@mv -f $@.new $@

# -----

tags:
	etags $P/*.py

test:
	$P/property.py
	$P/checklist.py
	$P/workspace.py
	$P/align.py --test
	$(MAKE) A=test1 B=test2 report

test-report: 
	$(MAKE) A=test1 B=test2 report
	@echo 'A: $(TEST1)'
	@echo 'B: $(TEST2)'
	$P/newick.py < work/test1-test2-merged.csv

TEST1="((a,b)d,c)f1"
TEST2="(a,(b,c)e)f2"

work/test1-raw.csv:
	$P/newick.py --newick $(TEST1)  >$@

work/test2-raw.csv:
	$P/newick.py --newick $(TEST2) >$@
