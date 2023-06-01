# $< = first prerequisite, $@ = target

# See lower part of file for examples (NCBI etc.)

# This has examples of the use of the list tools.
# Very much in flux!

all:
	@echo "Use make parameter syntax to specify the two checklists, e.g."
	@echo "  make A=work/ncbi201505-mammals B=work/ncbi202008-mammals demo"

outofdate:
	@echo "Please specify a target:"
	@echo "  demo       align and generate Euler/X form of alignment"
	@echo "  diff       show new/removed/changed records"
	@echo "  round      demo diff/patch round trip"
	@echo "  report     report on NCBI extensions differences 2015-2020"
	@echo "  mdd-report report on MDD extensions differences 1.6-1.7"

.SECONDARY:

# ----- General parameters

SHELL?=/usr/bin/bash
P?=src

# Primary key column
PRIMARY_KEY?=taxonID

# overridden to EOLid for EOL
DELTA_KEY?=$(PRIMARY_KEY)

# projection, to make the files smaller ... ?
KEEP?="taxonID,canonicalName,scientificName,tipe,canonicalStem,managed_id,parentNameUsageID,acceptedNameUsageID,taxanomicStatus,nomenclaturalStatus"

# not sure exactly what this means.  Add EOLid for EOL
MANAGED?=$(KEEP)

# This is for EOL.  Requires clone of 'plotter' repo
RAKE?=cd ../plotter && rake

# ----- General rules, top down

A=test1
B=test2
ANAME?=A
BNAME?=B
TAXON?=Mammalia
taxon?=mammals

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
SHORT=work/$(shell basename $A)-$(shell basename $B)-short.csv
EXEMPLARS=work/$(shell basename $A)-$(shell basename $B)-exemplars.csv

CODE=$P/demo.py $P/align.py $P/theory.py $P/simple.py $P/workspace.py \
     $P/checklist.py $P/rcc5.py $P/eulerx.py $P/linkage.py $P/lub.py \
     $P/exemplar.py $P/parse.py $P/util.py

demo: $(DEMO)

$(DEMO): $(CODE) $A.csv $B.csv
	@echo
	@echo "--- PREPARING DEMO ---"
	$P/demo.py --A $A.csv --B $B.csv --Aname $(ANAME) --Bname $(BNAME) \
	  --eulerx $(EULERX).new --short $(SHORT).new --long $@.new
	@mv -f $@.new $@
	@mv -f $(SHORT).new $(SHORT)
	@mv -f $(EULERX).new $(EULERX)

# Record matches, required for merging:

matches: $(MATCHES)

$(MATCHES): $A.csv $B.csv $P/match_records.py
	@echo
	@echo "--- COMPUTING RECORD MATCHES ---"
	$P/match_records.py --A $A.csv --B $B.csv --pk $(DELTA_KEY) \
		    	    > $@.new
	@mv -f $@.new $@

# Stale targets, relicts of some EOL work

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
		    --manage $(KEEP)
	| $P/sortcsv.py --key $(DELTA_KEY) \
	> $@.new
	@mv -f $@.new $@
	wc $@

# something:
# 	$P/scatter.py --dest $(basename $(DELTA)) < $(DELTA)

# ----- Generally useful little file transformation rules

# Download a DwCA zip file
in/%.zip: work/%.dwca-url
	wget -O $@.new $$(cat $<)
	touch $@.new
	mv -f $@.new $@

# Extract contents of a .zip file

work/%.dump: in/%.zip
	mkdir -p $@.new
	unzip -d $@.new $<
	rm -rf $@ && mv $@.new $@
	touch $@/*

# Normalize the DwCA taxon file (no managed id space).

%-raw.csv: %.dump $P/start.py
	$P/start.py --pk $(PRIMARY_KEY) --input `src/find_taxa.py $<` \
	  >$@.new
	@mv -f $@.new $@

# Adjoin 'tipe' and 'year' columns.  To skip this step, change the rule to
# cp -pf $< $@
%.csv: %-raw.csv $P/extract_names.py $P/use_gnparse.py
	$P/extract_names.py < $< \
	| gnparser -s \
	| $P/use_gnparse.py --source $< > $@.new
	@mv -f $@.new $@

# Subsetting:

%-$(taxon)-raw.csv: %-raw.csv $P/subset.py
	$P/subset.py < $< --hierarchy $< --root $(TAXON) > $@.new
	@mv -f $@.new $@

# Make the file smaller by eliminating unneeded columns

%-narrow.csv: %.csv $P/project.py
	$P/project.py --keep $(KEEP) <$< | \
	$P/sortcsv.py --key $(PRIMARY_KEY) >$@.new
	@mv -f $@.new $@

tags:
	etags $P/*.py

# Testing

test-demo: work/test1.csv work/test2.csv
	$(MAKE) A=work/test1 B=work/test2 demo
work/test1.csv: work/test1-raw.csv $P/extract_names.py $P/use_gnparse.py
work/test2.csv: work/test2-raw.csv $P/extract_names.py $P/use_gnparse.py

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

# -----------------------------------------------------------------------------
# Examples

# N.b. managed_id has a specific meaning in this system; it's an id in
#  some managed namespace, prefixed by something to identify the namespace.
# Example syntax:
#   ncbi:123, eol:4567 (for page id), gbif:8910, worms:2345   etc.

# ----- 1. NCBI example:

ncbi-demo:
	$(MAKE) A=work/ncbi201505-$(taxon) B=work/ncbi202008-$(taxon) \
	  ANAME=N2015 BNAME=N2020 demo

work/ncbi201505-$(taxon).csv: work/ncbi201505-$(taxon)-raw.csv \
   $P/extract_names.py $P/use_gnparse.py
work/ncbi202008-$(taxon).csv: work/ncbi202008-$(taxon)-raw.csv \
   $P/extract_names.py $P/use_gnparse.py
# raw-to-raw subsetting is implicit... 

work/ncbi201505-$(taxon)-raw.csv: work/ncbi201505-raw.csv
work/ncbi202008-$(taxon)-raw.csv: work/ncbi202008-raw.csv

# Convert NCBI taxdump to DwC form
work/ncbi%-raw.csv: work/ncbi%.dump/names.dmp $P/ncbi_to_dwc.py $P/start.py
	$P/ncbi_to_dwc.py `dirname $<` \
	| $P/start.py --pk $(PRIMARY_KEY) > $@.new
	@mv -f $@.new $@

# Extract files from NCBI .zip file
work/%.dump/names.dmp: in/%.zip
	rm -rf `dirname $@`
	mkdir -p `dirname $@`
	unzip -d `dirname $@` $<
	touch `dirname $@`/*

# Download and unpack some version of NCBI Taxonomy
# The theory of in/ is not well baked yet
in/%.zip: sources/%.ncbi-url
	wget -O $@.new $$(cat $<)
	mv -f $@.new $@

work/ncbi201505.ncbi-url:
	echo ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2015-05-01.zip \
	  >$@
work/ncbi202008.ncbi-url:
	echo ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2020-08-01.zip \
	  >$@

# ----- 2. GBIF examples:

gbif-demo:
	$(MAKE) A=work/gbif20190916-$(taxon) B=work/gbif20210303-$(taxon) \
	  ANAME=G2019 BNAME=G2021 demo

work/gbif20190916-$(taxon).csv: work/gbif20190916-$(taxon)-raw.csv \
   $P/extract_names.py $P/use_gnparse.py
work/gbif20210303-$(taxon).csv: work/gbif20210303-$(taxon)-raw.csv \
   $P/extract_names.py $P/use_gnparse.py
# raw-to-raw subsetting is implicit, don't need the following
#work/gbif20190916-$(taxon)-raw.csv: work/gbif20190916-raw.csv

# an instance of the DwC ingest rule.  shouldn't need
#work/gbif20190916-raw.csv: work/gbif20190916.dump

# make A=ncbi202008-mammals B=gbif20210303-mammals demo
# and so on.

# GBIF-specific rules

# Ingest GBIF dump, convert TSV to CSV, add managed_id column
work/gbif%-raw.csv: work/gbif%.dump $P/start.py
	$P/start.py --pk $(PRIMARY_KEY) --input `src/find_taxa.py $<` \
	  --managed gbif:taxonID >$@.new
	@mv -f $@.new $@
.PRECIOUS: work/gbif%-raw.csv

work/gbif20190916.dwca-url:
	echo https://rs.gbif.org/datasets/backbone/2019-09-06/backbone.zip >$@

work/gbif20210303.dwca-url:
	echo https://rs.gbif.org/datasets/backbone/2021-03-03/backbone.zip >$@

# ----- 3. BioKIC/ATCR examples:

# To obtain Darwin Core versions of MDD:
#   1. Download the main .csv files for MDD from Zenodo
#   2. Run them through Prashant's tool at
#   https://github.com/jar/MDD-DwC-mapping/ 
# Automation for all this is on the to-do list... don't know if I'll
# ever get around to it

mdd-demo:
	$(MAKE) A=work/msw3 B=work/mdd1.0 ANAME=MSW3 BNAME=MDD1_1 demo
	$(MAKE) A=work/msw3 B=work/mdd1.10 ANAME=MSW3 BNAME=MDD1_10 demo
	$(MAKE) A=work/mdd1.0 B=work/mdd1.10 ANAME=MDD1 BNAME=MDD1_10 demo

mdd-demo-67:
	$(MAKE) A=work/mdd1.6 B=work/mdd1.7 ANAME=MDD1_6 BNAME=MDD1_7 demo

mdd-demo-67p:
	$(MAKE) A=work/mdd1.6-primates B=work/mdd1.7-primates ANAME=MDD1_6 BNAME=MDD1_7 \
	  taxon=Primates TAXON=Primates demo

mdd-demo-01:
	$(MAKE) A=work/mdd1.10 B=work/mdd1.11 ANAME=MDD1_10 BNAME=MDD1_11 demo

work/mdd1.6.csv: work/mdd1.6-raw.csv
work/mdd1.7.csv: work/mdd1.7-raw.csv
work/mdd1.6-primates-raw.csv: work/mdd1.6-raw.csv
work/mdd1.7-primates-raw.csv: work/mdd1.7-raw.csv

# make A=mdd1.2-mammals B=mdd1.3 report
# make A=mdd1.2 B=mdd1.3 report
# make A=gbif20210303-mammals B=mdd1.0-mammals report
# Prashant's request, see slack on 5/2/2022:
# make A=mdd1.7 B=gbif20220317-mammals report


# Norway = artsnavnebase

hyg-demo:
	$(MAKE) A=work/nor-hyg B=work/swe-hyg demo

work/nor-raw.csv: work/nor.dump $P/start.py

# Lost track of the URL.  Go do artsnavnebase web site and look around
work/nor.dwca-url:
	echo "http://mumble.net/something/something/dwca-artsnavnebase-v1.128.zip" >$@

# Sweden = dyntaxa = artdatabanken

work/swe-hyg.csv: work/swe-hyg-raw.csv

work/swe-hyg-raw.csv: work/swe-raw.csv $P/subset.py
	$P/subset.py < $< --hierarchy $< --root Hygrophorus > $@

work/swe-raw.csv: work/swe.dump $P/start.py
	$P/start.py --pk $(PRIMARY_KEY) --input `src/find_taxa.py $<` | \
	  grep -v ",speciesAggregate," \
	  >$@.new
	@mv -f $@.new $@


# This doesn't work - requires an API token...
#   --header="Ocp-Apim-Subscription-Key: a300a2..etc..etc"
work/swe.dwca-url:
	echo "https://api.artdatabanken.se/taxonservice/v1/darwincore/download" >$@



# ----- 4. EOL examples:

# Requires clone of 'plotter' repo.

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

col-demo:
	$(MAKE) A=col2019-$(taxon) B=col2021-$(taxon) demo

# make A=col2021-mammals B=mdd1.7 demo
# and so on.

# ----- 3. ASU/BioKIC example

# Sources are in pgasu/MDD-DwC-mapping repo, based on original sources
# on zenodo
#  https://zenodo.org/record/7394529/files/MDD_v1.10_6615species.csv?download=1
#  https://zenodo.org/record/4139723/files/MDD_v1_6495species_JMamm.csv?download=1

in/m: sources/mdd1.10.url
	wget -O m.csv $$(cat $<)

# MDD

# 1.0 and 1.1 don't use the later managed ids
work/mdd1.0.csv: work/mdd1.0-raw.csv \
   $P/extract_names.py $P/use_gnparse.py
work/mdd1.0-raw.csv: work/mdd/mdd1.0.csv
	$P/start.py < $< --pk taxonID > $@.new
	@mv -f $@.new $@

# Need to clone the pgasu/MDD-DwC-mapping repo and put the clone sister to this repo
# Get later versions at https://zenodo.org/record/7394529#.Y-z1dOLMI1I
MAPPER?=../MDD-DwC-mapping
MDDSOURCE?=$(MAPPER)/data
MDDDWC?=$(MAPPER)/dwc
$(MDDDWC)/mdd1.0-dwc.csv: $(MDDSOURCE)/MDD_v1_6495species_JMamm.csv
$(MDDDWC)/mdd1.1-dwc.csv: $(MDDSOURCE)/MDD_v1.1_6526species.csv
$(MDDDWC)/mdd1.2-dwc.csv: $(MDDSOURCE)/MDD_v1.2_6485species.csv
$(MDDDWC)/mdd1.3-dwc.csv: $(MDDSOURCE)/MDD_v1.3_6513species.csv
$(MDDDWC)/mdd1.31-dwc.csv: $(MDDSOURCE)/MDD_v1.31_6513species.csv
$(MDDDWC)/mdd1.4-dwc.csv: $(MDDSOURCE)/MDD_v1.4_6533species.csv
$(MDDDWC)/mdd1.5-dwc.csv: $(MDDSOURCE)/MDD_v1.5_6554species.csv
$(MDDDWC)/mdd1.6-dwc.csv: $(MDDSOURCE)/MDD_v1.6_6557species.csv
$(MDDDWC)/mdd1.7-dwc.csv: $(MDDSOURCE)/MDD_v1.7_6567species.csv
$(MDDDWC)/mdd1.8-dwc.csv: $(MDDSOURCE)/MDD_v1.8_6591species.csv
$(MDDDWC)/mdd1.9-dwc.csv: $(MDDSOURCE)/MDD_v1.9_6596species.csv
$(MDDDWC)/mdd1.10-dwc.csv: $(MDDSOURCE)/MDD_v1.10_6615species.csv
$(MDDDWC)/mdd1.11-dwc.csv: $(MDDSOURCE)/MDD_v1.11_6649species.csv

work/mdd%-raw.csv: $(MDDDWC)/mdd%-dwc.csv $P/start.py
	$P/start.py < $< --pk taxonID > $@.new
	@mv -f $@.new $@

work/mdd1.0-raw.csv: $(MDDDWC)/mdd1.0-dwc.csv $P/start.py
work/mdd1.1-raw.csv: $(MDDDWC)/mdd1.1-dwc.csv $P/start.py
work/mdd1.6-raw.csv: $(MDDDWC)/mdd1.6-dwc.csv $P/start.py
work/mdd1.7-raw.csv: $(MDDDWC)/mdd1.7-dwc.csv $P/start.py
work/mdd1.10-raw.csv: $(MDDDWC)/mdd1.10-dwc.csv $P/start.py
work/mdd1.11-raw.csv: $(MDDDWC)/mdd1.11-dwc.csv $P/start.py

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
col%-raw.csv: col%.dump $P/start.py
	$P/start.py --pk $(PRIMARY_KEY) --input `src/find_taxa.py $<` \
	  >$@.new
	@mv -f $@.new $@
.PRECIOUS: col%-raw.csv

work/col2023-mammals-raw.csv: work/col2023-raw.csv
	$P/subset.py < $< --hierarchy $< --root "6224G" > $@.new
	@mv -f $@.new $@

work/col2021-mammals-raw.csv: work/col2021-raw.csv
	$P/subset.py < $< --hierarchy $< --root "6224G" > $@.new
	@mv -f $@.new $@

# ----- 6. ITIS
# Where do we get ITIS?  Need a subset step.

work/itis2022-mammals-raw.csv: work/itis2022.dump
work/itis2022-mammals.csv: work/itis2022-mammals-raw.csv

# ----- 7. MDD and other

# Where did I get the file sources/msw3-raw.csv ?
# How was it created?

work/msw3.csv: work/msw3-raw.csv \
   $P/extract_names.py $P/use_gnparse.py

work/msw3-raw.csv: sources/msw3-source.csv $P/start.py
	$P/start.py --pk $(PRIMARY_KEY) --input $< \
	  > $@.new
	@mv -f $@.new $@

# Markus's use case from checklistbank

# Sweden = dyntaxa = artdatabanken
# https://www.checklistbank.org/dataset/2041/download
# It would be better to get this from artdatabanken rather than GBIF

work/dyntaxa%-raw.csv: in/dyntaxa%.tsv $P/start.py
	$P/start.py --pk $(PRIMARY_KEY) --input $< \
	> $@.new
	@mv -f $@.new $@

# Norway = artsnavebasen
# https://www.checklistbank.org/dataset/2030/download

work/arts%-raw.csv: in/arts%.tsv $P/start.py
	$P/start.py --pk $(PRIMARY_KEY) --input $< \
	> $@.new
	@mv -f $@.new $@

# -----

