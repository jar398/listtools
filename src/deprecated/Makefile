# $< = first prerequisite, $@ = target

# This Makefile has general rules for invoking the list tools, as well
# as examples.  For the examples see the lower part of the file.

# This file is in flux as the current focus of development changes.
# Currently (November 2023) development is focusing on examplar.py
# and plugin.py, which do not currently have Makefile rules.  Older
# rules may not work right now.


all:
	@echo "Use make parameter syntax to specify the two checklists, e.g."
	@echo "  make A=work/ncbi201505-mammals B=work/ncbi202008-mammals demo"

outofdate:
	@echo "Please specify a target:"
	@echo "  demo       align and generate Euler/X form of alignment"
	@echo "  diff       show new/removed/changed records"
	@echo "  round      demo diff/patch round trip"
	@echo "  report     report on NCBI extensions differences 2015-2020"
	@echo "  mdd-demo   report on MDD extensions differences 1.0-1.10"
	@echo "  mdd-plugin report on MDD differences 1.0-1.10 (plugin form)"

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

CODE=$P/demo.py $P/align.py $P/theory.py $P/simple.py $P/workspace.py \
     $P/checklist.py $P/rcc5.py $P/eulerx.py $P/linkage.py $P/estimate.py \
     $P/exemplar.py $P/parse.py $P/util.py $P/typify.py

A=work/test1
B=work/test2
ANAME?=A
BNAME?=B
TAXON?=Mammalia
taxon?=mammals

include makefile-examples.mk

# ----- General rules, top down

# For checklistbank accessory

PLUGIN=work/$(shell basename $A)-$(shell basename $B)-plugin.csv
EXEMPLARS=work/$(shell basename $A)-$(shell basename $B)-exemplars.csv

# E.g. 
# make A=work/col23.1-mammals B=work/gbif20230902-mammals exemplars
# make A=work/col19-mammals B=work/col23.1-mammals plugin

exemplars: $(EXEMPLARS)
$(EXEMPLARS): $(CODE) $A.csv $B.csv $P/exemplar.py
	@echo
	$P/exemplar.py --A $A.csv --B $B.csv --Aname $(ANAME) --Bname $(BNAME) \
	  > $@.new
	@mv -f $@.new $@

plugin: $(PLUGIN)
$(PLUGIN): $(CODE) $A.csv $B.csv $(EXEMPLARS) $P/plugin.py
	@echo
	$P/plugin.py --exemplars work/$(shell basename $A)-$(shell basename $B)-exemplars.csv \
		   --A $A.csv --B $B.csv --Aname $(ANAME) --Bname $(BNAME) \
	  > $@.new
	@mv -f $@.new $@


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

eol_report: $(EOL_REPORT)
EOL_REPORT_OPTIONS?=

$(EOL_REPORT): $(MERGED) $P/eol_report.py $P/property.py
	@echo
	@echo "--- PREPARING EOL_REPORT ---"
	$P/eol_report.py $(EOL_REPORT_OPTIONS) --summary $(SUMMARY) < $(MERGED) > $@.new
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

# Extract contents of a .zip file (GBIF, NCBI, anything)

work/%.dump: in/%.zip
	rm -rf $@.new
	mkdir -p $@.new
	unzip -d $@.new $<
	touch $@.new/*
	rm -rf $@ && mv $@.new $@

# Normalize the DwCA taxon file (no managed id space).

%-clean.csv: %.dump $P/clean.py
	$P/clean.py --pk $(PRIMARY_KEY) --input `src/find_taxa.py $<` \
	  >$@.new
	@mv -f $@.new $@

# Adjoin 'tipe' and 'year' columns.  To skip this step, change the rule to
# cp -pf $< $@
%.csv: %-clean.csv $P/extract_names.py $P/use_gnparse.py
	$P/extract_names.py < $< \
	| gnparser -s \
	| $P/use_gnparse.py --source $< > $@.new
	@mv -f $@.new $@

# Subsetting:

%-$(taxon)-clean.csv: %-clean.csv $P/subset.py
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
work/test1.csv: work/test1-clean.csv $P/extract_names.py $P/use_gnparse.py
work/test2.csv: work/test2-clean.csv $P/extract_names.py $P/use_gnparse.py

test:
	$P/property.py
	$P/checklist.py
	$P/workspace.py
	$P/align.py --test
	$(MAKE) A=test1 B=test2 eol_report

test-eol_report: 
	$(MAKE) A=test1 B=test2 eol_report
	@echo 'A: $(TEST1)'
	@echo 'B: $(TEST2)'
	$P/newick.py < work/test1-test2-merged.csv

TEST1="((a,b)d,c)f1"
TEST2="(a,(b,c)e)f2"

work/test1-clean.csv:
	$P/newick.py --newick $(TEST1)  >$@

work/test2-clean.csv:
	$P/newick.py --newick $(TEST2) >$@
