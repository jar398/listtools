# List tools

# 'make' utilities.  The intent is that this gile gets 'included' in
# Makefiles using the 'include' command.
# N.b. file references in this makefile are to files in subdirectories
# of the listtools root directory.

# Prerequisite (used by python code):
# pip3 install regex

# Parameters:
#   A  e.g. A=$W/colmam_2023
#   B  e.g. B=$W/colmam_2024

# Github repo clone, where the python code is to be found
LISTTOOLS?=~/g/listtools

# Directory in which to put temporary 'working' files
W?=work

ANAME=$(notdir $A)
BNAME=$(notdir $B)
ABW=$W/$(ANAME)_$(BNAME)

# Location of .py files
P=$(LISTTOOLS)/src

all: $(ABW)_report.csv

# Needed for exemplars list
CODE1=$P/exemplar.py $P/simple.py \
     $P/parse.py $P/typify.py $P/redundant.py $P/specimen.py \
     $P/checklist.py $P/rcc5.py $P/util.py $P/workspace.py \
     $P/cluster.py

# Needed for report
CODE2=$P/plugin.py $P/jumble.py $P/theory.py $P/estimate.py \
      $P/eulerx.py $(CODE1)

# Combining/comparing two inputs:
#   get exemplars
#   run (checklistbank) plugin to make report(s)

$(ABW)_report.csv: $(ABW)_exemplars.csv $(CODE2)
	$P/plugin.py --exemplars $< \
		   --A $A.csv --B $B.csv --Aname $(ANAME) --Bname $(BNAME) \
	  > $@.new
	@mv -f $@.new $@

$(ABW)_exemplars.csv: $A.csv $B.csv $(CODE1)
	@echo
	$P/exemplar.py --A $A.csv --B $B.csv --Aname $(ANAME) --Bname $(BNAME) \
	  > $@.new
	@mv -f $@.new $@

$A.csv: $A-clean.csv
$B.csv: $B-clean.csv

# Adjoin 'tipe' and 'year' columns.  To skip this step, change the rule to
# cp -pf $< $@
%.csv: %-clean.csv $P/extract_names.py $P/use_gnparse.py
	$P/extract_names.py < $< \
	| gnparser -s \
	| $P/use_gnparse.py --source $< > $@.new
	@mv -f $@.new $@

%-clean.csv: %.dump $P/clean.py
	mkdir -p $W
	$P/clean.py --pk taxonID --input `$P/find_taxa.py $<` \
	  >$@.new
	@mv -f $@.new $@

# Extract contents of a .zip file (e.g. colmam_2023.dump)
# -d $@.new  = optional directory to which to extract files.
# Wart: This target is a directory and can't be 'rm'd by make.

%.dump: %.zip
	mkdir -p $W
	rm -rf $@.new
	mkdir -p $@.new
	unzip -d $@.new $<
	touch $@.new/*
	rm -rf $@ && mv $@.new $@
