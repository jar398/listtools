Example: J Liljeblad et al., TDWG conference 2022

This exercise looks at two fungal genera as covered in three different
checklists.  The checklists are:

 * _Hygrophorus_ Fr.
 * _Tricholoma_ (Fr.)

The checklists are:

 * Dyntaxa. Svensk taxonomisk databas 2022-08-22 (Sweden)
   Download DwC from: https://www.checklistbank.org/dataset/2041/download
 * Artsnavebasen (Norway)
   Download DwC from: https://www.checklistbank.org/dataset/2030/download
 * Catalogue of Life (COL)
   Download DwC from: https://www.checklistbank.org/dataset/3LR/download

After manual dowload (from either checklistbank or direct from the
source taxonomies using the subset tool) we have six files, which
we'll name as follows:

    in/dyntaxa-hygrophorus.tsv
    in/dyntaxa-tricholoma.tsv
    in/arts-hygrophorus.tsv
    in/arts-tricholoma.tsv
    in/col-hygrophorus.tsv
    in/col-tricholoma.tsv

Convert to CSV and do some normalizations:

    src/start.py --input in/dyntaxa-hyg.tsv > work/dyntaxa-hyg-raw.csv
    src/start.py --input in/dyntaxa-trich.tsv > work/dyntaxa-trich-raw.csv
    src/start.py --input in/arts-hyg.tsv > work/arts-hyg.csv

etc.

Create Euler/X format alignment:

    make demo A=work/dyntaxa-hyg B=work/arts-hyg ANAME=SWE BNAME=NOR

This goes through the following stages:

 1. Apply `gnparse` to the scientificNames
 1. Find name matches using scoring
 1. Figure out which matches are 'tipwards'
 1. For each tipwards match, create an 'exemplar' record via union/find
 1. Find 'cross MRCAs' linking previously unlinked records
 1. Find 'estimates' that are < or = between trees
    (usually same as cross MRCAs, just a few exceptions)
 1. Alignment generation code now has full RCC5 record comparison ability.

Enhancements: distance measure and two-pass estimation to link genus
changes and other boundary conditions.