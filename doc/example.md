Checklist comparison example

Let's put zip files downloaded from the internet in directory `in/`,
and all other files in `work/`.


First, obtain the two checklists.

E.g. from checklistbank.org, obtain COL19 and COL23.1.

Put them somewhere, for example:
 * `in/col19-mammals.zip`
 * `in/col23.1-mammals.zip`

Prepare them for use.  Replace 'A' with the particular checklist to be processed.
    mkdir -p work/A.dump
    unzip -d work/A.dump in/A.zip
    src/find_taxa.py work/A.dump

    src/clean.py --pk taxonID --input `src/find_taxa.py work/gbif20230902-mammals.dump` \
      --managed gbif:taxonID >work/gbif20230902-mammals-clean.csv
