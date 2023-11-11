## Checklist comparison example

Let's put zip (or tar) files downloaded from the internet in directory `in/`,
and all other files in `work/`.


First, obtain the two checklists.  E.g. from checklistbank.org, obtain
COL19 and COL23.1, so that:
 * the A checklist = COL19
 * the B checklist = COL23.1

Checklistbank lets you choose a taxon such as Mammalia for download.
If all you have is a complete archive, you can use the `src/subset.py`
command to pull out just the taxon you care about.

Put the zip files somewhere: (choose different file names if you like)
 * `in/A.zip`
 * `in/B.zip`

Prepare file `work/A.csv` for use starting with `in/A.zip`:

    mkdir -p work/A.dump
    unzip -d work/A.dump in/A.zip

    src/clean.py --pk taxonID --input `src/find_taxa.py work/A.dump` \
      --managed gbif:taxonID > work/A-clean.csv
    src/extract_names.py < work/A-clean.csv \
    | gnparser -s \
    | src/use_gnparse.py --source work/A-clean.csv > work/A.csv

[I suppose I could create a script `src/prepare.py` that bundles up
all these steps.]

Do similarly to create `work/B.csv` from `work/B.dump`.

Now compare the two checklists, producing a report:

    src/exemplar.py --A work/A.csv --B work/B.csv > work/AB-exemplars.csv
    src/plugin.py --A work/A.csv --B work/B.csv \
      --exemplars work/AB-exemplars.csv > work/AB-report.csv

See [guide.md](guide.md) for advice on how to interpret the report.
