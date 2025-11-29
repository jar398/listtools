# Might require bash and/or gnu make

# Do this only once and only if you have the .zip file
# mv /home/jar/Downloads/3f3bb6d7-703a-4004-8db3-be169ae96eb5.zip col24-mammals.zip

L=~/g/listtools
P=$L/src
unzip col24-mammals.zip -d col24-mammals
$P/clean.py --input `$P/find_taxa.py col24-mammals` >col24-mammals-clean.csv
$P/extract_names.py < col24-mammals-clean.csv \
	| gnparser -s \
	| $P/use_gnparse.py --source col24-mammals-clean.csv > col24-mammals.csv
