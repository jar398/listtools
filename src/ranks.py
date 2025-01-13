# Borrowed from Open Tree of Life

# Highest numbered ranks are for 'bigger' taxa

rank_names = [
         "never zero",
         "cluster",
         "subform",
         "forma",
         "form",                # 2016 GBIF
         "subvariety",
         "varietas",
         "variety",
         "natio",               # worms
         "subspecies",
         "infraspecificname",
         "species",
         "species subgroup",
         "species group",
         "subgenus",
         "genus",
         "subtribe",
         "tribe",
         "supertribe",			# worms
         "subfamily",
         "family",
         "superfamily",
         "subsection",			# worms
         "section",				# worms - note, there are two 'section' ranks
         "parvorder",
         "infraorder",
         "suborder",
         "order",
         "superorder",
         "subcohort",           # NCBI
         "cohort",              # NCBI Polyneoptera
         "subterclass",         # worms Colobognatha
         "infraclass",
         "subclass",
         "class",
         "superclass",
         "subdivision",			# worms
         "infraphylum",			# worms
         "subphylum",
         "phylum",
         "superphylum",
         "infrakingdom",		# worms
         "division",            # h2007
         "subkingdom",
         "kingdom",
         "superkingdom",
         "domain",
         "top",
         ]

ranks_dict = {name: i for (name, i) in zip(rank_names, range(0, len(rank_names)))}
