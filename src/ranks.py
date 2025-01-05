# Borrowed from Open Tree of Life, with no attribution since this is
# not 'creative expression'

ranks = ["top",
         "domain",
         "superkingdom",
         "kingdom",
         "subkingdom",
         "division",            # h2007
         "infrakingdom",		# worms
         "superphylum",
         "phylum",
         "subphylum",
         "infraphylum",			# worms
         "subdivision",			# worms
         "superclass",
         "class",
         "subclass",
         "infraclass",
         "subterclass",         # worms Colobognatha
         "cohort",              # NCBI Polyneoptera
         "subcohort",           # NCBI
         "superorder",
         "order",
         "suborder",
         "infraorder",
         "parvorder",
         "section",				# worms - note, there are two 'section' ranks
         "subsection",			# worms
         "superfamily",
         "family",
         "subfamily",
         "supertribe",			# worms
         "tribe",
         "subtribe",
         "genus",
         "subgenus",
         "species group",
         "species subgroup",
         "species",
         "infraspecificname",
         "subspecies",
         "natio",               # worms
         "variety",
         "varietas",
         "subvariety",
         "form",                # 2016 GBIF
         "forma",
         "subform",
         "cluster",
         ]

ranks_dict = {name: i for (name, i) in zip(ranks[::-1], range(0, len(ranks)))}
