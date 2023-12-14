  `Makefile`
assumes that `gnparser` can be located via your shell's `PATH`.  It is
easy to remove this dependency if desired, with a degradation in
effectiveness, by modifying the appropriate `Makefile` rules to elide
the `extract_names` and `use_gnparse` steps.


## The Makefile

The Makefile does not have rules for building programs, only for
running the tools.  It has a variety of canned pipelines
for creating artifacts that are interesting and illustrative.  For
example `make ncbi-report` downloads two versions of the NCBI
(Genbank) taxonomy and produces comparison reports.  It may be helpful
to consult the Makefile for examples of use of the various tools and
ideas on how to string them together.

Read the comments in the [`Makefile`](../Makefile) for other things to try.

