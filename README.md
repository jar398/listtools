# listtools
Tools for manipulating lists of things, and taxonomic checklists in particular.

## NCBI worked example

It is possible to compare a clade as treated in two different versions
of NCBI Taxonomy.  Here's how.



## Configuration for EOL

For EOL, the list tools are intended to be used in conjunction with
`plotter`, which provides access to the EOL content repositories.

Configure `plotter` according to its documentation.  Then, assuming
the `plotter` repository is cloned in directory `../plotter`, do the
following:

    ln -sfT ../plotter/config config

If it's elsewhere, adjust the above command, and also update
`Makefile` where it mentions `plotter`.
