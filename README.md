# listtools

Title: List Tools v1.alpha1

Description: Tools for manipulating lists of things, and in
particular, lists of taxonomic names (i.e. 'checklists'), which are
sometimes organized hierachically ('taxonomies').

Author: Jonathan A. Rees

* [Installation](INSTALL.md)
* [User guide](doc/guide.md)

As of January 2026, this software is located in 'jar398/listtools' on
github.com.  The location for any or all versions may change from time
to time in the future, so look around to find the latest version.

### Status of this software

The is research-grade software, by which I mean that I've done the
minimum to make it work on the examples in front of me.  I do have
some experience in software engineering, but it does not show here.
There have been frequent changes of algorithms, data structures, and
philosophy, and this fact bleeds through into the code and internal
interfaces.  Documentation is poor, there is dead code, there are no
unit tests or modules.  There is no web UI or graphical UI, just
command line (bash) or script (python3) use.

### History

The genesis for this software was a conversation around 2018 or 2019
at ASU with Nico Franz, who suggested that a generic tool for
comparing checklists might be widely useful.  Before this conversation
I had been working on
[Smasher](https://github.com/opentreeoflife/reference-taxonomy/) (see
[article](https://doi.org/10.3897/bdj.5.e12581)), which did something
similar in a different way (similar to the 'splits' of Redelings and
Holder), and i felt a desire to rewrite it on a firmer theoretical,
philosophical, and nomenclatural foundation to make it more "sound",
and to create a free-standing tool independent of the open tree of
life project or any other framework.

This tool suite replaces the older ['checklist
diff'](https://github.com/jar398/cldiff) tool.

### Reproducibility note

As of 2026 an article on listtools is in preparation.  We will
endeavor to make the supplementary material reproducible from source
material available on the Internet.  See [reproducibility
notes](doc/colmdd.sh).
