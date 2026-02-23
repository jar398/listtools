## Installation

The tools require no particular installation step.  They can be
invoked from the shell out of a clone of this `git` repository.

Python 3 is assumed.

Install the python3 `regex` module using `pip3` or equivalent.

Also required is `gnparser`, which has download instructions
[here](https://github.com/gnames/gnparser#installation).  I have
tested listtools against version 1.11.8 (July 2025).

The `gnparser` installation instructions are sensitive to local setup
(OS version and homebrew vs. macports vs. rebuild from sources), so
check its repository if you encounter problems
(e.g. [this](https://github.com/gnames/gnparser/issues/274) and ask
the maintainer for help if you get stuck.

If you are working with MDD (the Mammal Diversity Database), you may
want to get version 2 or beyond of the [MDD-DwC-mapping
tool](https://github.com/jar398/MDD-DwC-mapping/blob/main/src/explore_data.py),
a single python file located as of February 2026 in the
'MDD-DwC-mapping' repository, user jar398, at github.com.  You will
need this or something like it to convert MDD files to Darwin Core.
