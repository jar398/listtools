# listtools

Tools for manipulating lists of things, and in particular, lists of taxonomic
names (i.e. 'checklists'), which are sometimes organized hierachically.

* [Installation](doc/guide.md#installation)
* [User guide](doc/guide.md)
* [Method](doc/method.md)

<strong>Warning:</strong> The code is in transition from one method to
another.  You will see redundancy in functionality as the older code
has just been left in place rather than removed or integrated with the
new code.

As is customary in computer science, trees grow downwards, like a root
system, not upwards.  So:

 * The root is at the top of the tree; we go up
   the tree to get to the root.
 * As sets of individual specimens/observations, taxa are therefore
   bigger toward the top/root, and smaller toward the bottom/tips.
 * As an upper semilattice, the root is the top node, and
   is greater than all others; the tips are smallest.  A < B means that
   A is descended from B.

As is customary in biology, we say 'tip' rather than 'leaf' for
generality, in case someone is tempted to take 'leaf' too literally.
Oddly, the same is not done for 'root'.

The genesis for this software was a conversation around 2018 or 2019
with Nico Franz, who suggested that a generic tool for comparing
checklists might be widely useful.  Before this conversation I had been working on
[smasher](https://github.com/opentreeoflife/reference-taxonomy/) (see
[article](https://doi.org/10.3897/BDJ.5.e12581)), which did something
similar, and I felt a desire to rewrite it on a firmer theoretical
foundation to make it "more sound", and to make create a free-standing
tool independent of the Open Tree of Life project or any other framework.

This tool suite replaces the older ['checklist
diff'](https://github.com/jar398/cldiff) tool.

