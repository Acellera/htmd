Protein preparation with HTMD
===============

Content moved in a tutorial.


Technical notes
---------------

The version imported is based on a modification of the
[apbs-pdb2pqr](https://github.com/Electrostatics/apbs-pdb2pqr)
code hosted on GitHub.  Part of the conversion involved those steps:

 * Removed autoconf-related files
 * Apply `2to3`
 * Amend `string.xx()` -> `str.xx()` and several others
 * Ported to propka31

and is available at this [forked repository](https://github.com/tonigi/apbs-pdb2pqr).
