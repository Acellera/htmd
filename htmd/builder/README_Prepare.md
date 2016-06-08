Protein preparation with HTMD
===============

*Toni Giorgino*


The system preparation phase is based on the PDB2PQR software. It 
includes the following steps (from the
[PDB2PQR algorithm
description](http://www.poissonboltzmann.org/docs/pdb2pqr-algorithm-description/)):

 * Assign titration states at the user-chosen pH;
 * Flipping the side chains of HIS (including user defined HIS states), ASN, and GLN residues;
 * Rotating the sidechain hydrogen on SER, THR, TYR, and CYS (if available);
 * Determining the best placement for the sidechain hydrogen on neutral HIS, protonated GLU, and protonated ASP;
 * Optimizing all water hydrogens.

The hydrogen bonding network calculations are performed by the
[PDB2PQR](http://www.poissonboltzmann.org/) software package. The pKa
calculations are performed by the [PROPKA
3.1](https://github.com/jensengroup/propka-3.1) software packages.
Please see the copyright, license  and citation terms distributed with each.

Note that this version was modified in order to use an 
externally-supplied PROPKA **3.1**, whereas
the original had propKa 3.0 *embedded*!

The results of the function should be roughly equivalent of the system
preparation wizard's preprocessing and optimization steps
of Schrodinger's Maestro software.




Usage
----------

See the docstring for options. You need propka31 installed via conda.
    
    tryp = Molecule('sysprep/tests/3ptb.pdb')
    tryp_op = proteinPrepare(tryp)
    tryp_op.write('sysprep/tests/systempreparation-test-main-ph-7.pdb')



Modified residue names
----------------------

The molecule produced has residue names modified according to their protonation.
Later system-building functions assume these residue names.


Charge +1    |  Neutral   | Charge -1
-------------|------------|----------
 -           |  ASH       | ASP
 -           |  CYS       | CYM
 -           |  GLH       | GLU
HIP          |  HID/HIE   |  -
LYS          |  LYN       |  -
 -           |  TYR       | TYM
ARG          |  AR0       |  -


Forcefields vary in their coverage of alternative charge states.




Citations
---------

Please acknowledge your use of PDB2PQR by citing:

 *   Dolinsky TJ, Czodrowski P, Li H, Nielsen JE, Jensen JH, Klebe G, Baker NA. PDB2PQR: Expanding and upgrading automated preparation of biomolecular structures for molecular simulations. Nucleic Acids Res, 35, W522-5, 2007. 
 *   Dolinsky TJ, Nielsen JE, McCammon JA, Baker NA. PDB2PQR: an automated pipeline for the setup, execution, and analysis of Poisson-Boltzmann electrostatics calculations. Nucleic Acids Res, 32, W665-W667, 2004.
 
 
Please acknowledge your use of PROPKA by citing:

 *   Sondergaard, Chresten R., Mats HM Olsson, Michal Rostkowski, and Jan H. Jensen. "Improved Treatment of Ligands and Coupling Effects in Empirical Calculation and Rationalization of pKa Values." Journal of Chemical Theory and Computation 7, no. 7 (2011): 2284-2295.
 *   Olsson, Mats HM, Chresten R. Sondergaard, Michal Rostkowski, and Jan H. Jensen. "PROPKA3: consistent treatment of internal and surface residues in empirical pKa predictions." Journal of Chemical Theory and Computation 7, no. 2 (2011): 525-537.





Developer notes
===============



The version imported is based on a modification of the
[apbs-pdb2pqr](https://github.com/Electrostatics/apbs-pdb2pqr)
code hosted on GitHub.

Part of the conversion involved those steps:

 * Removed autoconf-related files
 * Apply `2to3`
 * Amend `string.xx()` -> `str.xx()` and several others
 * Ported to propka31

and is available at this [forked repository](https://github.com/tonigi/apbs-pdb2pqr).


