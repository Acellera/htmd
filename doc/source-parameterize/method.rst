Method
======

A method implemented in ``parameterize`` is inspired by GAAMP [#]_. A detailed description of our method
will be published soon [#]_.

.. rubric:: General procedure

#. Get initial parameters (from GAFF or CGenFF)
#. Minimize the geometry of molecule with QM
#. Get atomic charges:
    #. Compute ESP with QM
    #. Fit the atomic charges to reproduce the QM results
#. Get dihedral angle parameters:
    #. Detect rotatable dihedral angles
    #. Scan the dihedral angles with QM
    #. Fit the parameters to reproduce the QM results
#. Check parameterization quality:
    #. Plot rotamer QM and MM energies

.. [#]  L. Huang and B. Roux, Automated Force Field Parameterization for Nonpolarizable and Polarizable
        Atomic Models Based on Ab Initio Target Data, J. Chem. Theory Comput., 2013, 9 (8), pp 3543–3556.
        doi: `10.1021/ct4003477 <http://dx.doi.org/10.1021/ct4003477>`_
.. [#]  M. J. Harvey, J. M. Damas, G. Martínez and G. De Fabritiis, Small Molecule Forcefield Parameterization with
        HTMD, in preparation.