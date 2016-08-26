Simulations
==============================

HTMD allows to prepare molecular simulations and run them with minimal knowledge of the technical details 
required to do so, thanks to prepared protocols. Some of these functionalities are general, some are 
specific to the acemd molecular simulation package, however the MD engine should not matter. In particular,
there are functions to setup the configuration of a simulation into a specific directory and way to perform adaptive sampling
methods which is the primary method of sampling in HTMD.

References

*  S.Doerr and G. De Fabritiis, On-the-fly learning and sampling of ligand binding by high-throughput molecular simulations, J. Chem. Theory Comput. 10 (5), pp 2064–2069(2014).

*  M.Harvey, G. Giupponi and G. De Fabritiis, ACEMD: Accelerated molecular dynamics simulations in the microseconds timescale, J. Chem. Theory and Comput. 5, 1632 (2009). 

*  M.J. Harvey and G. De Fabritiis, Acecloud: Molecular Dynamics Simulations in the Cloud, J. Chem. Inf. Model., 2015, 55 (5), pp 909–914.

Contents:

.. toctree::
   :maxdepth: 1

   htmd.acemd
   htmd.adaptive  
   htmd.protocols
