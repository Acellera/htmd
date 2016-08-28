Molecular dynamics simulations
==============================

HTMD allows to prepare molecular simulations and run them with minimal knowledge of the technical details 
required to do so, thanks to prepared protocols. 

All functionalities of HTMD are general, it is possible to use HTMD with any MD engine. However, ACEMD the MD engine embedded in HTMD is easier to use and is supported.   In particular, there are functions to setup the configuration of a simulation into a specific directory and way to perform adaptive sampling methods which is the primary method of sampling in HTMD.

ACEMD a powerful and simple MD engine which pioneered  GPU computing since 2009 is distributed together with HTMD in a standard version. A profesional version is available by writing at info@acellera.com (http://www.acellera.com/acemd).

Nevertheless, 

References

*  S.Doerr and G. De Fabritiis, On-the-fly learning and sampling of ligand binding by high-throughput molecular simulations, J. Chem. Theory Comput. 10 (5), pp 2064–2069(2014).

*  M.Harvey, G. Giupponi and G. De Fabritiis, ACEMD: Accelerated molecular dynamics simulations in the microseconds timescale, J. Chem. Theory and Comput. 5, 1632 (2009). 

*  M.J. Harvey and G. De Fabritiis, Acecloud: Molecular Dynamics Simulations in the Cloud, J. Chem. Inf. Model., 2015, 55 (5), pp 909–914.

Contents:

.. toctree::
   :maxdepth: 1

   The ACEMD MD engine <htmd.acemd>
   Adaptive sampling <htmd.adaptive> 
   State-of-the-art protocols for molecular simulations <htmd.protocols>
