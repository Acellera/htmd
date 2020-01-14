MD simulations
==============

HTMD allows to prepare molecular simulations and run them with minimal knowledge of the technical details 
required to do so, thanks to prepared protocols. 

All functionalities of HTMD are general, it is possible to use HTMD with any MD engine. However, ACEMD the MD engine
embedded in HTMD is easier to use and is supported.   In particular, there are functions to setup the configuration of
a simulation into a specific directory and way to perform adaptive sampling methods which is the primary method of
sampling in HTMD.

ACEMD, a powerful and simple MD engine which has pioneered GPU computing since 2009, is distributed together with HTMD
in a standard version. A professional version is available by writing at `info@acellera.com`_ (for more information
visit `www.acellera.com/acemd`_).

.. _info@acellera.com: info@acellera.com
.. _www.acellera.com/acemd: http://www.acellera.com/acemd

References

*  S. Doerr and G. De Fabritiis, On-the-fly learning and sampling of ligand binding by high-throughput molecular
   simulations, J. Chem. Theory Comput., 2014, 10 (5), pp 2064–2069. doi: `10.1021/ct400919u`_

*  M. J. Harvey, G. Giupponi and G. De Fabritiis, ACEMD: Accelerated molecular dynamics simulations in the microseconds
   timescale, J. Chem. Theory and Comput., 2009, 5 (6), 1632. doi: `10.1021/ct9000685`_

*  M. J. Harvey and G. De Fabritiis, AceCloud: Molecular Dynamics Simulations in the Cloud, J. Chem. Inf. Model.,
   2015, 55 (5), pp 909–914. doi: `10.1021/acs.jcim.5b00086`_

.. _10.1021/ct400919u: http://dx.doi.org/10.1021/ct400919u
.. _10.1021/ct9000685: http://dx.doi.org/10.1021/ct9000685
.. _10.1021/acs.jcim.5b00086: http://dx.doi.org/10.1021/acs.jcim.5b00086

Contents:

.. toctree::
    :maxdepth: 1

    Applications <htmd.apps>
    Queues <../jobqueues/jobqueues.rst>
    Adaptive sampling <adaptive>
    State-of-the-art protocols for molecular simulations <htmd.protocols>
