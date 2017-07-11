Adaptive sampling
=================

HTMD is build around adaptive sampling, i.e. use on-the-fly information from the current data to decide where to restart new simulations. This in practice works really well with speed-up of over an order of magnitude compared to standard sampling methods. Adaptive sampling replaces to a large extent any need for biased sampling in a much safer way than biasing. Also it does not require to have very good reaction coordinates between it is building its own reaction space on-the-fly while sampling. It also integrate very well with Markov state models.

Relevant papers to read:

*  S.Doerr and G. De Fabritiis, On-the-fly learning and sampling of ligand binding by high-throughput molecular simulations, J. Chem. Theory Comput. 10 (5), pp 2064–2069(2014).

*  S.Doerr , M.J. Harvey, F. Noé ,G. De Fabritiis, HTMD: High-throughput molecular dynamics for molecular discovery, J. Chem. Theory Comput., 2016, 12 (4), pp 1845–1852


Contents:

.. toctree::
   :maxdepth: 2

   Adaptive sampling <htmd.adaptive.adaptiverun>
