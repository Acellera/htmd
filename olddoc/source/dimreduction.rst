Dimensionality Reduction
========================

On top of projection methods it is also highly recommended to use dimensionality reduction methods to further reduce
the space. HTMD provides, for instance, time independent component analysis (TICA) and Kmeans with triangle inequality.
These can be used on MetricData objects. TICA is recommended for Markov Model construction.

Contents:

.. toctree::
    :maxdepth: 1

    TICA - Time independent component analysis <htmd.projections.tica>
    KMeansTri - kmeans triangle inequality <htmd.projections.kmeanstri>
    GWPCA Principal Component Analysis <htmd.projections.gwpca>