Clustering
==========

Clustering is done using the `scikit-learn clustering library`_. Other clustering classes can be used as long as they
adhere to the same interface (Methods: fit; Attributes: cluster_centers\_, labels\_).

For example, `MiniBatchKMeans`_ can be directly passed to the cluster command of MetricData::

    metricdata.cluster(MiniBatchKMeans(n_clusters=1000), mergesmall=3)


.. _scikit-learn clustering library: http://scikit-learn.org/stable/modules/clustering.html#clustering
.. _MiniBatchKMeans: http://scikit-learn.org/stable/modules/generated/sklearn.cluster.MiniBatchKMeans.html

Contents:

.. toctree::
    :maxdepth: 1

    KCenters clustering method <htmd.clustering.kcenters>
    RegCluster regular sized clustering  <htmd.clustering.regular>
