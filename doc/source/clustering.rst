Clustering
==========

Clustering is done using the `scikit-learn clustering library <http://scikit-learn.org/stable/modules/clustering.html#clustering>`_ . 
But other clustering classes can be used as well as long as they adhere to the same interface (methods: fit, attributes: cluster_centers, labels).

For example `MiniBatchKMeans <http://scikit-learn.org/stable/modules/generated/sklearn.cluster.MiniBatchKMeans.html>`_ can be directly passed to the cluster command of MetricData::

	metricdata.cluster(MiniBatchKMeans(n_clusters=1000), mergesmall=3)


Contents:

.. toctree::
   :maxdepth: 2

   KCenters clustering method <htmd.clustering.kcenters>
   RegCluster regular sized clustering  <htmd.clustering.regular>
