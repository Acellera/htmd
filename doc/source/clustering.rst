Clustering
==========

Clustering is done using the scikit-learn clustering library(http://scikit-learn.org/stable/modules/clustering.html#clustering). 
But other clustering classes can be used as well as long as they adhere to the same interface (methods: fit, attributes: cluster_centers_, labels_)

For example MiniBatchKMeans (http://scikit-learn.org/stable/modules/generated/sklearn.cluster.MiniBatchKMeans.html) can be directly passed to the cluster command of MetricData
````
  metricdata.cluster(MiniBatchKMeans(n_clusters=1000), mergesmall=5)
````

Contents:

.. toctree::
   :maxdepth: 2

   MiniBatchKMeans   (http://scikit-learn.org/stable/modules/generated/sklearn.cluster.MiniBatchKMeans.html)
   KCenters clustering method <htmd.clustering.kcenters>
   Regular clustering  <htmd.clustering.regular>
