# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from moleculekit.projections.metricdihedral import MetricDihedral, Dihedral
from htmd.projections.metric import Metric
from moleculekit.molecule import Molecule
from joblib import Parallel, delayed
from sklearn.metrics import mutual_info_score
import numpy as np
import pandas as pd


class MutualInformation:
    def __init__(self, model, mol=None, fstep=0.1, skip=1):
        """ Class that calculates the mutual information of protein residues.

        Parameters
        ----------
        model : :class:`Model <htmd.model.Model>` object
            A Model object with a calculated MSM
        mol : :class:`Molecule <moleculekit.molecule.Molecule>` object
            A reference molecule from which to obtain structural information. By default model.data.simlist[0].molfile
            will be used.
        fstep : float
            The frame step of the simulations
        skip : int
            Frame skipping

        Examples
        --------
        >>> from htmd.mutualinformation import MutualInformation
        >>> from htmd.ui import *
        >>>
        >>> sims = simlist(glob('./filtered/*/'), './filtered/filtered.pdb')
        >>> mol = Molecule('./filtered/filtered.pdb')
        >>>
        >>> metr = Metric(sims)
        >>> metr.set(MetricDihedral())
        >>> datadih = metr.project()
        >>> datadih.fstep = 0.1
        >>> datadih.dropTraj()
        >>>
        >>> metr = Metric(datadih.simlist)
        >>> metr.set(MetricSelfDistance('protein and name CA', metric='contacts'))
        >>> dataco = metr.project()
        >>> dataco.fstep = 0.1
        >>> dataco.dropTraj()
        >>>
        >>> tica = TICA(datadih, 20, units='ns')
        >>> datatica = tica.project(3)
        >>> datatica.cluster(MiniBatchKMeans(n_clusters=1500))
        >>> model = Model(datatica)
        >>> model.markovModel(12, 4, units='ns')
        >>>
        >>> mu = MutualInformation(model)
        >>> mu.calculate()
        >>> mu.saveMI('./mi_matrix.npy')
        >>> mu.weightGraph(dataco, 0.005)
        >>> mu.save_graphml('./weightgraph.graphml')
        """
        self.model = model

        self.mol = mol
        if mol is None:
            self.mol = Molecule(self.model.data.simlist[0].molfile)

        self._computeChiDihedrals(fstep=fstep, skip=skip)

        self.bindihcat = self._histogram()  # Lasts two minutes
        self.resids = self.mol.get('resid', 'name CA')
        self.residmap = np.ones(self.mol.resid.max() + 1, dtype=int) * -1
        self.residmap[self.resids] = np.arange(len(self.resids))

        self.mi_matrix = None
        self.graph_array = None
        self.graph = None


    def calculate(self, njobs=1):
        """
        Parameters
        ----------
        njobs : int
            Number of parallel jobs to spawn for the calculation of MI
        """
        from htmd.config import _config
        from htmd.parallelprogress import ParallelExecutor
        numchi = self.chi.numDimensions
        statdist = self.model.msm.stationary_distribution
        stconcat = np.concatenate(self.model.data.St)
        microcat = self.model.micro_ofcluster[stconcat]

        aprun = ParallelExecutor(n_jobs=njobs)
        res = aprun(total=numchi, desc='Calculating MI')(delayed(self._parallelAll)(numchi, dih1, 4, self.model.micronum, self.bindihcat, microcat, statdist) for dih1 in range(numchi))
        MI_all = np.zeros((len(self.resids), len(self.resids)))
        for r in res:
            dihcounts = r[0]
            pairs = r[1]
            for dihc, p in zip(dihcounts, pairs):
                dih1, dih2 = p
                if dih1 == dih2:
                    continue
                resid1 = self.residmap[self.mol.resid[self.chi.description.atomIndexes[dih1][0]]]
                resid2 = self.residmap[self.mol.resid[self.chi.description.atomIndexes[dih2][0]]]
                MI_all[resid1][resid2] = self._calcMutualInfo(dihc)
        self.mi_matrix = self._cleanautocorrelations(MI_all)

    def _computeChiDihedrals(self, fstep=0.1, skip=1):
        chis = []
        protmol = self.mol.copy()
        protmol.filter('protein')
        caidx = self.mol.atomselect('protein and name CA')
        resids = self.mol.resid[caidx]
        resnames = self.mol.resname[caidx]
        for residue, resname in zip(resids, resnames):
            ch = Dihedral.chi1(protmol, residue)
            if ch is not None:
                chis.append(ch)

        metr = Metric(self.model.data.simlist, skip=skip)
        metr.set(MetricDihedral(chis, sincos=False))
        data = metr.project()
        data.fstep = fstep
        self.chi = data

    def _histogram(self):
        condata = np.concatenate(self.chi.dat)
        bins = np.array([-180, -150, -90, -30, 0, 30, 90, 150, 180])
        dic = {1: 3, 2: 0, 3: 1, 4: 0, 5: 0, 6: 2, 7: 0, 8: 3, 9:3}
        binneddih = np.zeros([condata.shape[0], condata.shape[1]])
        for dihedral in range(condata.shape[1]):
            binning = np.digitize(condata[:, dihedral], bins)
            binneddih[:, dihedral] = [dic[n] if n in dic.keys() else n for n in binning]
        return binneddih

    def _calcMutualInfo(self, contingency):
        # Ripped out of sklearn since it converts floats to integers without warning which breaks our use-case
        from math import log
        nzx, nzy = np.nonzero(contingency)
        nz_val = contingency[nzx, nzy]
        contingency_sum = contingency.sum()
        pi = np.ravel(contingency.sum(axis=1))
        pj = np.ravel(contingency.sum(axis=0))
        log_contingency_nm = np.log(nz_val)
        contingency_nm = nz_val / contingency_sum
        # Don't need to calculate the full outer product, just for non-zeroes
        outer = pi.take(nzx) * pj.take(nzy)
        log_outer = -np.log(outer) + log(pi.sum()) + log(pj.sum())
        mi = (contingency_nm * (log_contingency_nm - log(contingency_sum)) +
              contingency_nm * log_outer)
        return mi.sum()

    def _parallelAll(self, numdih, dih1, numbins, micronum, bindihcat, microcat, stat_dist):
        results = []
        resultpairs = []
        for dih2 in range(dih1, numdih):
            dihcounts = np.zeros((numbins, numbins, micronum))

            # Find pairs of dihedrals (keys) and which absolute frames they occur in (vals)
            df = pd.DataFrame({'a': list(zip(bindihcat[:, dih1], bindihcat[:, dih2]))})
            gg = df.groupby(by=df.a).groups

            for pair in gg:
                microsofpairs = microcat[gg[pair]]  # Get the microstates of all frames having given dihedral pair
                microsofpairs = np.delete(microsofpairs, np.where(microsofpairs == -1)[0])  # Delete dropped clusters
                counts = np.bincount(microsofpairs)  # Count number of frames with that pair for each microstate (0, max_micro_seen)
                dihcounts[int(pair[0]), int(pair[1]), :len(counts)] += counts  # Add the counts
            dihcounts /= dihcounts.sum(axis=0).sum(axis=0)  # Normalize all slices by total counts in each microstate
            dihcounts *= stat_dist  # Multiply by stationary distribution of each state
            dihcounts = np.sum(dihcounts, axis=2)  # Calculate the weighted sum over all states
            results.append(dihcounts)
            resultpairs.append((dih1, dih2))

        return results, resultpairs

    def _cleanautocorrelations(self, mi):
        # TODO: vectorize this
        # np.diag(np.ones(3), k=-1).astype(bool) | np.diag(np.ones(3), k=1).astype(bool) | np.diag(np.ones(4)).astype(bool)
        for i in range(mi.shape[0]):
            for j in range(mi.shape[1]):
                if abs(i - j) < 2:
                    mi[i][j] = 0
        return mi

    def saveMI(self, path):
        np.save(path, self.mi_matrix)

    def loadMI(self, path):
        self.mi_matrix = np.load(path)

    def weightGraph(self, datacontacts, mi_threshold, time_treshold=0.6):
        if len(self.mol.get('resid', 'name CA')) != len(self.resids):
            raise Exception('The length of the protein doesn\'t match the Mutual Information data')
        contactcat = np.concatenate(datacontacts.dat)
        contacts_matrix = np.zeros([len(self.resids), len(self.resids)])
        for i in range(contactcat.shape[1]):
            counter = np.count_nonzero(contactcat[:, i])
            resid1 = self.residmap[self.mol.resid[datacontacts.description.atomIndexes[i][0]]]
            resid2 = self.residmap[self.mol.resid[datacontacts.description.atomIndexes[i][1]]]
            contacts_matrix[resid1][resid2] = counter

        self.graph_array = np.zeros([contacts_matrix.shape[0], contacts_matrix.shape[0]])
        mask = (self.mi_matrix > mi_threshold) & (contacts_matrix > (time_treshold * contactcat.shape[0]))
        self.graph_array[mask] = self.mi_matrix[mask]

        intermed = []
        for source in range(self.graph_array.shape[0]):
            for target in range(source, self.graph_array.shape[1]):
                if self.graph_array[source, target] != 0 and target > source:
                    intermed.append(
                        [int(self.resids[source]), int(self.resids[target]), float(self.graph_array[source, target])])
        import pandas as pd
        import networkx as nx
        from sklearn.cluster.spectral import SpectralClustering

        pd = pd.DataFrame(intermed, columns=['source', 'target', 'weight'])
        pd[['source', 'target']] = pd[['source', 'target']].astype(type('int', (int,), {}))
        pd['weight'] = pd['weight'].astype(type('float', (float,), {}))
        G = nx.from_pandas_edgelist(pd, 'source', 'target', 'weight')
        ## setSegment
        segids = self.mol.get('segid', 'name CA')
        seg_res_dict = {key: value for (key, value) in zip(self.resids, segids) if
                        np.any(pd.loc[(pd['source'] == key)].index) or np.any(pd.loc[(pd['target'] == key)].index)}
        nx.set_node_attributes(G,  seg_res_dict, 'Segment')
        ## set
        if not nx.is_connected(G):
            G = max(nx.connected_component_subgraphs(G), key=len)
        flow_cent = nx.current_flow_betweenness_centrality(G, weight='weight')
        nx.set_node_attributes(G, flow_cent, 'flowcent')
        Spectre = SpectralClustering(n_clusters=10, affinity='precomputed')
        model = Spectre.fit_predict(self.graph_array)
        model = model.astype(type('float', (float,), {}))
        spectral_dict = {key: value for (key, value) in zip(self.resids, model) if key in G.nodes()}
        nx.set_node_attributes(G, spectral_dict, 'spectral')
        self.graph = G

    def save_graphml(self, path):
        import networkx as nx
        nx.write_graphml(self.graph, path)

    # def save_pdf(self, graphpath, outpath):
    #     print(self.graph_array)
    #     self.save_graphml(graphpath)
    #     compile_java('./GephiGraph.java')
    #     compile_java('./GraphTest.java')
    #     execute_java('GraphTest', graphpath, outpath)

# if __name__ == '__main__':
#     from htmd.mutualinformation import MutualInformation
#     from htmd import *
#
#     sims = simlist(glob('/shared/adria/2ov5/adaptive_amber/batches/1/filtered/*/'),
#                    '/shared/adria/2ov5/adaptive_amber/batches/1/filtered/filtered.pdb')
#     mol = Molecule('/shared/adria/2ov5/adaptive_amber/batches/1/filtered/filtered.pdb')
#
#     metr = Metric(sims[0:100])
#     metr.set(MetricDihedral())
#     datadih = metr.project()
#     datadih.fstep = 0.1
#     datadih.dropTraj()
#
#     metr = Metric(datadih.simlist)
#     metr.set(MetricSelfDistance('protein and name CA', metric='contacts', threshold=8))
#     dataco = metr.project()
#     dataco.fstep = 0.1
#     dataco.dropTraj()
#
#     tica = TICA(datadih, 20, units='ns')
#     datatica = tica.project(3)
#     datatica.cluster(MiniBatchKMeans(n_clusters=1500))
#     model = Model(datatica)
#     model.markovModel(12, 4, units='ns')
#
#     mu = MutualInformation(model)
#     mu.calculate()
#     mu.saveMI('/tmp/mi_matrix.npy')
#     mu.weightGraph(dataco, 0.005)
#     mu.save_graphml('/tmp/weightgraph_0-005.graphml')
