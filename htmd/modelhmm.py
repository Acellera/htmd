# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
import logging
logger = logging.getLogger(__name__)


class ModelHMM(object):
    def __init__(self, data=None):
        if data is None:
            return
        self.data = data

    def markovModel(self, lag, macronum):
        from pyemma.msm import estimate_hidden_markov_model
        self.hmm = estimate_hidden_markov_model(self.data.St.tolist(), macronum, lag, connectivity='largest')
        print('Active set includes macrostates: {}'.format(self.hmm.active_set))

    @property
    def macronum(self):
        return len(self.hmm.active_set)

    def plotTimescales(self, maxlag, numstates):
        import pyemma.msm as msm
        import pyemma.plots as mplt
        its = msm.timescales_hmsm(self.data.St.tolist(), numstates, lags=maxlag)
        mplt.plot_implied_timescales(its, ylog=True, units='ns', dt=self.data.fstep, linewidth=2)

    def eqDistribution(self, plot=True):
        if plot:
            from matplotlib import pylab as plt
            plt.bar(range(self.macronum), self.hmm.stationary_distribution)
            #ax = plt.gca()
            #ax.set_xticklabels([str(a) for a in self.hmm.active_set])
            plt.xticks(np.arange(self.macronum)+0.4, [str(a) for a in self.hmm.active_set])
        return self.hmm.stationary_distribution

    def viewStates(self, protein=None, ligand=None, nsamples=20):
        from htmd.projections.metric import _singleMolfile
        from moleculekit.molecule import Molecule
        from moleculekit.vmdviewer import getCurrentViewer
        (single, molfile) = _singleMolfile(self.data.simlist)
        if not single:
            raise RuntimeError('Can''t visualize states without unique molfile')

        viewer = getCurrentViewer()
        colors = [0, 1, 3, 4, 5, 6, 7, 9]

        print('Active set includes macrostates: {}'.format(self.hmm.active_set))

        # dtraj = np.vstack(self.hmm.discrete_trajectories_full)
        res = self.hmm.sample_by_observation_probabilities(nsamples)
        refmol = Molecule(molfile)

        for i, s in enumerate(self.hmm.active_set):
            mol = Molecule(molfile)
            mol.coords = []
            mol.box = []
            # idx = np.where(dtraj == i)[0]
            # samples = np.random.choice(idx, 20)
            # frames = self.data.abs2sim(samples)

            frames = self.data.rel2sim(res[i])
            for f in frames:
                mol._readTraj(f.sim.trajectory[f.piece], frames=[f.frame], append=True)
            mol.wrap('protein')
            mol.align('protein', refmol=refmol)
            viewer.loadMol(mol, name='hmm macro ' + str(s))
            if ligand is not None:
                viewer.rep('ligand', sel=ligand, color=colors[np.mod(i, len(colors))])
            if protein is not None:
                viewer.rep('protein')
            viewer.send('start_sscache')
