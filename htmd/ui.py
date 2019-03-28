# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from __future__ import print_function
import htmd.home
from htmd.session import htmdsave, htmdload
from htmd.simlist import simlist, simfilter, simmerge
from htmd.metricdata import MetricData
from moleculekit.projections.metricdistance import MetricDistance, MetricSelfDistance, reconstructContactMap
from moleculekit.projections.metricrmsd import MetricRmsd
from moleculekit.projections.metricfluctuation import MetricFluctuation
from moleculekit.projections.metriccoordinate import MetricCoordinate
from moleculekit.projections.metricdihedral import MetricDihedral, Dihedral
from moleculekit.projections.metricshell import MetricShell
from moleculekit.projections.metricsecondarystructure import MetricSecondaryStructure
from moleculekit.projections.metricsasa import MetricSasa
from moleculekit.projections.metrictmscore import MetricTMscore
from htmd.projections.metric import Metric
from htmd.projections.tica import TICA
from htmd.projections.kmeanstri import KMeansTri
from htmd.projections.gwpca import GWPCA
from moleculekit.molecule import Molecule
from htmd.adaptive.adaptiverun import AdaptiveMD
from htmd.adaptive.adaptivegoal import AdaptiveGoal
from htmd.adaptive.adaptive import reconstructAdaptiveTraj
from htmd.model import Model, getStateStatistic
from htmd.kinetics import Kinetics
from htmd.vmdviewer import viewer, getCurrentViewer
from htmd.builder.solvate import solvate
from htmd.mdengine.acemd.acemd import Acemd, Acemd2, AtomRestraint, GroupRestraint
from htmd.builder.builder import detectDisulfideBonds, autoSegment, embed, DisulfideBridge
import htmd.builder.charmm as charmm
import htmd.builder.amber as amber
from moleculekit.util import uniformRandomRotation
from moleculekit.util import rotationMatrix
from htmd.builder.preparation import proteinPrepare
from htmd.dock import dock
from htmd.util import tempname
from htmd.util import testDHFR
from htmd.clustering.kcenters import KCenter
from htmd.clustering.regular import RegCluster
from htmd.queues.localqueue import LocalGPUQueue, LocalCPUQueue
from htmd.queues.slurmqueue import SlurmQueue
from htmd.queues.lsfqueue import LsfQueue
from htmd.queues.pbsqueue import PBSQueue
from htmd.vmdgraphics import VMDConvexHull, VMDBox, VMDIsosurface, VMDSphere, VMDText
from htmd.builder.loopmodeler import loopModeller
from htmd.ffevaluation.ffevaluate import FFEvaluate
from htmdx.cli import check_registration, show_news
from htmd.latest import compareVersions
from htmd.config import config

# -------- Shortcuts ---------
import os
import numpy as np
import math
import shutil
from glob import glob
from sklearn.cluster import MiniBatchKMeans

if not (os.getenv("HTMD_NONINTERACTIVE")):
    check_registration(product='htmd')
    show_news()
    compareVersions()

# No longer import matplotlib here, as it breaks
# Parameterise's import and attmept to set alternate
# render back-end

# from matplotlib import pylab as plt
# ----------------------------
import warnings

# Get rid of pyemma version warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', category=UserWarning)
    import pyemma



