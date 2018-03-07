# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from __future__ import print_function
import htmd.home
from htmd.session import htmdsave, htmdload
from htmd.simlist import simlist, simfilter, simmerge
from htmd.metricdata import MetricData
from htmd.projections.metricdistance import MetricDistance, MetricSelfDistance, reconstructContactMap
from htmd.projections.metricrmsd import MetricRmsd
from htmd.projections.metricfluctuation import MetricFluctuation
from htmd.projections.metriccoordinate import MetricCoordinate
from htmd.projections.metricdihedral import MetricDihedral, Dihedral
from htmd.projections.metricshell import MetricShell
from htmd.projections.metricsecondarystructure import MetricSecondaryStructure
from htmd.projections.metricsasa import MetricSasa
from htmd.projections.metrictmscore import MetricTMscore
from htmd.projections.metric import Metric
from htmd.projections.tica import TICA
from htmd.projections.kmeanstri import KMeansTri
from htmd.projections.gwpca import GWPCA
from htmd.molecule.molecule import Molecule
from htmd.adaptive.adaptiverun import AdaptiveMD
from htmd.adaptive.adaptivegoal import AdaptiveGoal
from htmd.adaptive.adaptive import reconstructAdaptiveTraj
from htmd.model import Model, getStateStatistic
from htmd.kinetics import Kinetics
from htmd.vmdviewer import viewer, getCurrentViewer
from htmd.builder.solvate import solvate
from htmd.apps.acemd_v2 import Acemd
from htmd.builder.builder import detectDisulfideBonds, autoSegment, embed, DisulfideBridge
import htmd.builder.charmm as charmm
import htmd.builder.amber as amber
from htmd.molecule.util import uniformRandomRotation
from htmd.rotationmatrix import rotationMatrix
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
from htmdx.cli import check_registration, show_news
from htmd.latest import compareVersions
from htmd.config import config
import logging.config

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
try:
    logging.config.fileConfig(os.path.join(htmd.home.home(), 'logging.ini'), disable_existing_loggers=False)
except:
    print("HTMD: Logging setup failed")



