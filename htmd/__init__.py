# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from __future__ import print_function
from htmd.home import home
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
from htmd.userinterface import UserInterface
from htmd.molecule.molecule import Molecule
from htmd.adaptive.adaptiverun import AdaptiveMD
from htmd.adaptive.adaptivegoal import AdaptiveGoal
from htmd.adaptive.adaptive import reconstructAdaptiveTraj
from htmd.model import Model, getStateStatistic
from htmd.kinetics import Kinetics
from htmd.vmdviewer import viewer, getCurrentViewer
from htmd.builder.solvate import solvate
from htmd.apps.acemd import Acemd
from htmd.apps.pmemd import Pmemd
from htmd.apps.acemdlocal import AcemdLocal
from htmd.apps.pmemdlocal import PmemdLocal
from htmd.apps.lsf import LSF
from htmd.apps.aws import AWS
from htmd.builder.builder import detectDisulfideBonds, segmentgaps, autoSegment, embed, DisulfideBridge
import htmd.builder.charmm as charmm
import htmd.builder.amber as amber
from htmd.molecule.util import uniformRandomRotation
from htmd.rotationmatrix import rotationMatrix
from htmd.builder.preparation import proteinPrepare
from htmd.dock import dock
from htmdx.cli import check_registration, show_news
from htmd.latest import compareVersions
from htmd.util import tempname
from htmd.util import testDHFR
from htmd.config import config
from htmd.clustering.kcenters import KCenter
from htmd.clustering.regular import RegCluster
from htmd.queues.localqueue import LocalGPUQueue
from htmd.queues.slurmqueue import SlurmQueue
from htmd.queues.lsfqueue import LsfQueue
from htmd.queues.pbsqueue import PBSQueue
from htmd.vmdgraphics import VMDConvexHull, VMDBox, VMDIsosurface, VMDSphere, VMDText
from htmd.builder.loopmodeler import loopModeller
import logging.config

from htmd.version import version as _version
__version__ = _version()

import htmd

# -------- Shortcuts ---------
import os
import numpy as np
import math
import shutil
from glob import glob
from sklearn.cluster import MiniBatchKMeans

# No longer impoty matplotlib here, as it breaks
# Parameterise's import and attmept to set alternate
# render back-end

# from matplotlib import pylab as plt
# ----------------------------
try:
    logging.config.fileConfig(os.path.join(home(), 'logging.ini'), disable_existing_loggers=False)
except:
    print("HTMD: Logging setup failed")

if not (os.getenv("HTMD_NONINTERACTIVE")):
    check_registration(product='htmd')
    show_news()
    compareVersions()

import progress_reporter.bar.gui as __gui  # Disabling pyemma progress widgets
__gui.ipython_notebook_session = False

config()

if os.getenv('HTMD_CONFIG'):
    configfile = os.getenv('HTMD_CONFIG')
    if not os.path.isfile(configfile):
        raise FileNotFoundError('HTMD Config file {} does not exist (Check HTMD_CONFIG).'.format(configfile))
    elif not configfile.endswith('.py'):
        raise Warning('HTMD Config file {} may not be a python file (Check HTMD_CONFIG).'.format(configfile))
    else:
        try:
            with open(configfile) as f:
                code = compile(f.read(), configfile, 'exec')
                exec(code)
        except:
            raise RuntimeError('Failed to execute the HTMD Config file {}.'.format(configfile))
        else:
            print('\nHTMD Config file {} executed.'.format(configfile))
