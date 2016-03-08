# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from __future__ import print_function
#from htmd.projections.metricsecondarystructure import MetricSecondaryStructure
from htmd.home import home
from htmd.session import htmdsave, htmdload
from htmd.simlist import simlist, simfilter, simmerge
from htmd.metricdata import MetricData
from htmd.projections.metricdistance import MetricDistance, MetricSelfDistance
from htmd.projections.metricrmsd import MetricRmsd
from htmd.projections.metriccoordinate import MetricCoordinate
from htmd.projections.metricdiherdal import MetricDihedral
from htmd.projections.metricshell import MetricShell
from htmd.projections.metricsecondarystructure import MetricSecondaryStructure
from htmd.projections.metric import Metric
from htmd.projections.tica import TICA
from htmd.projections.kmeanstri import KMeansTri
from htmd.userinterface import UserInterface
from htmd.molecule.molecule import Molecule
from htmd.adaptive.adaptiverun import AdaptiveRun
from htmd.adaptive.adaptive import reconstructAdaptiveTraj
from htmd.model import Model, getStateStatistic
from htmd.kinetics import Kinetics
from htmd.vmdviewer import viewer, getCurrentViewer
from htmd.builder.solvate import solvate
from htmd.acemd.acemd import Acemd
from htmd.apps.acemdlocal import AcemdLocal
from htmd.apps.lsf import LSF
from htmd.apps.aws import AWS
from htmd.builder.builder import detectDisulfideBonds, segmentgaps, embed, DisulfideBridge
import htmd.builder.charmm as charmm
import htmd.builder.amber as amber
from htmd.molecule.util import uniformRandomRotation
from htmd.dock import dock
from htmdx.cli import check_registration, show_news
from htmd.latest import compareVersions
from htmd.util import tempname
from htmd.config import config
import logging.config
import htmd

# -------- Shortcuts ---------
import os
import numpy as np
import math
import shutil
import tempfile
import random
from glob import glob
from sklearn.cluster import MiniBatchKMeans
from matplotlib import pylab as plt
# ----------------------------

logging.config.fileConfig(os.path.join(home(), 'logging.ini'), disable_existing_loggers=False)

check_registration(product='htmd')
show_news()
compareVersions()

config()
