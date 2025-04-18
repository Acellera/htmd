# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import warnings
import htmd.home
from htmd.session import htmdsave, htmdload
from htmd.simlist import simlist, simfilter, simmerge, Sim
from htmd.metricdata import MetricData
from moleculekit.projections.metricdistance import (
    MetricDistance,
    MetricSelfDistance,
    reconstructContactMap,
)
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
from moleculekit.vmdviewer import viewer, getCurrentViewer
from htmd.builder.solvate import solvate
from htmd.builder.builder import detectDisulfideBonds, embed
import htmd.builder.charmm as charmm
import htmd.builder.amber as amber
from moleculekit.util import uniformRandomRotation
from moleculekit.util import rotationMatrix
from htmd.dock import dock
from htmd.util import tempname
from htmd.util import testDHFR
from htmd.clustering.kcenters import KCenter
from htmd.clustering.regular import RegCluster
from jobqueues.localqueue import LocalGPUQueue, LocalCPUQueue
from jobqueues.slurmqueue import SlurmQueue
from jobqueues.lsfqueue import LsfQueue
from jobqueues.pbsqueue import PBSQueue
from moleculekit.vmdgraphics import (
    VMDConvexHull,
    VMDBox,
    VMDIsosurface,
    VMDSphere,
    VMDText,
)
from moleculekit.tools.autosegment import autoSegment
from htmd.builder.loopmodeler import loopModeller

try:
    from moleculekit.tools.preparation import systemPrepare, proteinPrepare
except Exception as e:
    warnings.warn(f"{e}")


from htmd.license import htmd_registration, htmd_show_news, htmd_register
from htmd.latest import compareVersions
from htmd.config import config

# -------- Shortcuts ---------
import os
import numpy as np
import math
import shutil
from glob import glob
from sklearn.cluster import MiniBatchKMeans

if (
    not (os.getenv("HTMD_NONINTERACTIVE"))
    and (not os.getenv("CI"))
    and (not os.getenv("APPTAINER_CONTAINER"))
):
    htmd_show_news()
    compareVersions()
    htmd_registration()
