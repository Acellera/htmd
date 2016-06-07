# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from __future__ import print_function
#from htmd.projections.metricsecondarystructure import MetricSecondaryStructure
from htmd.home import home
from htmdx.cli import check_registration, show_news
from htmd.latest import compareVersions
from htmd.config import config
import logging.config
import htmd

# -------- Shortcuts ---------
import os

# No longer impoty matplotlib here, as it breaks
# Parameterise's import and attmept to set alternate
# render back-end

#from matplotlib import pylab as plt
# ----------------------------
try:
	logging.config.fileConfig(os.path.join(home(), 'logging.ini'), disable_existing_loggers=False)
except:
	print("HTMD: Logging setup failed")

check_registration(product='htmd')
show_news()
compareVersions()

config()
