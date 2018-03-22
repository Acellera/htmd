# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

from joblib import Parallel, delayed  # Import delayed as well for other modules

"""
A wrapper for joblib.Parallel to allow custom progress bars.
"""

from tqdm import tqdm

def ParallelExecutor(**joblib_args):
    def aprun(**tq_args):
        tqdm_f = lambda x, args: tqdm(x, **args)
        return lambda x: Parallel(**joblib_args)(tqdm_f(x, tq_args))
    return aprun