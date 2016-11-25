# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.progress.progress import ProgressBar
from joblib import Parallel, delayed  # Import delayed as well for other modules

"""
A wrapper for joblib.Parallel to allow custom progress bars.
"""


def progressbar(seq, description=None, total=None):
    p = ProgressBar(total, description=description)
    while True:
        try:
            yield next(seq)
            p.progress() # Had to put progress after yield because last call goes over the total and then I can't decrement in stop()
        except StopIteration:
            p.stop()
            raise


def ParallelExecutor(**joblib_args):
    def aprun(**tq_args):
        def tmp(op_iter):
            foo = lambda args: lambda x: progressbar(x, **args)
            bar_func = foo(tq_args)
            return Parallel(**joblib_args)(bar_func(op_iter))
        return tmp
    return aprun