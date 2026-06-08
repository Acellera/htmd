# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from typing import TYPE_CHECKING

import numpy as np
from htmd.units import convert as unitconvert
import logging

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from htmd.metricdata import MetricData


def vampScore(
    data: "MetricData",
    lag: float,
    dim: int | None = None,
    units: str = "frames",
    r: int = 2,
    nfolds: int = 10,
    blocksize: int | None = None,
    random_state: int | None = None,
    return_scores: bool = False,
):
    """Compute the cross-validated VAMP-2 score of projected data.

    The VAMP-2 score measures how much of the slow dynamics a featurization
    captures. Evaluated in a cross-validated manner it gives an objective way to
    compare featurizations and to choose the TICA lag time and number of
    dimensions: a higher score is better, the trivial baseline (no slow
    processes resolved) is approximately 1.

    Parameters
    ----------
    data : :class:`MetricData <htmd.metricdata.MetricData>` object
        The projected per-trajectory features to score.
    lag : float
        The VAMP lag time, in the units given by ``units``.
    dim : int, optional
        Number of VAMP dimensions to score on. If None, VAMP uses its default
        (variance cutoff). Vary this to choose how many dimensions to keep.
    units : str, optional
        The units of ``lag`` and ``blocksize``. Can be ``'frames'`` or a time
        unit given as a string.
    r : int, optional
        Which VAMP-r score to compute. ``2`` (default) is the VAMP-2 score.
    nfolds : int, optional
        Number of cross-validation folds.
    blocksize : int, optional
        If None (default), cross-validate over whole trajectories. If given,
        split trajectories into blocks of this many frames (in ``units``) and
        cross-validate over the blocks. Use this when you have few long
        trajectories.
    random_state : int, optional
        Seed for the cross-validation fold assignment, for reproducibility.
    return_scores : bool, optional
        If True, return the raw per-fold score array instead of (mean, std).

    Returns
    -------
    mean, std : float
        The mean and standard deviation of the per-fold scores. Returned unless
        ``return_scores`` is True.
    scores : np.ndarray
        The per-fold scores. Returned only if ``return_scores`` is True.

    Examples
    --------
    >>> from htmd.projections.vamp import vampScore
    >>> mean, std = vampScore(data, lag=20, dim=4, units="ns")
    """
    from deeptime.decomposition import VAMP, vamp_score_cv

    lag = unitconvert(units, "frames", lag, fstep=data.fstep)
    if lag == 0:
        raise ValueError(
            "Lag time conversion resulted in 0 frames. Please use a larger lag time."
        )

    estimator = VAMP(lagtime=lag, dim=dim)
    if blocksize is None:
        if len(data.dat) < 2:
            raise ValueError(
                "Cross-validating over whole trajectories needs at least 2 trajectories, "
                f"but the data has {len(data.dat)}. Provide more trajectories or set "
                "blocksize to split a trajectory into blocks."
            )
        scores = vamp_score_cv(
            estimator,
            data.dat,
            n=nfolds,
            r=r,
            blocksplit=False,
            random_state=random_state,
        )
    else:
        blocksize = unitconvert(units, "frames", blocksize, fstep=data.fstep)
        scores = vamp_score_cv(
            estimator,
            data.dat,
            n=nfolds,
            r=r,
            blocksize=blocksize,
            blocksplit=True,
            random_state=random_state,
        )

    scores = np.asarray(scores)
    if return_scores:
        return scores
    return float(np.mean(scores)), float(np.std(scores))
