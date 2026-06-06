# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from typing import TYPE_CHECKING
import numpy as np
import os
import re

if TYPE_CHECKING:
    from htmd.metricdata import MetricData

parentregex = re.compile(r"e\d+s\d+_(e\d+s\d+)p(\d+)f(\d+)")
epochregex = re.compile(r"^e(\d+)s")


def getEpochTrajectoryDictionary(simlist: list) -> dict:
    """Build a mapping from epoch number to simulation indexes.

    Parameters
    ----------
    simlist : list
        A simulation list created using the :func:`simlist <htmd.simlist.simlist>` function.

    Returns
    -------
    simepochs : dict
        A dictionary mapping epoch number (int) to a list of simulation indexes
        from ``simlist`` belonging to that epoch.
    """
    from collections import defaultdict

    simepochs = defaultdict(list)
    for i, sim in enumerate(simlist):
        name = os.path.basename(os.path.dirname(os.path.abspath(sim.trajectory[0])))
        e = int(epochregex.findall(name)[0])
        simepochs[e].append(i)

    return simepochs


def getEpochSimIdx(data: "MetricData", epoch: int) -> np.ndarray:
    """Return the simulation indexes in a MetricData object for a given epoch.

    Parameters
    ----------
    data : MetricData
        Projected simulation data whose ``simlist`` is searched.
    epoch : int
        The epoch number to select.

    Returns
    -------
    idx : np.ndarray
        Integer array of simulation indexes in ``data.simlist`` that belong to
        the requested epoch.
    """
    idx = []
    for i, sim in enumerate(data.simlist):
        name = os.path.basename(os.path.dirname(os.path.abspath(sim.trajectory[0])))
        e = int(epochregex.findall(name)[0])
        if e == epoch:
            idx.append(i)

    return np.array(idx)


def getParentSimIdxFrame(data: "MetricData", trajidx: int) -> tuple:
    """Find the parent simulation index and frame for a given trajectory.

    For a simulation that was spawned from a parent trajectory, returns the
    index of the parent simulation within ``data.simlist`` and the frame
    index within that parent simulation.

    Parameters
    ----------
    data : MetricData
        Projected simulation data containing the simlist to search.
    trajidx : int
        Index of the simulation in ``data.simlist`` for which to find the parent.

    Returns
    -------
    parentidx : int
        Index of the parent simulation in ``data.simlist``.
    frameidx : int
        Frame index within the parent simulation corresponding to the spawn point.
    """
    name = os.path.basename(
        os.path.dirname(os.path.abspath(data.simlist[trajidx].trajectory[0]))
    )
    res = parentregex.findall(name)
    if len(res) == 0:
        print(
            "Failed to find parent simulation of {}. Assuming it's the first epoch. Will use the first frame as the state".format(
                name
            )
        )
        return trajidx, 0

    res = res[0]
    parentname, parentpiece, parentframe = res
    parentpiece = int(parentpiece)
    parentframe = int(parentframe)

    parentidx = None
    for i, sim in enumerate(data.simlist):
        name = os.path.basename(os.path.dirname(sim.trajectory[0]))
        if name.startswith(parentname + "_"):
            if parentidx is None:
                parentidx = i
            else:
                raise RuntimeError("More than one simulation matches parent name!")

    pieceframe = np.array([parentpiece, parentframe])
    frameidx = np.where(
        np.all(data.trajectories[parentidx].reference == pieceframe, axis=1)
    )[0]
    if len(frameidx) != 1:
        raise RuntimeError("Only one frame should match simpiece/frame combo")
    frameidx = frameidx[0]

    return parentidx, frameidx


def updatingMean(oldmean: float, oldcount: int, newdata: np.ndarray) -> float:
    """Compute a running mean incorporating new observations.

    Parameters
    ----------
    oldmean : float
        The previously accumulated mean value.
    oldcount : int
        The number of observations that produced ``oldmean``.
    newdata : np.ndarray
        Array of new observations to incorporate.

    Returns
    -------
    mean : float
        Updated mean over all observations (old and new combined).
    """
    if oldcount == 0:
        assert oldmean == 0
        return np.mean(newdata)
    newcount = oldcount + len(newdata)
    return (oldmean + np.sum(newdata) / oldcount) * (oldcount / newcount)
