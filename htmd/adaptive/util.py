import numpy as np
import os
import re

parentregex = re.compile('e\d+s\d+_(e\d+s\d+)p(\d+)f(\d+)')
epochregex = re.compile('^e(\d+)s')


def getEpochTrajectoryDictionary(simlist):
    from collections import defaultdict

    simepochs = defaultdict(list)
    for i, sim in enumerate(simlist):
        name = os.path.basename(os.path.dirname(os.path.abspath(sim.trajectory[0])))
        e = int(epochregex.findall(name)[0])
        simepochs[e].append(i)

    return simepochs

def getEpochSimIdx(data, epoch):
    idx = []
    for i, sim in enumerate(data.simlist):
        name = os.path.basename(os.path.dirname(os.path.abspath(sim.trajectory[0])))
        e = int(epochregex.findall(name)[0])
        if e == epoch:
            idx.append(i)

    return np.array(idx)

def getParentSimIdxFrame(data, trajidx):
    name = os.path.basename(os.path.dirname(os.path.abspath(data.simlist[trajidx].trajectory[0])))
    res = parentregex.findall(name)
    if len(res) == 0:
        print('Failed to find parent simulation of {}. Assuming it\'s the first epoch. Will use the first frame as the state'.format(name))
        return trajidx, 0

    res = res[0]
    parentname, parentpiece, parentframe = res
    parentpiece = int(parentpiece)
    parentframe = int(parentframe)

    parentidx = None
    for i, sim in enumerate(data.simlist):
        name = os.path.basename(os.path.dirname(sim.trajectory[0]))
        if name.startswith(parentname+'_'):
            if parentidx is None:
                parentidx = i
            else:
                raise RuntimeError('More than one simulation matches parent name!')

    pieceframe = np.array([parentpiece, parentframe])
    frameidx = np.where(np.all(data.trajectories[parentidx].reference == pieceframe, axis=1))[0]
    if len(frameidx) != 1:
        raise RuntimeError('Only one frame should match simpiece/frame combo')
    frameidx = frameidx[0]

    return parentidx, frameidx

def updatingMean(oldmean, oldcount, newdata):
    if oldcount == 0:
        assert oldmean == 0
        return np.mean(newdata)
    newcount = oldcount + len(newdata)
    return (oldmean + np.sum(newdata) / oldcount) * (oldcount / newcount)


