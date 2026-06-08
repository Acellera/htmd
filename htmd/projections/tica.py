# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from typing import TYPE_CHECKING

import numpy as np
from htmd.projections.metric import Metric, _projectionGenerator
from htmd.units import convert as unitconvert
import logging

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from htmd.metricdata import MetricData


class TICA(object):
    """Class for calculating the TICA projections of a MetricData object.

    Time-based Independent Component Analysis projects your data on the slowest
    coordinates identified for a given lagtime.

    Parameters
    ----------
    data : :class:`MetricData <htmd.metricdata.MetricData>` object or :class:`Metric <htmd.projections.metric.Metric>` object
        The MetricData object whose data to project, or a Metric object for memory-efficient
        streaming TICA (projects trajectories on the fly).
    lag : float
        The correlation lagtime to use for TICA. Units are controlled by ``units``.
    units : str, optional
        The units of ``lag``. Can be ``'frames'`` or any time unit given as a string.
    dimensions : list, optional
        A list of dimensions of the original data on which to apply TICA. All other dimensions
        will stay unaltered. If None, TICA is applied on all dimensions.
    njobs : int, optional
        Number of jobs to spawn for parallel computation of TICA components. If None it will use
        the default from htmd.config.

    Examples
    --------
    >>> from htmd.projections.tica import TICA
    >>> metr = Metric(sims)
    >>> metr.set(MetricSelfDistance('protein and name CA'))
    >>> data = metr.project()
    >>> tica = TICA(data, 20)
    >>> datatica = tica.project(3)

    References
    ----------
    Perez-Hernandez, G. and Paul, F. and Giorgino, T. and de Fabritiis, G.
    and Noe, F. (2013) Identification of slow molecular order parameters
    for Markov model construction. J. Chem. Phys., 139 . 015102.
    """

    def __init__(
        self,
        data: "MetricData | Metric",
        lag: float,
        units: str = "frames",
        dimensions: list | range | np.ndarray | None = None,
        njobs: int | None = None,
    ):
        from deeptime.decomposition import TICA as TICAdt
        from tqdm import tqdm
        from htmd.util import _getNjobs

        self.data = data
        self.dimensions = dimensions
        self.njobs = njobs if njobs is not None else _getNjobs()

        if isinstance(
            data, Metric
        ):  # Memory efficient TICA projecting trajectories on the fly
            if units != "frames":
                raise RuntimeError(
                    "Cannot use delayed projection TICA with units other than frames for now. Report this to HTMD issues."
                )
            tic = TICAdt(lagtime=lag)
            metr = data

            pbar = tqdm(total=len(metr.simulations))
            for proj in _projectionGenerator(metr, self.njobs):
                for pro in proj:
                    if pro is None:
                        continue
                    if self.dimensions is None:
                        tic.partial_fit(pro[0])
                    else:  # Sub-select dimensions for fitting
                        tic.partial_fit(pro[0][:, self.dimensions])
                pbar.update(len(proj))
            pbar.close()
        else:  # In-memory TICA
            lag = unitconvert(units, "frames", lag, data.fstep)
            if lag == 0:
                raise RuntimeError(
                    "Lag time conversion resulted in 0 frames. Please use a larger lag-time for TICA."
                )

            tic = TICAdt(lagtime=lag)
            if self.dimensions is None:
                datalist = data.dat
            else:  # Sub-select dimensions for fitting
                datalist = [x[:, self.dimensions].copy() for x in data.dat]
            tic.fit(datalist)
        self.model = tic.fetch_model()

    def project(
        self,
        ndim: int | None = None,
        var_cutoff: float = 0.95,
    ) -> "MetricData":
        """Project the data object given to the constructor onto the top TICA dimensions.

        Parameters
        ----------
        ndim : int, optional
            The number of TICA dimensions to project the data on. If None, ``var_cutoff`` is used
            to determine the number of dimensions automatically.
        var_cutoff : float, optional
            Variance cutoff used for automatically determining the number of dimensions.

        Returns
        -------
        dataTica : :class:`MetricData <htmd.metricdata.MetricData>` object
            A new MetricData object containing the TICA projected data.

        Examples
        --------
        >>> from htmd.projections.tica import TICA
        >>> tica = TICA(data, 20)
        >>> dataTica = tica.project(5)
        """
        from tqdm import tqdm

        if ndim is not None:
            self.model.dim = ndim
            # Replace the following lines in future deeptime versions with self.model.var_cutoff = None
            self.model._var_cutoff = None
            self.model._update_output_dimension()
            # End of hack
        elif var_cutoff is not None:
            self.model.var_cutoff = var_cutoff

        keepdata = []
        keepdim = None
        keepdimdesc = None
        if isinstance(
            self.data, Metric
        ):  # Memory efficient TICA projecting trajectories on the fly
            proj = []
            ref = []
            fstep = None

            metr = self.data
            k = -1
            droppedsims = []
            pbar = tqdm(total=len(metr.simulations))
            for projecteddata in _projectionGenerator(metr, self.njobs):
                for pro in projecteddata:
                    k += 1
                    if pro is None:
                        droppedsims.append(k)
                        continue
                    if self.dimensions is not None:
                        numDimensions = pro[0].shape[1]
                        keepdim = np.setdiff1d(range(numDimensions), self.dimensions)
                        keepdata.append(pro[0][:, keepdim])
                        proj.append(
                            self.model.transform(pro[0][:, self.dimensions]).astype(
                                np.float32
                            )
                        )  # Sub-select dimensions for projecting
                    else:
                        proj.append(self.model.transform(pro[0]).astype(np.float32))
                    ref.append(pro[1])
                    if fstep is None:
                        fstep = pro[2]
                pbar.update(len(projecteddata))
            pbar.close()

            simlist = self.data.simulations
            simlist = np.delete(simlist, droppedsims)
            parent = None
            if self.dimensions is not None:
                from htmd.projections.metric import _singleMolfile
                from moleculekit.molecule import Molecule

                (single, molfile) = _singleMolfile(metr.simulations)
                if single:
                    keepdimdesc = metr.getMapping(Molecule(molfile))
                    keepdimdesc = keepdimdesc.iloc[keepdim]
        else:
            if ndim is not None and self.data.numDimensions < ndim:
                raise RuntimeError(
                    f"TICA cannot increase the dimensionality of your data. Your data has {self.data.numDimensions} dimensions and you requested {ndim} TICA dimensions"
                )

            if self.dimensions is not None:
                keepdim = np.setdiff1d(range(self.data.numDimensions), self.dimensions)
                keepdata = [x[:, keepdim] for x in self.data.dat]
                if self.data.description is not None:
                    keepdimdesc = self.data.description.iloc[keepdim]
                # The model was fitted only on the requested subset, so transform
                # only those dimensions (the rest are carried through unaltered).
                proj = [
                    self.model.transform(tr.projection[:, self.dimensions])
                    for tr in self.data.trajectories
                ]
            else:
                proj = [
                    self.model.transform(tr.projection)
                    for tr in self.data.trajectories
                ]
            simlist = self.data.simlist
            ref = self.data.ref
            fstep = self.data.fstep
            parent = self.data

        # If TICA is done on a subset of dimensions, combine non-projected data with projected data
        if self.dimensions is not None:
            newproj = []
            for k, t in zip(keepdata, proj):
                newproj.append(np.hstack((k, t)))
            proj = newproj

        if ndim is None:
            ndim = proj[0].shape[1]
            logger.info(
                f"Kept {ndim} dimension(s) to cover {var_cutoff * 100:.1f}% of kinetic variance."
            )

        from htmd.metricdata import MetricData

        datatica = MetricData(
            dat=proj, simlist=simlist, ref=ref, fstep=fstep, parent=parent
        )
        from pandas import DataFrame, concat

        types = []
        indexes = []
        description = []
        for i in range(ndim):
            types += ["tica"]
            indexes += [-1]
            description += [f"TICA dimension {i + 1}"]
        datatica.description = DataFrame(
            {"type": types, "atomIndexes": indexes, "description": description}
        )

        if (
            self.dimensions is not None and keepdimdesc is not None
        ):  # If TICA is done on a subset of dims
            datatica.description = concat(
                [keepdimdesc, datatica.description], ignore_index=True
            )

        return datatica


if __name__ == "__main__":
    from htmd.simlist import simlist
    from glob import glob
    from moleculekit.projections.metricdistance import MetricSelfDistance
    from htmd.home import home
    from os.path import join

    testfolder = home(dataDir="villin")

    sims = simlist(glob(join(testfolder, "*", "")), join(testfolder, "filtered.pdb"))
    met = Metric(sims[0:2])
    met.set(MetricSelfDistance("protein and name CA"))
    data = met.project()
    data.fstep = 0.1

    tica = TICA(data, 2, dimensions=range(2, 10))
    datatica = tica.project(2)
    tica5 = TICA(data, 0.2, units="ns", dimensions=range(2, 10))
    datatica5 = tica5.project(2)
    expected = [
        [3.69098878, -0.33862674, 0.85779184],
        [3.77816105, -0.31887317, 0.87724227],
        [3.83537507, -0.11878026, 0.65236956],
    ]
    assert np.allclose(
        np.abs(datatica.trajectories[0].projection[-3:, -3:]),
        np.abs(np.array(expected, dtype=np.float32)),
        rtol=0,
        atol=0.01,
    )
    assert np.allclose(
        np.abs(datatica5.trajectories[0].projection[-3:, -3:]),
        np.abs(np.array(expected, dtype=np.float32)),
        rtol=0,
        atol=0.01,
    )
    assert np.all(datatica.description.iloc[[587, 588]].type == "tica")
    assert np.all(datatica.description.iloc[range(587)].type == "distance")
    print("In-memory TICA with subset of dimensions passed test.")

    tica2 = TICA(met, 2, dimensions=range(2, 10))
    datatica2 = tica2.project(2)
    assert np.allclose(
        np.abs(datatica2.trajectories[0].projection[-3:, -3:]),
        np.abs(np.array(expected, dtype=np.float32)),
        rtol=0,
        atol=0.01,
    )
    assert np.all(datatica2.description.iloc[[587, 588]].type == "tica")
    assert np.all(datatica2.description.iloc[range(587)].type == "distance")
    print("Streaming TICA with subset of dimensions passed test.")

    # assert np.max(np.abs(datatica.dat[0][:, -2:]) - np.abs(datatica2.dat[0][:, -2:])) < 0.01, 'Streaming and memory subdim TICA inconsistent.'

    tica3 = TICA(data, 2)
    datatica3 = tica3.project(2)
    expected = [
        [-1.36328638, -0.35354128],
        [-1.35348749, -0.13028328],
        [-1.43249917, -0.31004715],
    ]
    assert np.allclose(
        np.abs(datatica3.trajectories[0].projection[-3:, :]),
        np.abs(np.array(expected, dtype=np.float32)),
        rtol=0,
        atol=0.01,
    )
    assert np.all(datatica3.description.iloc[[0, 1]].type == "tica")
    print("In-memory TICA passed test.")

    tica4 = TICA(met, 2)
    datatica4 = tica4.project(2)
    assert np.allclose(
        np.abs(datatica4.trajectories[0].projection[-3:, :]),
        np.abs(np.array(expected, dtype=np.float32)),
        rtol=0,
        atol=0.01,
    )
    assert np.all(datatica4.description.iloc[[0, 1]].type == "tica")
    print("Streaming TICA passed test.")

    assert (
        np.max(
            np.abs(datatica4.trajectories[0].projection)
            - np.abs(datatica3.trajectories[0].projection)
        )
        < 0.01
    ), "Streaming and memory TICA inconsistent."
