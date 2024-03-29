"""
pdbreporter.py: Outputs simulation trajectories in PDB format

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Peter Eastman
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
from __future__ import absolute_import

__author__ = "Peter Eastman"
__version__ = "1.0"

from simtk.openmm.app import PDBFile, PDBxFile


class PDBReporter(object):
    """PDBReporter outputs a series of frames from a Simulation to a PDB file.

    To use it, create a PDBReporter, then add it to the Simulation's list of reporters.
    """

    def __init__(self, file, reportInterval, enforcePeriodicBox=True):
        """Create a PDBReporter.

        Parameters
        ----------
        file : string
            The file to write to
        reportInterval : int
            The interval (in time steps) at which to write frames
        """
        self._reportInterval = reportInterval
        self._enforcePeriodicBox = enforcePeriodicBox
        self._out = open(file, "w")
        self._topology = None
        self._nextModel = 0

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        tuple
            A five element tuple. The first element is the number of steps
            until the next report. The remaining elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.
        """
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, False, False, False)

    def report(self, simulation, _):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        _ : State
            The current state of the simulation
        """
        state = simulation.context.getState(
            getPositions=True, enforcePeriodicBox=self._enforcePeriodicBox
        )
        if self._nextModel == 0:
            PDBFile.writeHeader(simulation.topology, self._out)
            self._topology = simulation.topology
            self._nextModel += 1
        PDBFile.writeModel(
            simulation.topology, state.getPositions(), self._out, self._nextModel
        )
        self._nextModel += 1
        if hasattr(self._out, "flush") and callable(self._out.flush):
            self._out.flush()

    def __del__(self):
        if self._topology is not None:
            PDBFile.writeFooter(self._topology, self._out)
        self._out.close()


class PDBxReporter(PDBReporter):
    """PDBxReporter outputs a series of frames from a Simulation to a PDBx/mmCIF file.

    To use it, create a PDBxReporter, then add it to the Simulation's list of reporters.
    """

    def report(self, simulation, _):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        _ : State
            The current state of the simulation
        """
        state = simulation.context.getState(
            getPositions=True, enforcePeriodicBox=self._enforcePeriodicBox
        )
        if self._nextModel == 0:
            PDBxFile.writeHeader(simulation.topology, self._out)
            self._nextModel += 1
        PDBxFile.writeModel(
            simulation.topology, state.getPositions(), self._out, self._nextModel
        )
        self._nextModel += 1
        if hasattr(self._out, "flush") and callable(self._out.flush):
            self._out.flush()

    def __del__(self):
        self._out.close()
        PDBReporter.__del__(self)
