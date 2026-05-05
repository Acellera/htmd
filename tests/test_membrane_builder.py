import os
import pytest

from htmd.builder.amber import _findTeLeap
from htmd.membranebuilder.build_membrane import buildMembrane


reason = "teLeap is not installed. Cannot test membrane builder build/equil paths"
tleap_installed = _findTeLeap() is not None


def _test_subrandom_positions_scale_each_axis():
    """Each XY axis must be scaled by its own box edge, not the last one.

    Regression test for a bug where _subrandom_particle_positions scaled all
    columns by box_vectors[-1][-1], so non-square membranes had X positions
    in [-Ly/2, Ly/2] instead of [-Lx/2, Lx/2].
    """
    from openmm import unit
    from htmd.membranebuilder.ljfluid import _subrandom_particle_positions

    Lx, Ly, Lz = 50.0, 100.0, 80.0
    box = [
        unit.Quantity((Lx * unit.angstrom, 0 * unit.angstrom, 0 * unit.angstrom)),
        unit.Quantity((0 * unit.angstrom, Ly * unit.angstrom, 0 * unit.angstrom)),
        unit.Quantity((0 * unit.angstrom, 0 * unit.angstrom, Lz * unit.angstrom)),
    ]
    pos = _subrandom_particle_positions(200, box, 2).value_in_unit(unit.angstrom)

    assert pos[:, 0].min() >= -Lx / 2 and pos[:, 0].max() <= Lx / 2
    assert pos[:, 1].min() >= -Ly / 2 and pos[:, 1].max() <= Ly / 2
    assert pos[:, 0].max() - pos[:, 0].min() > 0.5 * Lx
    assert pos[:, 1].max() - pos[:, 1].min() > 0.5 * Ly
    assert (pos[:, 2] == 0).all()


def _test_build_membrane(tmp_path):
    buildMembrane(
        [20, 20],
        ratioupper={"popc": 0.42, "pope": 0.4, "chl1": 0.18},
        ratiolower={"popc": 0.42, "pope": 0.4, "chl1": 0.18},
        equilibrate=0,
        minimize=0,
        outdir=str(tmp_path),
        platform="CPU",
    )
    assert os.path.exists(tmp_path / "structure.pdb")
    assert not os.path.exists(tmp_path / "starting_structure.pdb")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def _test_build_membrane_minimize(tmp_path):
    buildMembrane(
        [20, 20],
        ratioupper={"popc": 0.42, "pope": 0.4, "chl1": 0.18},
        ratiolower={"popc": 0.42, "pope": 0.4, "chl1": 0.18},
        equilibrate=0,
        minimize=100,
        outdir=str(tmp_path),
        platform="CPU",
    )
    assert os.path.exists(tmp_path / "structure.pdb")
    assert os.path.exists(tmp_path / "starting_structure.pdb")


@pytest.mark.skipif(not tleap_installed, reason=reason)
def _test_build_membrane_equil(tmp_path):
    buildMembrane(
        [20, 20],
        ratioupper={"popc": 0.42, "pope": 0.4, "chl1": 0.18},
        ratiolower={"popc": 0.42, "pope": 0.4, "chl1": 0.18},
        equilibrate=0.01,
        minimize=100,
        outdir=str(tmp_path),
        platform="CPU",
    )
    assert os.path.exists(tmp_path / "structure.pdb")
    assert os.path.exists(tmp_path / "starting_structure.pdb")
