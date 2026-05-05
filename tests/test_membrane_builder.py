import os
import tempfile
import pytest

from htmd.builder.amber import _findTeLeap
from htmd.membranebuilder.build_membrane import buildMembrane


reason = "teLeap is not installed. Cannot test membrane builder build/equil paths"
tleap_installed = _findTeLeap() is not None


def _test_build_membrane():
    with tempfile.TemporaryDirectory() as tmpdir:
        buildMembrane(
            [20, 20],
            ratioupper={"popc": 0.42, "pope": 0.4, "chl1": 0.18},
            ratiolower={"popc": 0.42, "pope": 0.4, "chl1": 0.18},
            equilibrate=0,
            minimize=0,
            outdir=tmpdir,
            platform="CPU",
        )
        assert os.path.exists(os.path.join(tmpdir, "structure.pdb"))
        assert not os.path.exists(os.path.join(tmpdir, "starting_structure.pdb"))


@pytest.mark.skipif(not tleap_installed, reason=reason)
def _test_build_membrane_minimize():
    with tempfile.TemporaryDirectory() as tmpdir:
        buildMembrane(
            [20, 20],
            ratioupper={"popc": 0.42, "pope": 0.4, "chl1": 0.18},
            ratiolower={"popc": 0.42, "pope": 0.4, "chl1": 0.18},
            equilibrate=0,
            minimize=100,
            outdir=tmpdir,
            platform="CPU",
        )
        assert os.path.exists(os.path.join(tmpdir, "structure.pdb"))
        assert os.path.exists(os.path.join(tmpdir, "starting_structure.pdb"))


@pytest.mark.skipif(not tleap_installed, reason=reason)
def _test_build_membrane_equil():
    with tempfile.TemporaryDirectory() as tmpdir:
        buildMembrane(
            [20, 20],
            ratioupper={"popc": 0.42, "pope": 0.4, "chl1": 0.18},
            ratiolower={"popc": 0.42, "pope": 0.4, "chl1": 0.18},
            equilibrate=0.01,
            minimize=100,
            outdir=tmpdir,
            platform="CPU",
        )
        assert os.path.exists(os.path.join(tmpdir, "structure.pdb"))
        assert os.path.exists(os.path.join(tmpdir, "starting_structure.pdb"))
