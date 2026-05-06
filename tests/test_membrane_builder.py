import os

import numpy as np

from htmd.membranebuilder.build_membrane import buildMembrane


PLATFORM = "CPU"


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
        platform=PLATFORM,
        seed=42,
    )
    assert os.path.exists(tmp_path / "structure.pdb")
    assert not os.path.exists(tmp_path / "starting_structure.pdb")


def _test_build_membrane_minimize(tmp_path):
    buildMembrane(
        [20, 20],
        ratioupper={"popc": 0.42, "pope": 0.4, "chl1": 0.18},
        ratiolower={"popc": 0.42, "pope": 0.4, "chl1": 0.18},
        equilibrate=0,
        minimize=100,
        outdir=str(tmp_path),
        platform=PLATFORM,
        seed=42,
    )
    assert os.path.exists(tmp_path / "structure.pdb")
    assert os.path.exists(tmp_path / "starting_structure.pdb")


def _test_build_membrane_equil(tmp_path):
    buildMembrane(
        [20, 20],
        ratioupper={"popc": 0.42, "pope": 0.4, "chl1": 0.18},
        ratiolower={"popc": 0.42, "pope": 0.4, "chl1": 0.18},
        equilibrate=0.005,
        minimize=100,
        outdir=str(tmp_path),
        platform=PLATFORM,
        seed=42,
    )
    assert os.path.exists(tmp_path / "structure.pdb")
    assert os.path.exists(tmp_path / "starting_structure.pdb")


# def _test_build_membrane_with_solute_2kdc(tmp_path):
#     """End-to-end build around OPM-aligned 2KDC.

#     Asserts:
#     - Total lipid count is reduced vs. the no-solute baseline (count
#       reduction in ``_createLipids`` plus the final overlap drop).
#     - No membrane heavy atom comes within 2 A of any solute heavy atom
#       (combined effect of Halton filter, LJ obstacles, and overlap drop).
#     - The ratio drift between requested and resulting per-resname counts
#       is small (no per-resname drop fraction exceeds 25%).
#     """
#     from moleculekit.opm import get_opm_pdb
#     from scipy.spatial import cKDTree

#     ref, _ = get_opm_pdb("2kdc", validateElements=False)

#     ratios = {"popc": 0.5, "pope": 0.3, "chl1": 0.2}

#     membrane = buildMembrane(
#         [70, 70],
#         ratioupper=ratios,
#         ratiolower=ratios,
#         equilibrate=0.5,
#         minimize=500,
#         outdir=str(tmp_path),
#         platform=PLATFORM,
#         seed=42,
#         solute=ref,
#     )

#     heavy = ref.element != "H"
#     solute_xyz = ref.coords[heavy, :, 0]
#     is_lipid = (membrane.resname != "TIP3") & (membrane.resname != "HOH")
#     is_heavy = membrane.element != "H"
#     memb_xyz = membrane.coords[is_lipid & is_heavy, :, 0]

#     tree = cKDTree(solute_xyz)
#     dists, _ = tree.query(memb_xyz, k=1)
#     assert (
#         dists.min() >= 2.0
#     ), f"membrane heavy atom within 2 A of solute: min={dists.min():.2f}"


def _test_build_membrane_with_solute_2kdc_lj_only(tmp_path):
    """Run buildMembrane around 2KDC up to and including the LJ packing
    (no minimization or equilibration). Verify the lj_packing.pdb debug
    file has both lipid heads and obstacles in the user's coordinate frame.
    """
    from moleculekit.molecule import Molecule
    from moleculekit.opm import get_opm_pdb

    ref, _ = get_opm_pdb("2kdc", validateElements=False)

    buildMembrane(
        [70, 70],
        ratioupper={"popc": 0.5, "pope": 0.3, "chl1": 0.2},
        ratiolower={"popc": 0.5, "pope": 0.3, "chl1": 0.2},
        equilibrate=0,
        minimize=0,
        outdir=str(tmp_path),
        platform=PLATFORM,
        seed=42,
        solute=ref,
    )

    pdb_path = tmp_path / "lj_packing.pdb"
    assert pdb_path.exists()

    debug = Molecule(str(pdb_path))
    assert (debug.resname == "OBS").sum() > 0, "no obstacles in lj_packing.pdb"
    assert (debug.resname != "OBS").sum() > 0, "no lipid heads in lj_packing.pdb"

    # Lipid heads and obstacles should both be in the solute's xy frame.
    heavy = ref.element != "H"
    sol_com_xy = ref.coords[heavy, :2, 0].mean(axis=0)
    head_com_xy = debug.coords[debug.resname != "OBS", :2, 0].mean(axis=0)
    obs_com_xy = debug.coords[debug.resname == "OBS", :2, 0].mean(axis=0)
    assert np.allclose(head_com_xy, sol_com_xy, atol=10.0), (
        f"lipid heads not aligned with solute: heads={head_com_xy}, solute={sol_com_xy}"
    )
    assert np.allclose(obs_com_xy, sol_com_xy, atol=10.0), (
        f"obstacles not aligned with solute: obs={obs_com_xy}, solute={sol_com_xy}"
    )


def _test_solute_footprint_2kdc():
    """Footprint and area-fraction of OPM-aligned 2KDC (M2 TM peptide).

    The OPM-distributed structure already has bilayer center at z=0 and the
    TM axis aligned with z, so we can sample the leaflet head planes directly
    at +-thickness/2. With the default offset (4 A) the slab samples the
    membrane-embedded TM bundle, missing both the upper extramembrane region
    and the lower amphipathic helix that lies flat at the membrane surface.
    The two leaflets should give comparable fractions for this protein.
    """
    from moleculekit.opm import get_opm_pdb
    from htmd.membranebuilder.build_membrane import (
        _solute_footprint,
        _solute_area_fraction,
    )

    np.random.seed(0)
    ref, thickness = get_opm_pdb("2kdc", validateElements=False)
    head_z = thickness / 2

    fp_u = _solute_footprint(ref, head_z)
    fp_l = _solute_footprint(ref, -head_z)
    assert fp_u is not None and fp_l is not None

    f_u = _solute_area_fraction(fp_u, [60, 60])
    f_l = _solute_area_fraction(fp_l, [60, 60])
    # The default buffer (5 A) is added to vdw radii so the forbidden disks
    # are large enough to keep lipid head centers outside protein pores.
    assert 0.45 < f_u < 0.70, f"upper fraction out of range: {f_u}"
    assert 0.60 < f_l < 0.80, f"lower fraction out of range: {f_l}"
    # The TM bundle is broadly symmetric but the per-disk overlap differs
    # between leaflets (lower atoms are more spread), so absolute fractions
    # differ more than at zero buffer.
    assert abs(f_l - f_u) < 0.20, f"unexpectedly asymmetric: u={f_u} l={f_l}"
