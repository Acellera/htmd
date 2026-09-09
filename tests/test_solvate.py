import itertools
import os

import numpy as np
import pytest

from moleculekit.molecule import Molecule

from htmd.builder.solvate import solvate

curr_dir = os.path.dirname(os.path.abspath(__file__))


def _shell_spectrum(vectors, nmax=4, rmax_factor=1.6):
    """Rotation-invariant lattice fingerprint: (distance, multiplicity) shells.

    Two unit cells describe the same lattice iff their shell spectra agree.
    Comparing box vectors directly does not work, because two representatives
    of one lattice differ by a rotation as well as a change of basis.
    """
    idx = np.array(list(itertools.product(range(-nmax, nmax + 1), repeat=3)))
    r = np.linalg.norm(idx @ vectors, axis=1)
    rmax = rmax_factor * min(np.linalg.norm(vectors, axis=1))
    r = np.sort(r[(r > 1e-9) & (r < rmax)])
    shells, cur = [], [r[0]]
    for x in r[1:]:
        if x - cur[-1] < 1e-4:
            cur.append(x)
        else:
            shells.append((round(float(np.mean(cur)), 3), len(cur)))
            cur = [x]
    shells.append((round(float(np.mean(cur)), 3), len(cur)))
    return shells


def _water_triples(mol):
    """Return the (nwat, 3, 3) coordinate array of the solvate-added waters.

    Waters are emitted as consecutive O, H1, H2 triples.
    """
    watidx = np.where(mol.atomselect("water"))[0]
    assert watidx.size % 3 == 0, "water atom count is not a multiple of 3"
    return mol.coords[watidx, :, 0].reshape(-1, 3, 3), watidx


def _solute():
    """A tiny solute so solvate() has something to build a box around."""
    mol = Molecule().empty(2)
    mol.name[:] = ["C1", "C2"]
    mol.element[:] = ["C", "C"]
    mol.resname[:] = "LIG"
    mol.resid[:] = [1, 1]
    mol.segid[:] = "S0"
    mol.coords = np.array(
        [[[0.0], [0.0], [0.0]], [[1.5], [0.0], [0.0]]], dtype=np.float32
    )
    return mol


def test_solvate_returns_bonded_water():
    """Every solvate-added water must carry its two O-H bonds.

    Without them, anything that groups atoms by connectivity (notably
    Molecule.wrap) treats each water atom as an independent molecule.
    """
    smol = solvate(_solute(), pad=8)

    _, watidx = _water_triples(smol)
    nwat = watidx.size // 3
    assert nwat > 0, "solvate added no water"

    watset = set(watidx.tolist())
    watbonds = smol.bonds[np.isin(smol.bonds, list(watset)).any(axis=1)]
    assert (
        watbonds.shape[0] == 2 * nwat
    ), f"expected {2 * nwat} water bonds for {nwat} waters, got {watbonds.shape[0]}"

    # No bond may join two different residues (i.e. two different waters)
    resids = smol.resid[watbonds]
    segids = smol.segid[watbonds]
    assert np.all(resids[:, 0] == resids[:, 1]), "water bond crosses residues"
    assert np.all(segids[:, 0] == segids[:, 1]), "water bond crosses segments"


def test_solvated_water_survives_wrapping():
    """Regression: waters straddling a periodic face must wrap as whole molecules.

    Bond-less water made Molecule.wrap() translate individual O/H atoms across
    the boundary, producing intramolecular "bonds" the length of a box edge.
    This is what corrupted membrane builds.
    """
    smol = solvate(_solute(), pad=8)

    coords = smol.coords[:, :, 0]
    smol.box = (
        (coords.max(axis=0) - coords.min(axis=0)).astype(np.float32).reshape(3, 1)
    )
    smol.boxangles = np.array([[90.0], [90.0], [90.0]], dtype=np.float32)

    # Wrap around a corner of the box so plenty of waters straddle a face.
    smol.wrap(wrapcenter=coords.min(axis=0))

    triples, _ = _water_triples(smol)
    oh1 = np.linalg.norm(triples[:, 1] - triples[:, 0], axis=1)
    oh2 = np.linalg.norm(triples[:, 2] - triples[:, 0], axis=1)
    worst = max(oh1.max(), oh2.max())
    assert (
        worst < 1.5
    ), f"water split across the periodic boundary: max O-H = {worst:.2f} A"


def test_custom_bondless_solvent_file_still_gets_bonds(tmp_path):
    """A user-supplied PDB water box has no bonds; solvate must supply them."""
    from htmd.home import home

    ref = Molecule(os.path.join(home(shareDir=True), "solvate", "wat.bcif.gz"))
    assert ref.bonds.shape[0] > 0, "the shipped water box should carry bonds"

    # A hand-made water box PDB carries no CONECT records, which is the case
    # a user hits when passing their own solvent file.
    ref.bonds = np.empty((0, 2), dtype=np.uint32)
    ref.bondtype = np.empty(0, dtype=object)
    custom = str(tmp_path / "mywat.pdb")
    ref.write(custom)
    assert Molecule(custom).bonds.shape[0] == 0, "fixture is expected to be bond-less"

    smol = solvate(_solute(), pad=8, spdb=custom)

    _, watidx = _water_triples(smol)
    nwat = watidx.size // 3
    watbonds = smol.bonds[np.isin(smol.bonds, watidx.tolist()).any(axis=1)]
    assert (
        watbonds.shape[0] == 2 * nwat
    ), f"expected {2 * nwat} water bonds, got {watbonds.shape[0]}"


@pytest.mark.parametrize(
    "shape,volume_factor,angle",
    [
        ("cube", 1.0, 90.0),
        ("octahedron", 0.7698, 109.4712206),
        ("dodecahedron", 0.7071, 60.0),
    ],
)
def test_cell_vectors_volume_and_angles(shape, volume_factor, angle):
    """Each shape is equilateral: a=b=c and alpha=beta=gamma.

    Equal angles are what lets the AMBER prmtop's single-angle BOX_DIMENSIONS
    field store the cell losslessly.
    """
    from htmd.builder.solvate import _cell_vectors, _cell_lengths_and_angles

    d = 30.0
    vectors = _cell_vectors(shape, d)
    assert vectors.shape == (3, 3)
    assert abs(np.linalg.det(vectors)) == pytest.approx(volume_factor * d**3, rel=1e-4)

    lengths, angles = _cell_lengths_and_angles(vectors)
    assert lengths == pytest.approx([d, d, d], rel=1e-9)
    assert angles == pytest.approx([angle, angle, angle], abs=1e-6)


def test_cell_vectors_rejects_rectangular_and_unknown():
    from htmd.builder.solvate import _cell_vectors

    with pytest.raises(ValueError, match="rectangular"):
        _cell_vectors("rectangular", 30.0)
    with pytest.raises(ValueError, match="banana"):
        _cell_vectors("banana", 30.0)


@pytest.mark.parametrize("shape", ["cube", "octahedron", "dodecahedron"])
def test_cell_vectors_are_accepted_by_openmm(shape):
    """Guard for the reduced-form float boundary.

    The 60/60/60 and 109.47 cells sit exactly on OpenMM's reduction limit
    (2*|b0| <= a0), so building them via a lengths-and-angles conversion
    produces 1.5000000000000004 and is rejected. This test fails if anyone
    swaps the exact construction for that conversion.
    """
    import openmm
    import openmm.unit as unit
    from openmm import Vec3

    from htmd.builder.solvate import _cell_vectors

    widths = list(np.arange(20.0, 200.0, 0.5)) + [37.317, 61.803, 88.8888, 99.999]
    for d in widths:
        vectors = _cell_vectors(shape, d)
        system = openmm.System()
        system.setDefaultPeriodicBoxVectors(
            *[Vec3(*(v / 10.0)) * unit.nanometer for v in vectors]
        )


@pytest.mark.parametrize(
    "shape,expected_shells",
    [
        ("octahedron", [(30.0, 8), (34.641, 6)]),
        ("dodecahedron", [(30.0, 12), (42.426, 6)]),
    ],
)
def test_cell_vectors_match_openmm_lattice(shape, expected_shells):
    """The equilateral representative must be the same lattice OpenMM uses.

    Octahedron is BCC (8 neighbours at d, 6 at 1.155d), dodecahedron is FCC
    (12 at d, 6 at 1.414d). These are the shells of OpenMM's own
    Modeller._computeBoxVectors forms, which use different angles.
    """
    import math

    from htmd.builder.solvate import _cell_vectors

    d = 30.0
    s2, s6 = math.sqrt(2), math.sqrt(6)
    openmm_form = {
        "octahedron": np.array(
            [[d, 0, 0], [d / 3, 2 * s2 * d / 3, 0], [-d / 3, s2 * d / 3, s6 * d / 3]]
        ),
        "dodecahedron": np.array([[d, 0, 0], [0, d, 0], [d / 2, d / 2, s2 * d / 2]]),
    }[shape]

    ours = _cell_vectors(shape, d)
    assert _shell_spectrum(ours) == expected_shells
    assert _shell_spectrum(ours) == _shell_spectrum(openmm_form)
    assert abs(np.linalg.det(ours)) == pytest.approx(abs(np.linalg.det(openmm_form)))


def _big_solute():
    """A solute large enough that 2*radius dominates the width floor.

    The tiny `_solute()` fixture has radius 0.75, so any useful pad triggers
    the `4*pad` floor and the resulting box is far too small for water density
    to be meaningful. Water density on HTMD's tiling path is strongly
    box-size dependent: measured on the unmodified rectangular path with a
    negligible solute, it runs 0.815 at width 16.5 A, 0.868 at 24.5 A, 0.903
    at 32.5 A and 0.951 at 60.5 A, asymptoting near 0.96. That is a skin
    effect of `_outOfBoundaries` dropping whole water residues at the cell
    faces, and it predates this work.
    """
    import numpy as np
    from moleculekit.molecule import Molecule

    mol = Molecule().empty(2)
    mol.name[:] = ["C1", "C2"]
    mol.element[:] = ["C", "C"]
    mol.resname[:] = "LIG"
    mol.resid[:] = [1, 1]
    mol.segid[:] = "S0"
    # 50 A apart, so radius is 25 and pad=5 gives width 60, not the 4*pad floor
    mol.coords = np.array(
        [[[0.0], [0.0], [0.0]], [[50.0], [0.0], [0.0]]], dtype=np.float32
    )
    return mol


@pytest.mark.parametrize(
    "shape,volume_factor",
    [("cube", 1.0), ("octahedron", 0.7698), ("dodecahedron", 0.7071)],
)
def test_solvate_shape_sets_cell_and_density(shape, volume_factor):
    """The returned cell must match the shape, and the water must be dense.

    The 0.90 floor needs a box wider than about 32 A, so this uses
    `_big_solute()` with pad=5 to land near 60 A. See that fixture's docstring
    for the measured size dependence.
    """
    import numpy as np

    from htmd.builder.solvate import solvate, _cell_vectors

    smol = solvate(_big_solute(), pad=5, shape=shape)

    lengths = smol.box[:, 0]
    assert np.allclose(lengths, lengths[0], rtol=1e-4), lengths
    assert float(lengths[0]) == pytest.approx(60.0, abs=1e-3)

    expected = _cell_vectors(shape, float(lengths[0]))
    volume = abs(np.linalg.det(expected))
    assert volume == pytest.approx(volume_factor * float(lengths[0]) ** 3, rel=1e-3)

    _, watidx = _water_triples(smol)
    nwat = watidx.size // 3
    density = nwat * 18.0153 / 6.02214076e23 / (volume * 1e-24)
    assert 0.90 < density < 1.05, f"{shape}: water density {density:.3f} g/cm3"

    assert smol.crystalinfo is not None
    for key, value in zip(
        ("a", "b", "c", "alpha", "beta", "gamma"),
        list(smol.box[:, 0]) + list(smol.boxangles[:, 0]),
    ):
        assert smol.crystalinfo[key] == pytest.approx(float(value), abs=1e-3)


def test_solvate_shape_saves_the_advertised_water():
    """The feature's actual claim: 77.0% and 70.7% of a cube's water.

    This is the size-independent assertion. Absolute density depends on box
    size through the boundary-cull skin effect, but the ratio between shapes
    at one width is fixed by the cell volumes.
    """
    from htmd.builder.solvate import solvate

    counts = {}
    for shape in ("cube", "octahedron", "dodecahedron"):
        smol = solvate(_big_solute(), pad=5, shape=shape)
        counts[shape] = int(smol.atomselect("water").sum()) // 3

    assert counts["octahedron"] / counts["cube"] == pytest.approx(0.7698, rel=0.04)
    assert counts["dodecahedron"] / counts["cube"] == pytest.approx(0.7071, rel=0.04)


@pytest.mark.parametrize("shape", ["cube", "octahedron", "dodecahedron"])
def test_solvate_shape_water_lies_inside_the_cell(shape):
    """Every water must sit inside the Wigner-Seitz cell.

    The center is taken from the solute atoms, not from the solvated bounding
    box, so the half-space test can be pinned exactly rather than with a loose
    tolerance.
    """
    import numpy as np

    from htmd.builder.solvate import solvate, _cell_vectors, _ws_halfspaces

    solute = _big_solute()
    coords = solute.coords[:, :, 0]
    center = 0.5 * (coords.min(axis=0) + coords.max(axis=0))

    smol = solvate(solute, pad=5, shape=shape)
    vectors = _cell_vectors(shape, float(smol.box[0, 0]))
    normals, offsets = _ws_halfspaces(vectors)

    triples, _ = _water_triples(smol)
    dist = (triples.reshape(-1, 3) - center) @ normals.T
    worst = float((dist - offsets).max())
    assert worst <= 1e-3, f"{shape}: water outside the cell by {worst:.4f} A"


@pytest.mark.parametrize("shape", ["cube", "octahedron", "dodecahedron"])
def test_ws_inradius_is_half_the_width(shape):
    """The sizing formula assumes the region's inscribed diameter is `width`.

    That holds for the Wigner-Seitz cell and not for the primitive
    parallelepiped, whose inradius is 0.40825*width for the two skewed shapes.
    If this ever fails, `width = 2*radius + 2*pad` stops delivering `pad`.
    """
    import numpy as np

    from htmd.builder.solvate import _cell_vectors, _ws_halfspaces

    width = 60.0
    normals, offsets = _ws_halfspaces(_cell_vectors(shape, width))
    inradius = float((offsets / np.linalg.norm(normals, axis=1)).min())
    assert inradius == pytest.approx(width / 2.0, rel=1e-9)


def _shell_solute(radius=25.0, natoms=400, seed=0):
    """A hollow sphere of atoms, the adversarial case for cell culling.

    `_big_solute()` is a straight line along the a-axis, so its atoms never
    reach the short faces of a skewed cell and it cannot detect a cull that
    carves the wrong region. A shell reaches every direction at once, and is
    closer to how a real protein's atoms are distributed.
    """
    import numpy as np
    from moleculekit.molecule import Molecule

    rng = np.random.default_rng(seed)
    vecs = rng.normal(size=(natoms, 3))
    vecs /= np.linalg.norm(vecs, axis=1)[:, None]

    mol = Molecule().empty(natoms)
    mol.name[:] = ["C"] * natoms
    mol.element[:] = ["C"] * natoms
    mol.resname[:] = "LIG"
    mol.resid[:] = np.arange(natoms)
    mol.segid[:] = "S0"
    mol.coords = (vecs * radius).astype(np.float32).reshape(natoms, 3, 1)
    return mol


@pytest.mark.parametrize("shape", ["cube", "octahedron", "dodecahedron"])
def test_solvate_shape_never_waters_the_solutes_own_image(shape):
    """Regression: culling the parallelepiped put water on the solute's image.

    The primitive parallelepiped has the same volume as the Wigner-Seitz cell
    but a smaller inscribed diameter, so a solute with radius > 4.45*pad had
    atoms outside the watered region. Because `_overlapWithOther` only checks
    the primary image, water was then placed exactly where the solute's
    periodic image sits. Measured before the fix, at radius 25 and pad 5:
    12 of 400 atoms outside for the octahedron with a closest minimum-image
    solute-water contact of 0.22 A, and 16 of 400 at 0.89 A for the
    dodecahedron. The cube was clean, which is why no rectangular test and no
    line-shaped fixture could see it.
    """
    import numpy as np

    from htmd.builder.solvate import solvate, _cell_vectors, _ws_halfspaces

    solute = _shell_solute()
    pts = solute.coords[:, :, 0]
    center = 0.5 * (pts.min(axis=0) + pts.max(axis=0))
    radius = float(np.linalg.norm(pts - center, axis=1).max())

    smol = solvate(solute, pad=5, shape=shape)
    width = float(smol.box[0, 0])
    assert radius > 4.45 * 5.0, "fixture must exercise the failing regime"

    vectors = _cell_vectors(shape, width)
    normals, offsets = _ws_halfspaces(vectors)

    # Every solute atom must be inside the cell that receives water
    dist = (pts - center) @ normals.T
    outside = int(np.any(dist > offsets + 1e-6, axis=1).sum())
    assert outside == 0, f"{shape}: {outside} solute atoms outside the watered cell"

    # And no water may approach the solute under the minimum image convention
    triples, _ = _water_triples(smol)
    inv = np.linalg.inv(vectors)
    frac = (pts[None, :, :] - triples[:, 0:1, :]) @ inv
    frac -= np.round(frac)
    closest = float(np.linalg.norm(frac @ vectors, axis=2).min())
    assert (
        closest > 1.5
    ), f"{shape}: closest minimum-image solute-water contact {closest:.2f} A"


@pytest.mark.parametrize("shape", ["cube", "octahedron", "dodecahedron"])
def test_solvate_shape_survives_compact_wrapping(shape):
    """moleculekit's triclinic wrapping must accept the cell solvate built.

    wrap(unitcell="compact") re-images whole molecules into the minimum-volume
    representation. If solvate's cell and moleculekit's cell conventions
    disagree, this either raises or leaves water split across a face.
    """
    import numpy as np

    from htmd.builder.solvate import solvate

    smol = solvate(_big_solute(), pad=5, shape=shape)
    smol.wrap(unitcell="compact")

    triples, _ = _water_triples(smol)
    oh1 = np.linalg.norm(triples[:, 1] - triples[:, 0], axis=1)
    oh2 = np.linalg.norm(triples[:, 2] - triples[:, 0], axis=1)
    worst = max(oh1.max(), oh2.max())
    assert worst < 1.5, f"{shape}: water split by wrapping, max O-H = {worst:.2f} A"


def test_solvate_pad_is_per_side():
    """pad stays per-side, so the minimum image distance is 2*pad.

    This matches GROMACS editconf -d. OpenMM's addSolvent padding is the image
    distance itself, which is why openmm.build doubles pad before forwarding.

    Both branches of `width = max(2*radius + 2*pad, 4*pad)` are asserted. The
    `4*pad` floor is OpenMM's rule translated into the per-side convention: it
    stops a solute much smaller than the padding from coming within pad of two
    different periodic copies. A tiny solute triggers the floor; a large one
    is governed by the primary term.
    """
    import numpy as np

    from htmd.builder.solvate import solvate

    def _radius(mol):
        coords = mol.coords[:, :, 0]
        center = 0.5 * (coords.min(axis=0) + coords.max(axis=0))
        return float(np.linalg.norm(coords - center, axis=1).max())

    # Primary term: radius 25 dominates pad 5, so width is 2*25 + 2*5 = 60
    big = _big_solute()
    r_big = _radius(big)
    assert 2 * r_big > 2 * 5.0, "fixture must exercise the primary term"
    smol = solvate(big, pad=5.0, shape="cube")
    assert float(smol.box[0, 0]) == pytest.approx(2 * r_big + 2 * 5.0, abs=1e-3)

    # Floor: radius 0.75 is far smaller than pad 12, so width is 4*12 = 48
    tiny = _solute()
    r_tiny = _radius(tiny)
    assert 2 * r_tiny + 2 * 12.0 < 4 * 12.0, "fixture must exercise the floor"
    smol = solvate(tiny, pad=12.0, shape="cube")
    assert float(smol.box[0, 0]) == pytest.approx(4 * 12.0, abs=1e-3)


@pytest.mark.parametrize("shape", ["octahedron", "dodecahedron"])
def test_solvate_shape_rejects_incompatible_arguments(shape):
    from htmd.builder.solvate import solvate

    with pytest.raises(ValueError, match="single cell width"):
        solvate(_solute(), minmax=[[-20, -20, -20], [20, 20, 20]], shape=shape)
    with pytest.raises(ValueError, match="single cell width"):
        solvate(_solute(), centersel="all", boxsize=[40, 40, 60], shape=shape)
    with pytest.raises(ValueError, match="single cell width"):
        solvate(_solute(), pad=10, negz=5, shape=shape)


@pytest.mark.parametrize("shape", ["cube", "octahedron", "dodecahedron"])
def test_solvate_shape_rejects_nonpositive_width(shape):
    """A zero or negative width must raise, not build a degenerate cell.

    _cell_lengths_and_angles divides by the vector norms, so a zero width
    yields nan angles with only a RuntimeWarning.
    """
    from htmd.builder.solvate import solvate

    with pytest.raises(ValueError, match="positive cell width"):
        solvate(_solute(), centersel="all", boxsize=0, shape=shape)
    with pytest.raises(ValueError, match="positive cell width"):
        solvate(_solute(), centersel="all", boxsize=-10, shape=shape)


@pytest.mark.parametrize("shape", ["cube", "octahedron", "dodecahedron"])
def test_solvate_shape_works_with_centersel(shape):
    """centersel combined with a shape must actually build, not just not-raise.

    This is what the shape block's placement is load-bearing for. The block has
    to run before the pre-existing "centersel and boxsize must both be
    specified together." check, or this combination is unreachable while every
    other test still passes. Nothing else guards that ordering.
    """
    import numpy as np

    from htmd.builder.solvate import solvate

    solute = _big_solute()
    smol = solvate(solute, pad=5, shape=shape, centersel="all")

    assert smol.box is not None and np.all(smol.box > 0)
    assert int(smol.atomselect("water").sum()) > 0
    lengths = smol.box[:, 0]
    assert np.allclose(lengths, lengths[0], rtol=1e-4), lengths


@pytest.mark.parametrize("shape", ["cube", "octahedron", "dodecahedron"])
def test_solvate_shape_needs_pad_or_boxsize(shape):
    from htmd.builder.solvate import solvate

    with pytest.raises(ValueError, match="either pad or boxsize"):
        solvate(_solute(), shape=shape)


def test_solvate_unknown_shape_raises():
    from htmd.builder.solvate import solvate

    with pytest.raises(ValueError, match="banana"):
        solvate(_solute(), pad=8, shape="banana")


def test_solvate_exclude_z_leaves_a_gap():
    """No water residue may have an atom inside the excluded z range."""
    import numpy as np

    from htmd.builder.solvate import solvate

    zlo, zhi = -6.0, 6.0
    smol = solvate(_solute(), pad=14, exclude_z=(zlo, zhi))

    triples, _ = _water_triples(smol)
    z = triples[:, :, 2]
    inside = (z > zlo) & (z < zhi)
    assert not inside.any(), f"{int(inside.any(axis=1).sum())} waters inside the gap"

    # Water still present on both sides of the gap
    assert (z.max(axis=1) <= zlo).sum() > 0, "no water below the gap"
    assert (z.min(axis=1) >= zhi).sum() > 0, "no water above the gap"


@pytest.mark.parametrize(
    "bad",
    [
        (6.0, -6.0),
        (6.0, 6.0),
        (6.0,),
        (-6.0, 0.0, 6.0),
        (float("nan"), 6.0),
        (-6.0, float("nan")),
    ],
)
def test_solvate_exclude_z_validation(bad):
    """Every malformed pair must raise, NaN included.

    NaN is the dangerous case: every comparison against it is False, so an
    ordering-only check would let it through and then silently disable the
    filtering instead of raising. Callers derive these values from an atom
    selection, and an empty selection produces NaN, which would silently
    solvate straight through a membrane.
    """
    from htmd.builder.solvate import solvate

    with pytest.raises(ValueError, match="exclude_z"):
        solvate(_solute(), pad=8, exclude_z=bad)


@pytest.mark.parametrize("shape", ["octahedron", "dodecahedron"])
def test_solvate_exclude_z_applies_to_every_shape(shape):
    """exclude_z must apply to the Wigner-Seitz branch too, as documented.

    The fold-in sits after both branches build their mask. A refactor that
    moved it inside the rectangular branch would regress silently, since the
    docstring promises it works for every shape.
    """
    import numpy as np

    from htmd.builder.solvate import solvate

    zlo, zhi = -8.0, 8.0
    smol = solvate(_big_solute(), pad=5, shape=shape, exclude_z=(zlo, zhi))

    triples, _ = _water_triples(smol)
    z = triples[:, :, 2]
    inside = (z > zlo) & (z < zhi)
    assert (
        not inside.any()
    ), f"{shape}: {int(inside.any(axis=1).sum())} waters inside the excluded slab"
    assert int(smol.atomselect("water").sum()) > 0, "no water placed at all"


def test_solvate_exclude_z_matches_two_call_equivalent():
    """One call with exclude_z must match the two-slab pattern it replaces."""
    import numpy as np

    from htmd.builder.solvate import solvate

    mm = [[-25.0, -25.0, -25.0], [25.0, 25.0, 25.0]]
    zlo, zhi = -6.0, 6.0

    single = solvate(_solute(), minmax=mm, exclude_z=(zlo, zhi))
    n_single = int(single.atomselect("water").sum()) // 3

    upper = solvate(_solute(), minmax=[[-25.0, -25.0, zhi], [25.0, 25.0, 25.0]])
    both = solvate(upper, minmax=[[-25.0, -25.0, -25.0], [25.0, 25.0, zlo]])
    n_two = int(both.atomselect("water").sum()) // 3

    assert n_single == pytest.approx(
        n_two, rel=0.05
    ), f"single call {n_single} waters, two calls {n_two}"


def test_solvate_rectangular_sets_the_cell_it_built():
    """Rectangular now reports its own box instead of leaving it unset.

    amber.build used to re-measure the cell with tleap's setBox "vdw", which
    never matched solvate's region.

    The values are pinned exactly, because the old code already returned a
    zero-filled ``(3, 1)`` box inherited from appending the water box. A test
    that only checked the shape and dtype would have passed before the change.
    `_solute()` spans 1.5 A along x and nothing on y or z, so pad=10 gives
    21.5 by 20 by 20.
    """
    import numpy as np

    from htmd.builder.solvate import solvate

    smol = solvate(_solute(), pad=10)
    assert smol.box is not None and smol.box.shape == (3, 1)
    assert np.allclose(
        smol.box.ravel(), [21.5, 20.0, 20.0], atol=1e-3
    ), smol.box.ravel()
    assert np.allclose(smol.boxangles.ravel(), [90.0, 90.0, 90.0])
    for key, value in zip(
        ("a", "b", "c", "alpha", "beta", "gamma"), [21.5, 20.0, 20.0, 90.0, 90.0, 90.0]
    ):
        assert smol.crystalinfo[key] == pytest.approx(value, abs=1e-3)
