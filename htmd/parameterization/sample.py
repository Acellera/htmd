# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
from htmd.molecule.molecule import Molecule
import htmd.builder.charmm as charmm
import numpy.linalg as npla
import numpy as np
from optparse import OptionParser
import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pylab as plt
from htmd.molecule.util import maxDistance
from htmd.parameterization.ffmolecule import FFMolecule
from htmd.protocols.equilibration_v2 import Equilibration


class Sample:
    def __init__(self, mol, rtf, prm, outdir, solvated):
        out = os.path.join(outdir, "equil")
        traj = os.path.join(out, "output.xtc")
        ftraj = os.path.join(out, "filtered.xtc")
        pdb = os.path.join(out, "structure.pdb")

        if not os.access(traj, os.R_OK):
            print("Preparing a MD run of 50 ns:")
            self._prep_and_run(mol, rtf, prm, outdir, solvated)

        print("Analysing trajectory")
        self._analyse(mol, pdb, rtf, prm, traj, ftraj)

    def _analyse(self, mol, pdb, rtf, prm, traj, ftraj):
        t = Molecule(pdb)
        t.read(traj)
        t.filter('not water')
        t.write(ftraj)
        m = FFMolecule(filename=mol, rtf=rtf, prm=prm)
        m.read(ftraj)
        torsions = m.getParameterizableDihedrals()
        # For each torsion
        for i in range(len(torsions)):
            # Create title
            title = '{}-{}-{}-{}'.format(m.name[torsions[i][0]], m.name[torsions[i][1]], m.name[torsions[i][2]],
                                         m.name[torsions[i][3]])
            # Measure
            (r, theta) = self._measure_torsion(torsions[i], m.coords)

            self._plot_scatter(r, theta, title)
            self._plot_hist(theta, title)

    def _measure_torsion(self, aidx, coords):
        r = []
        theta = []
        # nprint(aidx)
        for i in range(coords.shape[2]):
            phi = self._measure_phi(aidx, coords, i)
            theta.append(phi)
            r.append(i / 10)  # in ns # shouldn't really hard code
            #  print(coords.shape)
            #  print(r)
        #      print(theta)
        return r, theta

    def _measure_phi(self, aidx, coords, frame):
        pos1 = coords[aidx[0] - 1, :, frame]
        pos2 = coords[aidx[1] - 1, :, frame]
        pos3 = coords[aidx[2] - 1, :, frame]
        pos4 = coords[aidx[3] - 1, :, frame]

        r12 = pos1 - pos2
        r23 = pos2 - pos3
        r34 = pos3 - pos4

        a = np.cross(r12, r23)
        b = np.cross(r23, r34)
        c = np.cross(r23, a)

        ra = 1. / npla.norm(a)
        rb = 1. / npla.norm(b)
        rc = 1. / npla.norm(c)

        b *= rb

        cos_phi = (np.dot(a, b)) * ra
        sin_phi = (np.dot(c, b)) * rc

        phi = - np.arctan2(sin_phi, cos_phi)
        return phi
        # phi = phi * 180./ np.pi;
        # print(phi)
        # return phi

    def _plot_scatter(self, r, theta, title):
        import matplotlib as mpl
        mpl.use('Agg')
        area = np.ones(len(r)) * 10
        colors = theta
        ax = plt.subplot(111, polar=True)
        c = plt.scatter(theta, r, c=colors, s=area, cmap=plt.cm.get_cmap('hsv'))
        c.set_alpha(0.75)
        ax.set_title(title)
        plt.savefig("torsion-polar-" + title + ".svg")

    def _plot_hist(self, theta, title):
        import matplotlib as mpl
        mpl.use('Agg')
        for i in range(len(theta)):
            theta[i] *= 180 / np.pi

        # area = np.ones(len(r)) * 10
        ax = plt.subplot(111)
        plt.hist(theta, 90, normed=1, alpha=0.75, range=[-180, 180])
        plt.xlim(-180, 180)
        plt.xticks([-180, -135, -90, -45, 0, 45, 90, 135, 180])
        ax.set_title(title)
        plt.savefig("torsion-hist-" + title + ".svg")

    def _prep_and_run(self, mol, rtf, prm, outdir, solvated):
        from htmd.builder.solvate import solvate
        from htmd.apps.acemdlocal import AcemdLocal
        # Do a simple solvation then run for 50ns
        ionize = True

        mol = Molecule(mol)
        mol.center()
        mol.set("segid", "L")
        d = maxDistance(mol, 'all') + 6

        if solvated:
            mol = solvate(mol, minmax=[[-d, -d, -d], [d, d, d]])

        if not solvated:
            ionize = False

        build_dir = os.path.join(outdir, "build")
        equil_dir = os.path.join(outdir, "equil")

        rtfs = ['top/top_water_ions.rtf', rtf]
        prms = ['par/par_water_ions.prm', prm]
        charmm.build(mol, topo=rtfs, param=prms, outdir=build_dir, ionize=ionize)
        md = Equilibration()
        md.runtime = 50
        md.timeunits = 'ns'
        md.temperature = 300
        md.write(build_dir, equil_dir)
        mdx = AcemdLocal()
        mdx.submit(equil_dir)
        mdx.wait()


def sample_main():
    parser = OptionParser()
    parser.add_option("-m", "--mol", action="store", dest="mol", default="mol.pdb",
                      help="Molecule in mol2 or PDB format")
    parser.add_option("-r", "--rtf", action="store", dest="rtf", default="mol.rtf", help="RTF parameter file")
    parser.add_option("-p", "--prm", action="store", dest="prm", default="mol.prm", help="PRM parameter file")
    parser.add_option("-o", "--out", action="store", dest="out", default="./", help="Output directory")
    #   parser.add_option( "-s", action="store_true", dest="solvated", default=False, help="Solvate model" )
    # TODO: the non-solvated run needs more work than just changing solvation (vacuum simulation MD parameters)

    options, args = parser.parse_args()
    testfile(options.mol)
    testfile(options.rtf)
    testfile(options.prm)

    Sample(mol=options.mol, rtf=options.rtf, prm=options.prm, outdir=options.out, solvated=True)


def testfile(filename):
    if not os.access(filename, os.R_OK):
        print("Cannot access file [" + filename + "]")
        sys.exit(1)


if __name__ == "__main__":
    # sample_main()
    sys.exit(0)
