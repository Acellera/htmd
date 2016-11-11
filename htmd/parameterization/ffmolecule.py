# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.molecule.molecule import Molecule
from htmd.molecule.vmdparser import guessAnglesAndDihedrals
from htmd.parameterization.detectsoftdihedrals import detectSoftDihedrals
from htmd.parameterization.detectequivalents import detectEquivalents
from htmd.parameterization.fftype import FFTypeMethod, FFType
from htmd.parameterization.ff import RTF, PRM
from htmd.parameterization.ffevaluate import FFEvaluate
from htmd.qm.qmcalculation import *
from htmd.parameterization.phi import setPhi, getPhi
from htmd.progress.progress import ProgressBar
import re
import math
import os
import scipy.optimize as optimize
import numpy as np
from copy import deepcopy


class QMFittingSet:
    # phi_coords  = []
    def __init__(self):
        self.phi = []
        self.qm = []
        self.mm_original = []
        self.mm_zeroed = []
        self.mm_delta = []
        self.mm_fitted = []
        self.coords = []


class FFMolecule(Molecule):
    def __init__(self, filename=None, name=None, rtf=None, prm=None, netcharge=None, method=FFTypeMethod.CGenFF_2b6,
                 basis=BasisSet._6_31G_star, solvent=True, theory=Theory.B3LYP, execution=Execution.Inline, qmcode=Code.PSI4, outdir="./"):
        # filename -- a mol2 format input geometry
        # rtf, prm -- rtf, prm files
        # method  -- if rtf, prm == None, guess atom types according to this method ( of enum FFTypeMethod )
        self.basis = basis
        self.theory = theory
        self.solvent= solvent

        self.solvent_name = "vacuum"
        if solvent: self.solvent_name="water"
      
        if theory == Theory.RHF:   self.theory_name="rhf" 
        if theory == Theory.B3LYP: self.theory_name="b3lyp" 

        if basis == BasisSet._6_31G_star:
            self.basis_name = "6-31g-star"
        elif basis == BasisSet._cc_pVDZ:
            self.basis_name = "cc-pVDZ"
        else:
            raise ValueError("Unknown Basis Set")

        self.execution = execution
        self.qmcode = qmcode
        self.method = method
        self.outdir = outdir

        if not (filename.endswith(".mol2")):
            raise ValueError("Input file must be mol2 format")

        super().__init__(filename=filename, name=name)
        (a, b) = guessAnglesAndDihedrals(self.bonds)
        self.natoms = self.serial.shape[0]
        self.angles = a
        self.dihedrals = b
        ee = detectEquivalents(self)
        self._soft_dihedrals = detectSoftDihedrals(self, ee)
        self._equivalent_atom_groups = ee[0]  # list of groups of equivalent atoms
        self._equivalent_atoms = ee[1]  # list of equivalent atoms, indexed by atom
        self._equivalent_group_by_atom = ee[2]  # mapping from atom index to equivalent atom group
        if netcharge is None:
            self.netcharge = int(round(np.sum(self.charge)))
        else:
            self.netcharge = int(round(netcharge))

        # Canonicalise the atom naming.
        self._rename_mol()

        if rtf and prm:
            # If the user has specified explicit RTF and PRM files go ahead and load those
            self._rtf = RTF(rtf)
            self._prm = PRM(prm)
        else:
            # Otherwise make atom types using the specified method
            # (Right now only MATCH)
            fftype = FFType(self, method=self.method)
            self._rtf = fftype._rtf
            self._prm = fftype._prm
        if not self._rtf or not self._prm:
            raise ValueError("RTF and PRM not defined")

        self.report()

    def report(self):
        print("Net Charge: {}".format(self.netcharge))
        print("Equivalent atom groups:")
        for i in self._equivalent_atom_groups:
            for j in i:
                print(" {}".format(self.name[j]), end="")
            print("")

        print("Soft dihedrals:")
        for i in self._soft_dihedrals:
            for j in i.atoms:
                print(" {}".format(self.name[j]), end="")
            print("")

    def _rename_mol(self):
        # This fixes up the atom naming and reside name to be consistent
        # NB this scheme matches what MATCH does. Don't change it
        # Or the naming will be inconsistent with the RTF
        import re

        hh = dict()

        for i in range(len(self.name)):
            # Remove any character that isn't alpha
            t = re.sub('[^A-Z]*', "", self.name[i].upper())
            # print("RENAMED %s to %s" %(self.name[i], t ) )
            idx = 0

            if not t in hh:
                hh[t] = idx

            idx = hh[t] + 1
            hh[t] = idx

            t += str(idx)
            self.name[i] = t
            self.resname[i] = "MOL"

    def output_directory_name( self ):
        return self.theory_name + "-" + self.basis_name + "-" + self.solvent_name 

    def minimize(self):
        mindir = os.path.join(self.outdir, "minimize", self.output_directory_name()) #self.theory_name + "-" + self.basis_name + "-" + self.solvent_name )
        try:
            os.makedirs(mindir, exist_ok=True)
        except:
            raise OSError('Directory {} could not be created. Check if you have permissions.'.format(mindir))

        # Kick off a QM calculation -- unconstrained geometry optimization
        qm = QMCalculation(self, charge=self.netcharge, optimize=True,
                           directory=mindir, basis=self.basis, theory=self.theory, solvent=self.solvent,
                           execution=self.execution, code=self.qmcode)
        results = qm.results()
        if results[0].errored:
            raise RuntimeError("QM Optimization failed")
        # Replace coordinates with the minimized set
        self.coords = np.atleast_3d(results[0].coords)

    def _fitCharges_map_back_to_charges(self, x):
        charges = np.zeros((self.natoms))

        qsum = 0.
        for i in range(len(x)):
            charges[self._equivalent_atom_groups[i]] = x[i]
            qsum += x[i] * len(self._equivalent_atom_groups[i])
            #  diff = self.netcharge - qsum;
            #  diff = diff / len(self._equivalent_atom_groups[ len(x) ])
            #  print( self._equivalent_atom_groups[ len(x) ] )
            #  print( diff )
            #  charges[ self._equivalent_atom_groups[ len(x) ] ] = diff
        return charges

    def _fitCharges_con(self, x):
        charges = self._fitCharges_map_back_to_charges(x)
        s = np.sum(charges) - self.netcharge
        return s

    def _fitCharges_objective(self, x):
        # Map the fit variables back to per-atom charges
        chisq = 0.
        charges = self._fitCharges_map_back_to_charges(x)

        #    if( range_penalty == 1 ): chisq = 1000.

        for i in range(self._fitCharges_grid.shape[0]):
            ee = np.sum(charges * self._fitCharges_distances[i, :])
            delta_ee = self._fitCharges_esp[i] - ee
            chisq = chisq + (delta_ee * delta_ee)

        return chisq

    def _fitDihedral_objective(self, x):
        inv = math.pi / 180.

        # evaluate the torsion with the input params
        # for each of the phi's poses
        chisq = 0.
        for t in range(self._fitDihedral_results.N):
            e = .0  # FFEvaluate.evaluateTorsion( self._fitDihedral_results["phi_coords"][t], phi )
            for s in range(len(self._fitDihedral_results.phis[t])):
                for j in range(6):
                    e += x[j] * (1. + cos((j + 1) * (self._fitDihedral_results.phis[t][s] * inv) - x[6 + j] * inv))

            e = e + x[12]
            diff = self._fitDihedral_results.mm_delta[t] - e
            chisq += diff * diff

        return chisq

    def _removeCOM(self):
        # Relocate centre of mass to the origin
        for f in range(self.coords.shape[2]):
            com = np.zeros(3)
            mass = 0.
            for i in range(self.coords.shape[0]):
                m = self._rtf.mass_by_type[self._rtf.type_by_index[i]]
                mass = mass + m
                com = com + self.coords[i, :, f] * m
            com = com / mass
            self.coords[:, :, f] = self.coords[:, :, f] - com

    def _try_load_pointfile(self):
        # Load a point file if one exists from a previous job
        pointfile = os.path.join(self.outdir, "esp", self.output_directory_name()   , "00000", "grid.dat")
        if os.path.exists(pointfile):
            f = open(pointfile, "r")
            fl = f.readlines()
            f.close()
            ret = np.zeros((len(fl), 3))
            for i in range(len(fl)):
                s = fl[i].split()
                ret[i, 0] = float(s[0])
                ret[i, 1] = float(s[1])
                ret[i, 2] = float(s[2])
            print("Reusing previously-generated point cloud")
            return ret
        return True

    def fitCharges(self):
        # Remove the COM from the coords, or PSI4 does it and then the grid is incorrectly centred
        self._removeCOM()
        # Kick off a QM calculation -- unconstrained single point with grid
        points = self._try_load_pointfile()
        espdir = os.path.join(self.outdir, "esp", self.output_directory_name() )
        try:
            os.makedirs(espdir, exist_ok=True)
        except:
            raise OSError('Directory {} could not be created. Check if you have permissions.'.format(espdir))

        qmcode = self.qmcode
        if self.qmcode == Code.TeraChem: 
           print("Charge-fitting requires a feature TeraChem doesn't have yet. Using PSI4 instead")
           qmcode = Code.PSI4

        qm = QMCalculation(self, charge=self.netcharge, optimize=False, esp=points, theory=self.theory, solvent=self.solvent,
                           directory=espdir, basis=self.basis, execution=self.execution,
                           code=qmcode)
        results = qm.results()
        if results[0].errored:
            raise RuntimeError("QM Calculation failed")
        esp_grid = results[0].esp_points
        esp = results[0].esp_scalar
        self.coords = results[0].coords

        #    print(results[0].dipole )
        #    print(results[0].quadrupole )
        #    print(results[0].mulliken )

        self._fitCharges_grid = esp_grid
        self._fitCharges_esp = esp

        # set up the restraints to fit

        N = len(self._equivalent_atom_groups)  # - 1
        lb = np.ones((N)) * -1.25
        ub = np.ones((N)) * +1.25
        # If the restraint relates to an H, set the lower bound to 0
        for i in range(N):
            if "H" == self.element[self._equivalent_atom_groups[i][0]]:
                lb[i] = 0.001

        bounds = []
        for a in range(len(lb)):
            bounds.append((lb[a], ub[a]))
            # Start off by equally distributing the mol's charge
        start = np.zeros(N)

        # Precompute the 1/r distances
        self._fitCharges_distances = np.zeros((self._fitCharges_grid.shape[0], self.coords.shape[0]))

        for i in range(self._fitCharges_grid.shape[0]):
            p1 = self._fitCharges_grid[i, :]
            for j in range(self.coords.shape[0]):
                p2 = self.coords[j, :, 0]
                r = np.linalg.norm(p1 - p2)
                self._fitCharges_distances[i, j] = 1. / r
            #    initial_chisq = self._fitCharges_objective( start )

        xopt = optimize.minimize(self._fitCharges_objective, start, method="SLSQP", bounds=bounds,
                                 options={"disp": False}, constraints={'type': 'eq', 'fun': self._fitCharges_con})
        #    xopt = optimize.minimize( self._fitCharges_objective, start, method="L-BFGS-B",
        # bounds = bounds, options={"disp":False} )

        charges = self._fitCharges_map_back_to_charges(xopt.x)

        # Calculate the dipole from the fitted charges
        dpx = dpy = dpz = 0.
        nc = 0.
        for i in range(len(charges)):
            dpx = dpx + charges[i] * self.coords[i, 0, 0]
            dpy = dpy + charges[i] * self.coords[i, 1, 0]
            dpz = dpz + charges[i] * self.coords[i, 2, 0]
            nc = nc + charges[i]
        fac = (2.541766 / 0.529177249)
        dpx *= fac
        dpy *= fac
        dpz *= fac
        dp = math.sqrt(dpx * dpx + dpy * dpy + dpz * dpz)

        fit_chisq = self._fitCharges_objective(xopt.x)

        self.charges = charges
        self._rtf.updateCharges(charges)

        return fit_chisq, results[0].dipole, [dpx, dpy, dpz, dp]

    def getSoftDihedrals(self):
        dd = []
        for d in self._soft_dihedrals:
            dd.append(d.atoms.copy())
        return dd

    #  def scanSoftDihedral(self, phi, directory = "dihedral", step=10):
    #    found=False
    #    phi_to_fit = None
    #    frozens=[]
    #    dih_index=0
    #    i=0
    #    for d in self._soft_dihedrals:
    #      if (d.atoms == phi).all():
    #         phi_to_fit = d
    #         dih_index=i
    #         frozens.append(d.atoms)
    #      else:
    #         pass
    #      i=i+1
    #    if not phi_to_fit: raise ValueError( "specified phi is not a recognised soft dihedral" )
    #
    #    atoms = phi_to_fit.atoms
    #    left  = phi_to_fit.left
    #    right = phi_to_fit.right
    #    equivs= phi_to_fit.equivalents
    #
    ##    step  = 10 # degrees
    #    nstep = (int)(360/step)
    #    cset  = np.zeros( ( self.natoms, 3, nstep ) )
    #
    #    i=0
    #    for phi in range( -180, 180, step ):
    #      cset[:,:,i] = setPhi( self.coords[:,:,0], atoms, left, right, phi )
    #      i=i+1
    #
    #    mol        = self.copy()
    #    mol.coords = cset
    #    try:
    #      os.mkdir( directory )
    #    except:
    #      pass
    #    dih_name = "%s-%s-%s-%s" % ( self.name[atoms[0]], self.name[atoms[1]], self.name[atoms[2]], self.name[atoms[3]] )
    #    qmset   = QMCalculation( mol, charge=self.netcharge, directory="%s/%s" % (directory, dih_name), frozen=frozens, optimized=True )
    #    r = qmset.results()
    #    x=0
    #    ret=[]
    #    for phi in range( -180, 180, step ):
    #      r[x].phi = phi
    #      if r[x].errored == False:
    #        ret.append(r[x])
    #      x=x+1
    #    return ret

    def fitSoftDihedral(self, phi, geomopt=True):
        found = False
        phi_to_fit = None
        frozens = []
        dih_index = 0
        i = 0
        bkp_coords = self.coords.copy()

        for d in self._soft_dihedrals:
            if (d.atoms == phi).all():
                phi_to_fit = d
                dih_index = i
                frozens.append(d.atoms)
            else:
                if not geomopt:
                    frozens.append(d.atoms)
            i += 1
        if not phi_to_fit:
            raise ValueError("specified phi is not a recognised soft dihedral")
        self._makeDihedralUnique(phi_to_fit)

        atoms = phi_to_fit.atoms
        left = phi_to_fit.left
        right = phi_to_fit.right
        equivs = phi_to_fit.equivalents

        step = 10  # degrees
        nstep = int(360 / step)
        cset = np.zeros((self.natoms, 3, nstep))

        i = 0
        for phi in range(-180, 180, step):
            cset[:, :, i] = setPhi(self.coords[:, :, 0], atoms, left, right, phi)
            i += 1

        mol = self.copy()
        mol.coords = cset

        dirname = "dihedral-single-point"
        if geomopt:
            dirname = "dihedral-opt"

        dih_name = "%s-%s-%s-%s" % (self.name[atoms[0]], self.name[atoms[1]], self.name[atoms[2]], self.name[atoms[3]])

        fitdir = os.path.join(self.outdir, dirname, dih_name, self.output_directory_name()) 

        try:
            os.makedirs(fitdir, exist_ok=True)
        except:
            raise OSError('Directory {} could not be created. Check if you have permissions.'.format(fitdir))

        qmset = QMCalculation(mol, charge=self.netcharge, directory=fitdir, frozen=frozens, optimize=geomopt, theory=self.theory, solvent=self.solvent,
                              basis=self.basis, execution=self.execution, code=self.qmcode)

        ret = self._makeDihedralFittingSetFromQMResults(atoms, qmset.results())

        # Get the initial parameters of the dihedral we are going to fit

        param = self._prm.dihedralParam(self._rtf.type_by_index[atoms[0]],
                                        self._rtf.type_by_index[atoms[1]],
                                        self._rtf.type_by_index[atoms[2]],
                                        self._rtf.type_by_index[atoms[3]])

        # Save these parameters as the best fit (fit to beat)
        best_param = np.zeros((13))
        for t in range(6):
            best_param[t] = param[t].k0
            best_param[t + 6] = param[t].phi0
        best_param[12] = 0.

        #    print(param)

        # Evalaute the mm potential with this dihedral zeroed out
        # The objective function will try to fit to the delta between
        # The QM potential and the this modified mm potential

        for t in param:
            t.k0 = t.phi0 = 0.
            t.e14 = 1.  # Always fit with e14 scaling of 1. per CHARMM
        self._prm.updateDihedral(param)

        ffe = FFEvaluate(self)
        #  print(ffe.evaluate( ret.coords[0] ) )
        #  input
        # Now evaluate the ff without the dihedral being fitted
        for t in range(ret.N):
            mm_zeroed = (ffe.evaluate(ret.coords[t])["total"])
            ret.mm_delta.append(ret.qm[t] - mm_zeroed)
            ret.mm_zeroed.append(mm_zeroed)

        mmin1 = min(ret.mm_zeroed)
        mmin2 = min(ret.mm_delta)
        for t in range(ret.N):
            ret.mm_zeroed[t] = ret.mm_zeroed[t] - mmin1
            ret.mm_delta[t] = ret.mm_delta[t] - mmin2

        self._fitDihedral_results = ret
        self._fitDihedral_phi = param

        # Now measure all of the soft dihedrals phis that are mapped to this dihedral
        ret.phis = []
        for i in range(ret.N):
            ret.phis.append([ret.phi[i]])
            for e in equivs:
                ret.phis[i].append(getPhi(ret.coords[i], e))
            #    print ("EQUIVALENT DIHEDRALS FOR THIS DIHEDRAL" )
            #    print(equivs)
            #    print ("PHI VALUES TO FIT")
            #    print (ret.phis)
        # Set up the NOLOPT fit
        #  There are 13 parameters, k,phi for n=1,2,3,4,5,6 and a shift
        N = 13
        # initial guess,
        st = np.zeros(13)
        # bounds

        best_chisq = self._fitDihedral_objective(best_param)
        #    print("CHISQ of initial = %f" % ( best_chisq ) )

        # Now zero out the terms of the dihedral we are going to fit
        bar = ProgressBar(64, description="Fitting")
        for i in range(64):

            (bounds, start) = self._fitDihedral_make_bounds(i)

            xopt = optimize.minimize(self._fitDihedral_objective, start, method="L-BFGS-B", bounds=bounds,
                                     options={'disp': False})

            chisq = self._fitDihedral_objective(xopt.x)
            #      print( "CHISQ of fit = %f " % (chisq) )
            if (chisq < best_chisq):
                best_chisq = chisq
                best_param = xopt.x
            bar.progress()
        bar.stop()
        #    print("Best ChiSQ = %f" %(best_chisq) )

        # Update the target dihedral with the optimized parameters
        # print(param)
        # print(best_param )
        for i in range(6):
            param[i].k0 = best_param[0 + i]
            param[i].phi0 = best_param[6 + i]

        self._prm.updateDihedral(param)
        # print(param)
        param = self._prm.dihedralParam(self._rtf.type_by_index[atoms[0]],
                                        self._rtf.type_by_index[atoms[1]],
                                        self._rtf.type_by_index[atoms[2]],
                                        self._rtf.type_by_index[atoms[3]])
        # print(param)

        # Finally evaluate the fitted potential
        ffe = FFEvaluate(self)
        for t in range(ret.N):
            ret.mm_fitted.append(ffe.evaluate(ret.coords[t])["total"])
        mmin = min(ret.mm_fitted)
        chisq = 0.

        #    print( "QM energies" )
        #    print( ret.qm )

        for t in range(ret.N):
            ret.mm_fitted[t] = ret.mm_fitted[t] - mmin
            delta = ret.mm_fitted[t] - ret.qm[t]
            chisq = chisq + (delta * delta)
        ret.chisq = chisq

        # TODO Score it
        self.coords = bkp_coords

        return ret

    def _fitDihedral_make_bounds(self, i):
        lb = np.zeros(13)
        ub = np.zeros(13)
        start = np.zeros(13)

        bounds = []

        for j in range(6):
            start[j] = 0.
            bounds.append((-20., 20.))

        for j in range(6):
            if i & (2 ** j):
                bounds.append((180., 180.))
                start[6 + j] = 180.
            else:
                bounds.append((0., 0.))
                start[6 + j] = 0.

        bounds.append((-10., 10.))
        return bounds, start

    def _makeDihedralFittingSetFromQMResults(self, atoms, results):
        # Extract the valid QM poses and energies from the QM result set
        # Evaluate the MM on those poses
        ffe = FFEvaluate(self)

        ret = QMFittingSet()
        ret.name = "%s-%s-%s-%s" % (
            self._rtf.names[atoms[0]], self._rtf.names[atoms[1]], self._rtf.names[atoms[2]], self._rtf.names[atoms[3]])

        completed = 0

        qmin = 1.e100
        for q in results:
            if q.completed and not q.errored:
                if q.energy < qmin:
                    qmin = q.energy

        completed = 0
        for q in results:
            if q.completed and not q.errored:
                if (q.energy - qmin) < 20.:  # Only fit against QM points < 20kcal above the minimum
                    completed += 1
                    ret.phi.append(getPhi(q.coords, atoms))

                    ret.qm.append(q.energy - qmin)
                    ret.mm_original.append(ffe.evaluate(q.coords)['total'])
                    ret.coords.append(q.coords)
                    if ret.phi[0] > 175.:
                        ret.phi[0] -= 360.
        mmin = min(ret.mm_original)
        # roughly align the qm with the mm
        for q in range(completed):
            ret.mm_original[q] = ret.mm_original[q] - mmin
        ret.N = completed

        if completed < 5:
            raise RuntimeError("Fewer than 5 valid QM points. Not enough to fit!")

        return ret

    def _makeDihedralUnique(self, phi_to_fit):
        #    (number_of_uses, uses) = self._countUsesOfDihedral( phi_to_fit.atoms )
        #    if( number_of_uses > 1 ):
        # Create a new type for (arbitrarily) a middle atom of the dihedral
        # So that the dihedral we are going to modify is unique
        # TODO -- check symmetry
        #    print( "Dihedral term is not unique. Copying type.." ) # Used %d times, by:" % ( number_of_uses ) )
        # print( uses )

        # Duplicate the dihedrals types so this modified term is unique
#        print("Duplicating types..")
        for i in range(4):
            if not ("x" in self._rtf.type_by_index[phi_to_fit.atoms[i]]):
                self._duplicateTypeOfAtom(phi_to_fit.atoms[i])

        (number_of_uses, uses) = self._countUsesOfDihedral(phi_to_fit.atoms)
        if number_of_uses > 1:
            print(phi_to_fit.atoms)
            print(number_of_uses)
            print(uses)
            raise ValueError("Dihedral term still not unique after duplication")

    def _countUsesOfDihedral(self, aidx):

        # Return the number of uses of the dihedral
        # specified by the types of the 4 atom indices in the aidx list
        #

        #    print( "countUsesOfDihedral in " )
        t1 = self._rtf.type_by_index[aidx[0]]
        t2 = self._rtf.type_by_index[aidx[1]]
        t3 = self._rtf.type_by_index[aidx[2]]
        t4 = self._rtf.type_by_index[aidx[3]]

        count = 0
        uses = []
        for d in self.dihedrals:
            s1 = self._rtf.type_by_index[d[0]]
            s2 = self._rtf.type_by_index[d[1]]
            s3 = self._rtf.type_by_index[d[2]]
            s4 = self._rtf.type_by_index[d[3]]
            if s1 == t1 and s2 == t2 and s3 == t3 and s4 == t4:
                count += 1
                uses.append(d)
            elif s1 == t4 and s2 == t3 and s3 == t2 and s4 == t1:
                count += 1
                uses.append(d)
            #    return(count, uses )
            #    print(uses)

        # Now for each of the uses, remove any which are equivalent
        c = 1
        unique_uses = [aidx]
        g1 = self._equivalent_group_by_atom[aidx[0]]
        g2 = self._equivalent_group_by_atom[aidx[1]]
        g3 = self._equivalent_group_by_atom[aidx[2]]
        g4 = self._equivalent_group_by_atom[aidx[3]]
        for u in uses:
            h1 = self._equivalent_group_by_atom[u[0]]
            h2 = self._equivalent_group_by_atom[u[1]]
            h3 = self._equivalent_group_by_atom[u[2]]
            h4 = self._equivalent_group_by_atom[u[3]]
            equiv = False
            if g1 == h1 and g2 == h2 and g3 == h3 and g4 == h4:
                equiv = True
            if g1 == h4 and g2 == h3 and g3 == h2 and g4 == h1:
                equiv = True
            if equiv is False:
                c += 1
                unique_uses.append(u)
            else:
#                print(" Dih %s-%s-%s-%s and %s-%s-%s-%s are equivalent " % (
#                    self._rtf.names[aidx[0]], self._rtf.names[aidx[1]], self._rtf.names[aidx[2]],
#                    self._rtf.names[aidx[3]], self._rtf.names[u[0]], self._rtf.names[u[1]], self._rtf.names[u[2]],
#                    self._rtf.names[u[3]]))
                pass
                #  return(count, uses )
            #    print( c )
            #    print( unique_uses )
        return c, unique_uses

    def _duplicateTypeOfAtom(self, aidx):
        # This duplicates the type of the specified atom
        # First get the type
        atype = self._rtf.type_by_index[aidx]

        # perhaps the type is already a duplicate? if so
        # remove the duplicated suffix
        atype = re.sub("x[0123456789]+$", "", atype)
        i = 0
        # make the new type name
        while ("%sx%d" % (atype, i)) in self._rtf.types:
            i += 1

        newtype = "%sx%d" % (atype, i)
        print("Creating new type %s from %s for atom %s" % (newtype, atype, self._rtf.names[aidx]))

        # duplicate the type in the fields RTF -- todo: move to a method in the RTF
        self._rtf.type_by_index[aidx] = newtype
        self._rtf.mass_by_type[newtype] = self._rtf.mass_by_type[atype]
        self._rtf.types.append(newtype)
        self._rtf.type_by_name[self._rtf.names[aidx]] = newtype
        self._rtf.type_by_index[aidx] = newtype
        self._rtf.typeindex_by_type[newtype] = self._rtf.typeindex_by_type[atype] + 1000
        self._rtf.element_by_type[newtype] = self._rtf.element_by_type[atype]

        #    # Now also reset the type of  any atoms that share equivalency
        for bidx in self._equivalent_atoms[aidx]:
            if aidx != bidx:
                if "x" in self._rtf.type_by_index[bidx]:
                    raise RuntimeError(
                        "Equivalent atom already has a duplicated type: {} {}".format(bidx, self._rtf.type_by_index[bidx]))
                self._rtf.type_by_index[bidx] = newtype
                self._rtf.type_by_name[self._rtf.names[bidx]] = newtype

        # the PRM parameters will be automatically duplicated by forcing an ff evaluation
        ffe = FFEvaluate(self)
        ffe.evaluate(self.coords)

    def plotDihedralFit(self, fit, show=True):
        import matplotlib as mpl
        if not show:
            mpl.use('Agg')
        import matplotlib.pyplot as plt

        fh = plt.figure()
        ax1 = fh.gca()
        ax1.set_xlim(-180., 180.)
        ax1.set_xticks([-180, -135, -90, -45, 0, 45, 90, 135, 180])
        ax1.set_xlabel("Phi")
        ax1.set_ylabel("kcal/mol")
        ax1.set_title(fit.name)

        x = sorted(fit.phi)
        plotdata=[]
        for i in range(len(fit.phi) ):
           plotdata.append( (fit.phi[i], fit.qm[i] ) )
        plotdata=sorted(plotdata)
        plotdatax = [float(i[0]) for i in plotdata]
        plotdatay = [float(i[1]) for i in plotdata]
        ax1.plot(plotdatax, plotdatay , label="QM", color="r", marker="o")

        plotdata=[]
        for i in range(len(fit.phi) ):
           plotdata.append( (fit.phi[i], fit.mm_original[i] ) )
        plotdata=sorted(plotdata)
        plotdatax = [float(i[0]) for i in plotdata]
        plotdatay = [float(i[1]) for i in plotdata]
        ax1.plot(plotdatax, plotdatay, label="MM Original", color="b", marker="d")

        plotdata=[]
        for i in range(len(fit.phi) ):
           plotdata.append( (fit.phi[i], fit.mm_fitted[i] ) )
        plotdata=sorted(plotdata)
        plotdatax = [float(i[0]) for i in plotdata]
        plotdatay = [float(i[1]) for i in plotdata]
        ax1.plot(plotdatax, plotdatay, label="MM Fitted", color="g", marker="s")



        #ax1.plot(fit.phi, fit.qm, label="QM", color="r", marker="o")
        #ax1.plot(fit.phi, fit.mm_original, label="MM Original", color="b", marker="d")
        #ax1.plot(fit.phi, fit.mm_fitted, label="MM Fitted", color="g", marker="s")
        ##    ax1.plot( fit.phi , fit.mm_zeroed  , label="MM With phi zeroed", color="black", marker="x" )
        ##    ax1.plot( fit.phi , fit.mm_delta   , label="QM-MM target", color="magenta", marker="x" )
        ax1.legend(prop={'size': 8})
        if show:
            plt.show()
        else:
            plotdir = os.path.join(self.outdir, "parameters", self.method.name, self.output_directory_name(), "plots")
            try:
                os.makedirs(plotdir, exist_ok=True)
            except:
                raise OSError('Directory {} could not be created. Check if you have permissions.'.format(plotdir))
            tf = os.path.join(plotdir, fit.name) + ".svg"
            plt.savefig(tf, format="svg")
            return tf

    def write(self, filename, sel=None, type=None, typemap=None):
        if hasattr(self, "_rtf"):  # Update base Molecule's attributes so write() works correctly
            for i in range(self.charge.shape[0]):
                self.segid[i] = "L"
                self.charge[i] = self._rtf.charge_by_name[self.name[i]]
                self.atomtype[i] = self._rtf.type_by_name[self.name[i]]

        ref_atomtype = deepcopy(self.atomtype)
        if typemap:
            for i in range(self.charge.shape[0]):
                self.atomtype[i] = typemap[self.atomtype[i]]

        super().write(filename, sel=sel, type=type)

        self.atomtype = ref_atomtype
