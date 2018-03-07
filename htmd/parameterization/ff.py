# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

import math
import re
from copy import deepcopy
from htmd.version import version as htmdversion


class BondPrm:
    def __init__(self, types, r0=0., k0=0.):
        self.types = types
        self.r0 = r0
        self.k0 = k0


class AnglePrm:
    def __init__(self, types, theta0=0., k0=0., rUB=0., kUB=0.):
        self.types = types
        self.theta0 = theta0
        self.k0 = k0
        self.rUB = rUB
        self.kUB = kUB


class TorsPrm:
    def __init__(self, types, n=0, k0=0., phi0=0., e14=1., improper=False):
        self.types = types
        self.n = n
        self.k0 = k0
        self.phi0 = phi0
        self.e14 = e14
        self.improper = improper


class NBPrm:
    def __init__(self, types, emin=0., rmin=0., emin_14=None, rmin_14=None):
        self.types = types
        self.emin = emin
        self.rmin = rmin
        self.emin_14 = emin_14
        self.rmin_14 = rmin_14


class PRM:
    def __init__(self, filename):
        f = open(filename, "r")
        lines = f.readlines()
        f.close()

        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.impropers = []
        self.nonbonded = []
        mode = None
        skip = 0
        for l in lines:
            l = l.strip()
            l = re.sub("!.*$", "", l)  # Remove any comment

            if l == "BONDS":
                mode = "BONDS"
                skip = 1
            elif l == "ANGLES":
                mode = "ANGLES"
                skip = 1
            elif l == "DIHEDRALS":
                mode = "DIHEDRALS"
                skip = 1
            elif l == "IMPROPER":
                mode = "IMPROPER"
                skip = 1
            elif l.startswith("NONBONDED "):
                mode = "NONBONDED"
                skip = 2
            elif l is "":
                mode = None
            if skip > 0:
                skip -= 1
            elif mode:
                ll = l.split()
                if mode.startswith("BONDS"):
                    self.bonds.append(BondPrm([ll[0], ll[1]], r0=float(ll[3]), k0=float(ll[2])))
                elif mode.startswith("ANGLES"):
                    if len(ll) <= 5:
                        self.angles.append(AnglePrm([ll[0], ll[1], ll[2]], theta0=float(ll[4]), k0=float(ll[3])))
                    else:
                        self.angles.append(
                            AnglePrm([ll[0], ll[1], ll[2]], theta0=float(ll[4]), k0=float(ll[3]), rUB=float(ll[6]),
                                     kUB=float(ll[5])))
                elif mode.startswith("DIHEDRALS"):
                    self.dihedrals.append(
                        TorsPrm([ll[0], ll[1], ll[2], ll[3]], n=int(ll[5]), k0=float(ll[4]), phi0=float(ll[6])))
                elif mode.startswith("IMPROPER"):
                    self.impropers.append(
                        TorsPrm([ll[0], ll[1], ll[2], ll[3]], n=int(ll[5]), k0=float(ll[4]), phi0=float(ll[6])))
                elif mode.startswith("NONBONDED"):
                    # Charm prm stores  rmin/2 for some pointless reason. remember to multiply by 2
                    if len(ll) < 5:
                        self.nonbonded.append(NBPrm([ll[0]], emin=float(ll[2]), rmin=float(ll[3]) * 2.))
                    else:
                        self.nonbonded.append(
                            NBPrm([ll[0]], emin=float(ll[2]), rmin=float(ll[3]) * 2., emin_14=float(ll[5]),
                                  rmin_14=float(ll[6])) * 2.)

    def write(self, filename):
        for i in self.dihedrals:
            if i.e14 != 1.0:
                raise ValueError("Can't express 1-4 electrostatic scaling in Charmm file format")

        f = open(filename, "w")
        print("* prm file built by HTMD parameterize version {}".format(htmdversion()), file=f)
        print("*\n", file=f)
        print("BONDS", file=f)
        for a in self.bonds:
            print("%-6s %-6s %8.2f %8.4f" % (a.types[0], a.types[1], a.k0, a.r0), file=f)
        print("\nANGLES", file=f)
        for a in self.angles:
            if a.kUB > 0.:
                print("%-6s %-6s %-6s %8.2f %8.2f %8.2f %8.2f" % (
                    a.types[0], a.types[1], a.types[2], a.k0, a.theta0, a.kUB, a.rUB), file=f)
            else:
                print("%-6s %-6s %-6s %8.2f %8.2f" % (a.types[0], a.types[1], a.types[2], a.k0, a.theta0), file=f)
        print("\nDIHEDRALS", file=f)
        for a in self.dihedrals:
            print("%-6s %-6s %-6s %-6s %12.8f %d %12.8f" % (
                a.types[0], a.types[1], a.types[2], a.types[3], a.k0, a.n, a.phi0), file=f)
        print("\nIMPROPER", file=f)
        for a in self.impropers:
            print("%-6s %-6s %-6s %-6s %12.8f %d %12.8f" % (
                a.types[0], a.types[1], a.types[2], a.types[3], a.k0, a.n, a.phi0), file=f)
        print("\nNONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -", file=f)
        print("cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5", file=f)
        for a in self.nonbonded:
            if a.emin_14 is not None:
                # Charmm prm stores rmin/2
                print("%-6s 0.0000 %8.4f %8.4f 0.0000 %8.4f %8.4f" % (
                    a.types[0], a.emin, a.rmin / 2., a.emin_14, a.rmin_14 / 2.), file=f)
            else:
                print("%-6s 0.0000 %8.4f %8.4f" % (a.types[0], a.emin, a.rmin / 2.), file=f)
        f.close()

    def writeFrcmod(self, rtf, filename):

        # Make a type mapping for any name != 2 chars in length
        # Because Amber file formats are horrid
        atom_type_map = dict()
        idx = 97
        for atom_type in rtf.types:
            atom_type_alias = atom_type
            if len(atom_type) != 2:
                atom_type_alias = "z%c" % idx
                idx += 1
            atom_type_map[atom_type] = atom_type_alias

        # check to see whether the parameters can be expressed in Amber FRCMOD format
  #      if len(self.impropers) != 0:
  #          raise ValueError("Can't express CHARMM-style impropers in Amber file format")
        for i in self.angles:
            if i.rUB != 0. or i.kUB != 0.:
                raise ValueError("Can't express Urey-Bradley terms in Amber file format")
        for i in self.nonbonded:
            if i.rmin_14 is None or i.emin_14 is None:
                raise ValueError("Can't express unset 1-4 VdW terms in Amber file format")
            eps = math.fabs(i.emin - i.emin_14 * 2.)
            if i.rmin != i.rmin_14 or eps > 1.e-6:
                raise ValueError("Can't express 1-4 VdW terms that aren't 0.5x scaled in Amber file format")

        f = open(filename, "w")
        print("Frcmod generated by HTMD parameterize version {}".format(htmdversion()), file=f)
        print("MASS", file=f)
        for i in rtf.types:
            print("%s %f %f" % (atom_type_map[i], rtf.mass_by_type[i], 0.), file=f)

        print("\nBOND", file=f)
        for i in self.bonds:
            print("%s-%s %f %f" % (atom_type_map[i.types[0]], atom_type_map[i.types[1]], i.k0, i.r0), file=f)

        print("\nANGL", file=f)
        for i in self.angles:
            print("%s-%s-%s %f %f" % (atom_type_map[i.types[0]], atom_type_map[i.types[1]], atom_type_map[i.types[2]], i.k0, i.theta0), file=f)

        print("\nDIHE", file=f)
        output = dict()
        for i in self.dihedrals:
            assert i.improper == False
            name = "%s-%s-%s-%s" % (atom_type_map[i.types[0]], atom_type_map[i.types[1]], atom_type_map[i.types[2]], atom_type_map[i.types[3]])
            if not (name in output):
                output[name] = 1
                prmx = self.dihedralParam(i.types[0], i.types[1], i.types[2], i.types[3])

                # Prune prms that are zero
                prm=list()
                for pi in range(len(prmx)):
                    dih = prmx[pi]
                    if dih.k0 != 0.: prm.append(dih)
                # HACK: leave at least one dihedral, even if the force constant is 0,
                #       otherwise "tleap" is not happy!
                if len(prm) == 0:
                    prm.append(prmx[0])

                for p, dih in enumerate(prm):
                    scee = 1. / dih.e14
                    scnb = 2.
                    if p == len(prm)-1:
                        per = dih.n
                    else:
                        per = -dih.n  # All terms of the same dihedral except the last one should be negative. http://ambermd.org/formats.html#frcmod
                    print("%2s-%2s-%2s-%2s 1 %12.6f %12.6f %12.6f %12.6f %12.6f" %
                          (atom_type_map[i.types[0]], atom_type_map[i.types[1]], atom_type_map[i.types[2]], atom_type_map[i.types[3]], dih.k0, dih.phi0, per, scee, scnb),
                          file=f)

        print("\nIMPR", file=f)
        output = dict()
        for i in self.impropers:
          if(i.improper == True):
            name = "%s-%s-%s-%s" % (atom_type_map[i.types[0]], atom_type_map[i.types[1]], atom_type_map[i.types[2]], atom_type_map[i.types[3]])
            if not (name in output):
                output[name] = 1
                prm = self.improperParam(i.types[0], i.types[1], i.types[2], i.types[3])
                for pi in range(len(prm)):
                    p = prm[pi]
                    sign = 1
                    if p.k0 != 0.:
                       print("%2s-%2s-%2s-%2s     %f %f %f" %
                          (atom_type_map[i.types[0]], atom_type_map[i.types[1]], atom_type_map[i.types[2]], atom_type_map[i.types[3]], p.k0, p.phi0, sign * p.n ), file=f)

        print("\nNONB", file=f)
        # Have to iterate over the types in use, which include cloned types, an map them back
        # to original type (which has the same vdw params), because a copy of a copy won't be in self.nonbonded.
        for atom_type in rtf.types:
            original_atom_type = re.sub('x[0123456789]$', '', atom_type)
            for nonbonded in self.nonbonded:
                if nonbonded.types[0] == original_atom_type:
                    print("%s %f %f" % (atom_type_map[atom_type], 0.5 * nonbonded.rmin, -nonbonded.emin), file=f)

        print("", file=f)
        f.close()

        return atom_type_map

    def vdwParam(self, n1, n2, s14):
        p1 = None
        p2 = None
        for b in self.nonbonded:
            if b.types[0] == n1:
                p1 = b
            if b.types[0] == n2:
                p2 = b
                # not found, maybe it's a duplicate that needs new params
        if (not p1) and ("x" in n1):
            xn1 = re.sub("x[0123456789]$", "", n1)
            for b in self.nonbonded:
                if b.types[0] == xn1:
                    p1 = b
                    b = deepcopy(b)
                    b.types[0] = n1
                    self.nonbonded.append(b)
            if n1 == n2:
                p2 = b

        if (not p2) and ("x" in n2):
            xn2 = re.sub("x[0123456789]$", "", n2)
            for b in self.nonbonded:
                if b.types[0] == xn2:
                    p2 = b
                    b = deepcopy(b)
                    b.types[0] = n2
                    self.nonbonded.append(b)

        if (not p1) or (not p2):
            raise ValueError("Could not find nb parameters for %s - %s" % (n1, n2))

        emin_1 = p1.emin
        emin_2 = p2.emin
        rmin_1 = p1.rmin
        rmin_2 = p2.rmin

        if s14 and (p1.emin_14 is not None):
            emin_1 = p1.emin_14
            rmin_1 = p1.rmin_14

        if s14 and (p2.emin_14 is not None):
            emin_2 = p2.emin_14
            rmin_2 = p2.rmin_14
        neg_emin = math.sqrt(emin_1 * emin_2)
        rmin = 0.5 * (rmin_1 + rmin_2)
        A = neg_emin * math.pow(rmin, 12)
        B = 2.0 * neg_emin * math.pow(rmin, 6)
        #    print( "emin=%f,%f rmin=%f,%f A=%f B=%f" % (emin_1, emin_2, rmin_1, rmin_2, A,B ) )

        return A, B

    def bondParam(self, n1, n2):
        for b in self.bonds:
            if b.types[0] == n1 and b.types[1] == n2:
                return b
            if b.types[1] == n1 and b.types[0] == n2:
                return b

        # not found, maybe it's a duplicate that needs new params
        if ("x" in n1) or ("x" in n2):
            xn1 = re.sub("x[0123456789]$", "", n1)
            xn2 = re.sub("x[0123456789]$", "", n2)
            b = self.bondParam(xn1, xn2)
            c = BondPrm([n1, n2], r0=b.r0, k0=b.k0)
            self.bonds.append(c)
            return b

        raise ValueError("Could not find bond parameters for %s-%s" % (n1, n2))

    def angleParam(self, n1, n2, n3):
        for b in self.angles:
            if b.types[0] == n1 and b.types[1] == n2 and b.types[2] == n3:
                return b
            if b.types[2] == n1 and b.types[1] == n2 and b.types[0] == n3:
                return b

        # not found, maybe it's a duplicate that needs new params
        if ("x" in n1) or ("x" in n2) or ("x" in n3):
            xn1 = re.sub("x[0123456789]$", "", n1)
            xn2 = re.sub("x[0123456789]$", "", n2)
            xn3 = re.sub("x[0123456789]$", "", n3)
            b = self.angleParam(xn1, xn2, xn3)
            c = AnglePrm([n1, n2, n3], theta0=b.theta0, k0=b.k0, rUB=b.rUB, kUB=b.kUB)
            self.angles.append(c)
            return b

        raise ValueError("Could not find angle parameters for %s-%s-%s" % (n1, n2, n3))

    def dihedralParam(self, n1, n2, n3, n4):
        ret = []
        found = False
        for i in range(6):
            ret.append(TorsPrm([n1, n2, n3, n4], phi0=0., n=i + 1, k0=0.))

        for b in self.dihedrals:
            if b.types[0] == n1 and b.types[1] == n2 and b.types[2] == n3 and b.types[3] == n4:
                ret[b.n - 1] = b
                found = True
            elif b.types[3] == n1 and b.types[2] == n2 and b.types[1] == n3 and b.types[0] == n4:
                # print(b)
                ret[b.n - 1] = b
                found = True

        # not found, maybe it's a duplicate that needs new params
        if not found:
            if ("x" in n1) or ("x" in n2) or ("x" in n3) or ("x" in n4):
                xn1 = re.sub("x[0123456789]$", "", n1)
                xn2 = re.sub("x[0123456789]$", "", n2)
                xn3 = re.sub("x[0123456789]$", "", n3)
                xn4 = re.sub("x[0123456789]$", "", n4)
                b = self.dihedralParam(xn1, xn2, xn3, xn4)
                r = []
                # print(b)
                for c in b:
                    c = deepcopy(c)
                    c.types[0] = n1
                    c.types[1] = n2
                    c.types[2] = n3
                    c.types[3] = n4
                    self.dihedrals.append(c)
                    r.append(c)
                return r

        if not found:
            raise ValueError("Could not find dihedral parameters for %s-%s-%s-%s" % (n1, n2, n3, n4))
        return ret

    def improperParam(self, n1, n2, n3, n4):
        ret = []
        for b in self.impropers:
            if b.types[0] == n1 and b.types[1] == n2 and b.types[2] == n3 and b.types[3] == n4:
                ret.append(b)
            elif b.types[3] == n1 and b.types[2] == n2 and b.types[1] == n3 and b.types[0] == n4:
                ret.append(b)

        # not found, maybe it's a duplicate that needs new params
        if len(ret) == 0:
            if ("x" in n1) or ("x" in n2) or ("x" in n3) or ("x" in n4):
                xn1 = re.sub("x[0123456789]$", "", n1)
                xn2 = re.sub("x[0123456789]$", "", n2)
                xn3 = re.sub("x[0123456789]$", "", n3)
                xn4 = re.sub("x[0123456789]$", "", n4)
                b = self.improperParam(xn1, xn2, xn3, xn4)
                r = []
                for c in b:
                    c = deepcopy(c)
                    c.types[0] = n1
                    c.types[1] = n2
                    c.types[2] = n3
                    c.types[3] = n4
                    self.impropers.append(c)
                    r.append(c)
                return r

        if len(ret) == 0:
            raise ValueError("Could not find improper parameters for %s-%s-%s-%s" % (n1, n2, n3, n4))
        return ret

    def updateDihedral(self, phi):
        pp = deepcopy(phi)
        for p in self.dihedrals:
            if phi[0].types[0] != p.types[0] or phi[0].types[1] != p.types[1] or phi[0].types[2] != p.types[2] or \
                            phi[0].types[3] != p.types[3]:
                pp.append(p)
        self.dihedrals = pp


class AmberPRM(PRM):
    def __init__(self, prepi, frcmod):

        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.impropers = []
        self.nonbonded = []

        import math
        f = open(frcmod, "r")
        lines = f.readlines()
        f.close()
        section = None
        for ff in lines:
            ff = ff.strip()
            if not len(ff):
                continue
            if ff == "MASS":
                section = "MASS"
            elif ff == "BOND":
                section = "BOND"
            elif ff == "ANGLE":
                section = "ANGLE"
            elif ff == "DIHE":
                section = "DIHE"
            elif ff == "IMPROPER":
                section = "IMPROPER"
            elif ff == "NONBON":
                section = "NONBON"
            else:
                if section == "BOND":
                    x = ff[5:].split()
                    y = ff[0:5].split("-")
                    self.bonds.append(BondPrm([y[0].strip(), y[1].strip()], r0=float(x[1]), k0=float(x[0])))
                elif section == "ANGLE":
                    x = ff[8:].split()
                    y = ff[0:8].split("-")
                    self.angles.append(
                        AnglePrm([y[0].strip(), y[1].strip(), y[2].strip()], theta0=float(x[1]), k0=float(x[0])))
                    pass
                elif section == "DIHE":
                    x = ff[11:].split()
                    y = ff[0:11].split("-")
                    self.dihedrals.append(TorsPrm([y[0].strip(), y[1].strip(), y[2].strip(), y[3].strip()],
                                                  n=int(math.fabs(int(float(x[3])))), k0=float(x[1]) / float(x[0]),
                                                  phi0=float(x[2]), e14=1. / 1.2))
                elif section == "IMPROPER":

                    x = ff[11:].split()
                    y = ff[0:11].split("-")
                    # Amber impropers have the same potential as dihedrals, except the scaling factor is different
                    self.impropers.append(TorsPrm([y[0].strip(), y[1].strip(), y[2].strip(), y[3].strip()],
                                                  n=int(math.fabs(int(float(x[2])))), k0=float(x[0]), phi0=float(x[1]),
                                                  e14=1. / 1.2, improper=True))
                elif section == "NONBON":
                    x = ff.split()
                    y = x[0].split("-")
                    rmin = float(x[1]) * 2.
                    emin = - float(x[2])
                    # Amber always scales 1-4 VDW interactions by 0.5
                    self.nonbonded.append(NBPrm([y[0].strip()], emin=emin, rmin=rmin, emin_14=0.5 * emin, rmin_14=rmin))
                    pass

    def dihedralParam(self, n1, n2, n3, n4):
        ret = []
        found = False
        for i in range(6):
            ret.append(TorsPrm([n1, n2, n3, n4], phi0=0., n=i + 1, k0=0., e14=1.0 / 1.2))

        for b in self.dihedrals:
            # print(b.types)
            if b.types[0] == n1 and b.types[1] == n2 and b.types[2] == n3 and b.types[3] == n4:
                ret[b.n - 1] = b
                found = True
            elif b.types[3] == n1 and b.types[2] == n2 and b.types[1] == n3 and b.types[0] == n4:
                # print(b)
                ret[b.n - 1] = b
                found = True

        # not found, maybe it's a duplicate that needs new params
        if not found:
            if ("x" in n1) or ("x" in n2) or ("x" in n3) or ("x" in n4):
                xn1 = re.sub("x[0123456789]$", "", n1)
                xn2 = re.sub("x[0123456789]$", "", n2)
                xn3 = re.sub("x[0123456789]$", "", n3)
                xn4 = re.sub("x[0123456789]$", "", n4)
                b = self.dihedralParam(xn1, xn2, xn3, xn4)
                r = []
                for c in b:
                    c = deepcopy(c)
                    c.types[0] = n1
                    c.types[1] = n2
                    c.types[2] = n3
                    c.types[3] = n4
                    self.dihedrals.append(c)
                    r.append(c)
                return r

        if not found:
            raise ValueError("Could not find dihedral parameters for %s-%s-%s-%s" % (n1, n2, n3, n4))
        return ret

    def improperParam(self, n1, n2, n3, n4):

        ret = []
        for b in self.impropers:
            # The 3rd atom is central
            if b.types[0] == n1 and b.types[1] == n2 and b.types[2] == n3 and b.types[3] == n4:
                ret.append(b)
            elif b.types[0] == n1 and b.types[3] == n2 and b.types[2] == n3 and b.types[1] == n4:
                ret.append(b)
            elif b.types[1] == n1 and b.types[0] == n2 and b.types[2] == n3 and b.types[3] == n4:
                ret.append(b)
            elif b.types[1] == n1 and b.types[3] == n2 and b.types[2] == n3 and b.types[0] == n4:
                ret.append(b)
            elif b.types[3] == n1 and b.types[0] == n2 and b.types[2] == n3 and b.types[1] == n4:
                ret.append(b)
            elif b.types[3] == n1 and b.types[1] == n2 and b.types[2] == n3 and b.types[0] == n4:
                ret.append(b)

        # Check if it is a duplicate that needs new parameters
        if len(ret) == 0:
            if ('x' in n1) or ('x' in n2) or ('x' in n3) or ('x' in n4):
                n1x = re.sub('x\d$', '', n1)
                n2x = re.sub('x\d$', '', n2)
                n3x = re.sub('x\d$', '', n3)
                n4x = re.sub('x\d$', '', n4)
                ret = self.improperParam(n1x, n2x, n3x, n4x)
                for param in ret:
                    param = deepcopy(param)
                    param.types[0] = n1
                    param.types[1] = n2
                    param.types[2] = n3
                    param.types[3] = n4
                    self.impropers.append(param)

        if len(ret) == 0:
            raise ValueError("Could not find improper parameters for %s-%s-%s-%s" % (n1, n2, n3, n4))

        return ret

class RTF:
    def __init__(self, filename):
        f = open(filename, "r")
        lines = f.readlines()
        f.close()

        self.types = []
        self.mass_by_type = dict()
        self.element_by_type = dict()
        self.type_by_name = dict()
        self.type_by_index = []
        self.index_by_name = dict()
        self.names = []
        self.charge_by_name = dict()
        self.bonds = []
        self.impropers = []
        self.typeindex_by_type = dict()
        self.netcharge = 0.

        aidx = 0

        for l in lines:
            if l.startswith("MASS "):
                k = l.split()
                at = k[2]
                self.mass_by_type[at] = float(k[3])
                self.element_by_type[at] = k[4]
                self.typeindex_by_type[at] = int(k[1])
                self.types.append(at)
            elif l.startswith("RESI "):
                k = l.split()
                self.netcharge = float(k[2])
            elif l.startswith("ATOM "):
                k = l.split()
                self.names.append(k[1])
                self.index_by_name[k[1]] = aidx
                self.type_by_index.append(k[2])
                self.type_by_name[k[1]] = k[2]
                self.charge_by_name[k[1]] = float(k[3])
                aidx += 1
            elif l.startswith("BOND "):
                k = l.split()
                self.bonds.append([self.index_by_name[k[1]], self.index_by_name[k[2]]])

            elif l.startswith("IMPR "):
                k = l.split()
                self.impropers.append([self.index_by_name[k[1]], self.index_by_name[k[2]], self.index_by_name[k[3]],
                                       self.index_by_name[k[4]]])

            self.natoms = aidx
        # if there weren't any "MASS" lines, we need to guess them
        typeindex=4000
        for idx in range(len(self.names)):
            atype = self.type_by_index[idx]
            name = self.names[idx]
            if atype not in self.element_by_type:
                self.element_by_type[atype] = self._guessElement(name)
                print("Guessing element %s for atom %s type %s" % (self.element_by_type[atype], name, atype))
            if atype not in self.mass_by_type:
                self.mass_by_type[atype] = self._guessMass(self.element_by_type[atype])

            if atype not in self.typeindex_by_type:
                self.typeindex_by_type[atype] = typeindex
                typeindex += 1
            if atype not in self.types:
                self.types.append(atype)

    def write(self, filename):
        f = open(filename, "w")
        print("* Charmm RTF built by HTMD parameterize version {}".format(htmdversion()), file=f)
        print("* ", file=f)
        print("  22     0", file=f)
        for a in self.types:
            print(
                "MASS %5d %s %8.5f %s" % (self.typeindex_by_type[a], a, self.mass_by_type[a], self.element_by_type[a]),
                file=f)
        print("\nAUTO ANGLES DIHE\n", file=f)
        print("RESI  MOL %8.5f" % self.netcharge, file=f)
        print("GROUP", file=f)
        for a in self.names:
            print("ATOM %4s %6s %8.6f" % (a, self.type_by_name[a], self.charge_by_name[a]), file=f)
        for a in self.bonds:
            print("BOND %4s %4s" % (self.names[a[0]], self.names[a[1]]), file=f)
        for a in self.impropers:
            print("IMPR %4s %4s %4s %4s" % (self.names[a[0]], self.names[a[1]], self.names[a[2]], self.names[a[3]]),
                  file=f)
        print("PATCH FIRST NONE LAST NONE", file=f)
        print("\nEND", file=f)
        f.close()

    def updateCharges(self, charges):
        if charges.shape[0] != self.natoms:
            raise ValueError("charge array not natoms in length")
        for i in range(self.natoms):
            name = self.names[i]
            self.charge_by_name[name] = charges[i]

    @staticmethod
    def _guessElement( name):
        import re
        name = re.sub('[0-9]*$', '', name)
        name = name.lower().capitalize()
        return name

    @staticmethod
    def _guessMass(element):
        from htmd.molecule import vdw
        return vdw.massByElement(element)


class AmberRTF(RTF):
    def __init__(self, mol, prepi, frcmod):
        f = open(prepi, "r")
        lines = f.readlines()
        f.close()
        f = lines

        # the prepi has the atoms re-ordered. Reorder the info based on the order in the mol

        self.names = []
        self.index_by_name = dict()
        idx = 0
        for i in mol.name:
            self.names.append(i)
            self.index_by_name[i] = idx
            idx += 1

        self.types = []
        self.mass_by_type = dict()
        self.element_by_type = dict()
        self.type_by_name = dict()
        self.type_by_index = dict()
        self.charge_by_name = dict()
        self.bonds = []
        self.impropers = []
        self.typeindex_by_type = dict()
        self.netcharge = 0.
        self.natoms = 0

        if f[4].split()[1] != "INT":
            raise ValueError("Invalid prepi format line 5")
        if f[5].strip() != "CORRECT     OMIT DU   BEG":
            raise ValueError("Invalid prepi format line 6")

        # Netcharge is 3rd term on 5th line
        cc = int(f[4].split()[2])

        ctr = 10
        while f[ctr].strip() != "":
            self.natoms += 1
            ff = f[ctr].split()
            ff[1] = ff[1].upper()
            idx = self.index_by_name[ff[1]]
            self.names.append(ff[1])
            self.type_by_name[ff[1]] = ff[2]
            self.type_by_index[idx] = ff[2]
            self.charge_by_name[ff[1]] = float(ff[10])
            if not (ff[2] in self.types):
                self.types.append(ff[2])
                self.typeindex_by_type[ff[2]] = 900 + len(
                    self.types) - 1  # add a big offset so it doesn't collide with real charm types
            self.element_by_type[ff[2]] = self._guessElement(ff[1])
            ctr += 1

        # Read improper section
        with open(prepi) as file:
            text = file.read()
        impropers = re.search('^IMPROPER\n(.+)\n\n', text, re.MULTILINE | re.DOTALL)  # extract improper section
        if impropers:
            impropers = impropers.group(1).split('\n')  # array of improper lines
            impropers = [improper.split() for improper in impropers]  # impropers by names
            for improper in impropers:
                improper_indices = [self.index_by_name[name.upper()] for name in improper]  # conv atom name to indices
                self.impropers.append(improper_indices)

        f = open(frcmod, "r")
        lines = f.readlines()
        f.close()
        section = None
        for ff in lines:
            ff = ff.strip()
            if ff == "MASS":
                section = "MASS"
            elif ff == "BOND":
                section = "BOND"
            elif ff == "ANGLE":
                section = "ANGLE"
            elif ff == "DIHE":
                section = "DIHE"
            elif ff == "IMPROPER":
                section = "IMPROPER"
            elif ff == "NONBON":
                section = "NONBON"
            else:
                if section == "MASS":
                    x = ff.split()
                    if len(x) >= 2:
                        self.mass_by_type[x[0]] = float(x[1])

# if __name__ == "__main__":
#
#  prm=PRM("ethanol.prm")
#  prm.write("ethanol_out.prm")
#  rtf=RTF("ethanol.rtf")
#  rtf.write("ethanol_out.rtf")
