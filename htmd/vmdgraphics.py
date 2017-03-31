# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import htmd.vmdviewer
import io
from htmd.vmdviewer import getCurrentViewer


class VMDGraphicObject(object):
    """A superclass from which VMD graphic objects (e.g. the convex hull) inherit."""
    counter = 1

    def __init__(self, data):
        """Generic creation method. Not useful for the user."""
        self.n = self.counter
        self.data = data
        self.script = io.StringIO()
        self.valid = True
        VMDGraphicObject.counter += 1

    def delete(self):
        """Undisplay and delete a graphic object."""
        if not self.valid:
            raise Exception("The object has been deleted already.")

        n = self.n
        vmd = htmd.vmdviewer.getCurrentViewer()
        cmd = """
            set htmd_tmp $htmd_graphics_mol({0:d})
            foreach i $htmd_graphics({0:d}) {{ graphics $htmd_tmp delete $i }}
            unset htmd_graphics({0:d})
            unset htmd_graphics_mol({0:d})
        """.format(n)
        vmd.send(cmd)
        self.valid = False

    def _remember(self, s):
        # We can't get data back from VMD, so let it remember what's to be deleted
        n = self.n
        self.script.write("lappend htmd_graphics({:d}) [{:s}]\n".format(n, s))

    @staticmethod
    def tq(v):
        """Quote a numpy 3-vector to a TCL list."""
        return '{ ' + str(v).strip('[]') + ' }'


class VMDConvexHull(VMDGraphicObject):
    def __init__(self, m, style="", preamble="", solid=False):
        """Display the convex hull of the given molecule.

        For preamble and color, see http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node129.html

        The function returns an instance of VMDGraphicsObject. To delete it, use the delete() method.

        Parameters
        ----------
        m: Molecule
            The object of which to show the hull (only 1 frame)
        style: str
            Style for wireframe lines
        preamble: str
            Commands (material, color) to be prefixed to the output.
            E.g.: "draw color red; graphics top materials on; graphics top material Transparent".
            Note that this affects later preamble-less commands.
        solid: bool
            Solid or wireframe

        Examples
        --------
        >>> m=Molecule("3PTB")
        >>> m.view()
        >>> mf=m.copy()
        >>> mf.filter("protein ")
        >>> gh=htmd.vmdgraphics.VMDConvexHull(mf)  # doctest: +ELLIPSIS
        >>> gh.delete()
        """

        from scipy.spatial import ConvexHull

        if m.coords.shape[2] != 1:
            raise Exception("Only one frame is supported")

        cc = m.coords[:, :, 0]
        hull = ConvexHull(cc)

        super().__init__(hull)
        self.script.write(preamble + "\n")

        for i in range(hull.nsimplex):
            v1, v2, v3 = hull.simplices[i, :]
            c1s = VMDGraphicObject.tq(cc[v1, :])
            c2s = VMDGraphicObject.tq(cc[v2, :])
            c3s = VMDGraphicObject.tq(cc[v3, :])
            if solid:
                self._remember("draw triangle {:s} {:s} {:s}".format(c1s, c2s, c3s))
            else:
                self._remember("draw line {:s} {:s} {:s}".format(c1s, c2s, style))
                self._remember("draw line {:s} {:s} {:s}".format(c1s, c3s, style))
                self._remember("draw line {:s} {:s} {:s}".format(c2s, c3s, style))
                self.script.write("\n")

        self.script.write("set htmd_graphics_mol({:d}) [molinfo top]".format(self.n))
        cmd = self.script.getvalue()
        vmd = getCurrentViewer()
        vmd.send(cmd)


class VMDBox(VMDGraphicObject):
    def __init__(self, box, color='red'):
        """ Displays a box in VMD as lines of its edges.

        The function returns an instance of VMDGraphicsObject. To delete it, use the delete() method.

        Parameters
        ----------
        box : list
            The min and max positions of the edges. Given as [xmin, xmax, ymin, ymax, zmin, zmax]
        color : str
            Color of the lines
        """
        super().__init__(box)
        xmin, xmax, ymin, ymax, zmin, zmax = box
        mi = [xmin, ymin, zmin]
        ma = [xmax, ymax, zmax]

        self._remember('draw materials off')
        self._remember('draw color {}'.format(color))

        self._remember('draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], mi[1], mi[2], ma[0], mi[1], mi[2]))
        self._remember('draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], mi[1], mi[2], mi[0], ma[1], mi[2]))
        self._remember('draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], mi[1], mi[2], mi[0], mi[1], ma[2]))

        self._remember('draw line "{} {} {}" "{} {} {}"\n'.format(ma[0], mi[1], mi[2], ma[0], ma[1], mi[2]))
        self._remember('draw line "{} {} {}" "{} {} {}"\n'.format(ma[0], mi[1], mi[2], ma[0], mi[1], ma[2]))

        self._remember('draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], ma[1], mi[2], ma[0], ma[1], mi[2]))
        self._remember('draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], ma[1], mi[2], mi[0], ma[1], ma[2]))

        self._remember('draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], mi[1], ma[2], ma[0], mi[1], ma[2]))
        self._remember('draw line "{} {} {}" "{} {} {}"\n'.format(mi[0], mi[1], ma[2], mi[0], ma[1], ma[2]))

        self._remember('draw line "{} {} {}" "{} {} {}"\n'.format(ma[0], ma[1], ma[2], ma[0], ma[1], mi[2]))
        self._remember('draw line "{} {} {}" "{} {} {}"\n'.format(ma[0], ma[1], ma[2], mi[0], ma[1], ma[2]))
        self._remember('draw line "{} {} {}" "{} {} {}"\n'.format(ma[0], ma[1], ma[2], ma[0], mi[1], ma[2]))
        self.script.write("\n")
        self.script.write("set htmd_graphics_mol({:d}) [molinfo top]".format(self.n))
        cmd = self.script.getvalue()
        vmd = getCurrentViewer()
        vmd.send(cmd)


class VMDSphere(VMDGraphicObject):
    def __init__(self, xyz, color='red', radius=1):
        """ Displays a sphere in VMD.

        The function returns an instance of VMDGraphicsObject. To delete it, use the delete() method.

        Parameters
        ----------
        xyz : list
            The center of the sphere
        color : str
            Color of the sphere
        radius : float
            The radius of the sphere
        """
        super().__init__(xyz)
        #self._remember('draw materials off')
        self._remember('draw color {}'.format(color))

        self._remember('draw sphere "{} {} {}" radius {}\n'.format(xyz[0], xyz[1], xyz[2], radius))
        self.script.write("\n")
        self.script.write("set htmd_graphics_mol({:d}) [molinfo top]".format(self.n))
        cmd = self.script.getvalue()
        vmd = getCurrentViewer()
        vmd.send(cmd)


class VMDText(VMDGraphicObject):
    def __init__(self, text, xyz, color='red'):
        """ Displays a text in VMD.

        The function returns an instance of VMDGraphicsObject. To delete it, use the delete() method.

        Parameters
        ----------
        text : str
            The text
        xyz : list
            The position of the text
        color : str
            Color of the text
        """
        super().__init__(xyz)
        self._remember('draw materials off')
        self._remember('draw color {}'.format(color))

        self._remember('draw text "{} {} {}" "{}"\n'.format(xyz[0], xyz[1], xyz[2], text))
        self.script.write("\n")
        self.script.write("set htmd_graphics_mol({:d}) [molinfo top]".format(self.n))
        cmd = self.script.getvalue()
        vmd = getCurrentViewer()
        vmd.send(cmd)


class VMDIsosurface(VMDGraphicObject):
    def __init__(self, arr, vecMin, vecRes, color='red'):
        """ Displays an isosurface in VMD

        The function returns an instance of VMDGraphicsObject. To delete it, use the delete() method.

        Parameters
        ----------
        arr: np.ndarray
                array with volumetric data
        filename: str
                string with the name of the cubefile
        vecMin: np.ndarray
                3D vector denoting the minimal corner of the grid
        vecRes: np.ndarray
                3D vector denoting the resolution of the grid in each dimension
        """
        super().__init__(arr)
        from htmd.util import tempname
        import numpy as np
        import os
        filename = tempname(suffix='.cube')
        outFile = open(filename, 'w')

        # conversion to gaussian units
        L = 1 / 0.52917725
        gauss_bin = vecRes * L
        # minCorner = 0.5*L*(vecMin - vecMax + vecRes)
        minCorner = L * (vecMin + 0.5 * vecRes)

        ngrid = arr.shape

        # write header
        outFile.write("CUBE FILE\n")
        outFile.write("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n")
        outFile.write("%5d %12.6f %12.6f %12.6f\n" % (1, minCorner[0], minCorner[1], minCorner[2]))
        outFile.write("%5d %12.6f %12.6f %12.6f\n" % (ngrid[0], gauss_bin[0], 0, 0))
        outFile.write("%5d %12.6f %12.6f %12.6f\n" % (ngrid[1], 0, gauss_bin[1], 0))
        outFile.write("%5d %12.6f %12.6f %12.6f\n" % (ngrid[2], 0, 0, gauss_bin[2]))
        outFile.write("%5d %12.6f %12.6f %12.6f %12.6f\n" % (1, 0, minCorner[0], minCorner[1], minCorner[2]))

        # main loop
        cont = 0

        for i in range(ngrid[0]):
            for j in range(ngrid[1]):
                for k in range(ngrid[2]):
                    outFile.write("%13.5g" % arr[i][j][k])
                    if np.mod(cont, 6) == 5:
                        outFile.write("\n")
                    cont += 1

        outFile.close()
        vmd = getCurrentViewer()

        vmd.send('mol new {} type cube first 0 last -1 step 1 waitfor 1 volsets {{0 }}'.format(filename))
        vmd.send('mol modstyle 0 top Isosurface -0.4 0 2 0 1 1')
        vmd.send('set htmd_graphics_mol({:d}) [molinfo top]'.format(self.n))

        os.unlink(filename)

    def delete(self):
        vmd = getCurrentViewer()
        vmd.send('mol delete htmd_graphics_mol({:d})'.format(self.n))


if __name__ == "__main__":
    from htmd import *
    from htmd import vmdgraphics

    """
    from htmd import *
    import htmd.vmdgraphics

    m=Molecule("3PTB")
    m.view()
    mf=m.copy()
    mf.filter("protein")
    mgh=htmd.vmdgraphics.VMDConvexHull(mf)

    n=Molecule("1JNO")
    n.view()
    nf=n.copy()
    nf.filter("protein")
    ngh=htmd.vmdgraphics.VMDConvexHull(nf,solid=True)

    """

    import doctest
    doctest.testmod()
