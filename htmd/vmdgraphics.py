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
        <htmd.vmdgraphics.VMDConvexHull object at ...
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
