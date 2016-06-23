import htmd.vmdviewer
import io
from htmd.vmdviewer import getCurrentViewer


class VMDGraphicObject(object):
    """A superclass from which VMD graphic objects (e.g. the convex hull) inherit."""
    counter = 1

    def __init__(self, data):
        """Generic creation method. Not useful for the user."""
        self.data = data
        self.n = self.counter
        self.script = io.StringIO()
        self.valid = True
        VMDGraphicObject.counter += 1

    def delete(self):
        """Undisplay and delete a graphic object."""
        if not self.valid:
            raise Exception("The object has been deleted already.")

        id = self.n
        vmd = htmd.vmdviewer.getCurrentViewer()
        cmd = "foreach i $htmd_graphics({:d}) {{ graphics top delete $i }}".format(id)
        vmd.send(cmd)
        cmd = "unset htmd_graphics({:d})".format(id)
        vmd.send(cmd)
        self.valid = False

    def _remember(self, s):
        n = self.n
        self.script.write("lappend htmd_graphics({:d}) [{:s}]\n".format(n, s))


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
            c1 = cc[v1, :]
            c1s = '{ ' + str(c1).strip('[]') + ' }'
            c2 = cc[v2, :]
            c2s = '{ ' + str(c2).strip('[]') + ' }'
            c3 = cc[v3, :]
            c3s = '{ ' + str(c3).strip('[]') + ' }'
            if solid:
                self._remember("draw triangle {:s} {:s} {:s}".format(c1s, c2s, c3s))
            else:
                self._remember("draw line {:s} {:s} {:s}".format(c1s, c2s, style))
                self._remember("draw line {:s} {:s} {:s}".format(c1s, c3s, style))
                self._remember("draw line {:s} {:s} {:s}".format(c2s, c3s, style))
                self.script.write("\n")

        cmd = self.script.getvalue()
        vmd = getCurrentViewer()
        vmd.send(cmd)


if __name__ == "__main__":
    from htmd import *
    from htmd import vmdgraphics

    """
    m=Molecule("3PTB")
    m.view()
    mf=m.copy()
    mf.filter("protein ")
    gh0=htmd.vmdgraphics.VMDConvexHull(mf)
    gh1=htmd.vmdgraphics.VMDConvexHull(mf,solid=True)

    """

    import doctest
    doctest.testmod()
