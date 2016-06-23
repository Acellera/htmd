import htmd.vmdviewer
import io
from htmd.vmdviewer import getCurrentViewer


class VMDGraphicObject(object):
    """A superclass from which VMD graphic objects (e.g. the convex hull) inherit."""
    counter = 1

    def __init__(self, data):
        """Generic creation method. Not useful for the user."""
        self.data = data
        self.id = self.counter
        self.valid = True
        self.counter += 1

    def delete(self):
        """Undisplay and delete a graphic object."""
        if not self.valid:
            raise Exception("The object has been deleted already.")

        id = self.id
        vmd = htmd.vmdviewer.getCurrentViewer()
        cmd = "foreach i $htmd_graphics({:d}) {{ graphics top delete $i }}".format(id)
        vmd.send(cmd)
        cmd = "unset htmd_graphics({:d})".format(id)
        vmd.send(cmd)
        self.valid = False

    def _remember(self, r, s):
        id = self.id
        r.write("lappend htmd_graphics({:d}) [{:s}]\n".format(id, s))


class VMDConvexHull(VMDGraphicObject):
    def __init__(self, m, style="", preamble="", solid=False):
        """Display the convex hull of the given molecule.

        For preamble and color, see http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.2/ug/node129.html

        The function returns an ID. To delete the objects, use the delete3DGraphics(<ID>) function.

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
        >>> htmd.vmdgraphics.VMDConvexHull(mf)
        """

        from scipy.spatial import ConvexHull

        if m.coords.shape[2] != 1:
            raise Exception("Only one frame is supported")

        cc = m.coords[:, :, 0]
        hull = ConvexHull(cc)

        super().__init__(hull)
        # VMDGraphicObject.__init__(self, hull)

        r = io.StringIO()
        r.write(preamble + "\n")

        for i in range(hull.nsimplex):
            v1, v2, v3 = hull.simplices[i, :]
            c1 = cc[v1, :]
            c1s = '{ ' + str(c1).strip('[]') + ' }'
            c2 = cc[v2, :]
            c2s = '{ ' + str(c2).strip('[]') + ' }'
            c3 = cc[v3, :]
            c3s = '{ ' + str(c3).strip('[]') + ' }'
            if solid:
                self._remember(r, "draw triangle {:s} {:s} {:s}".format(c1s, c2s, c3s))
            else:
                self._remember(r, "draw line {:s} {:s} {:s}".format(c1s, c2s, style))
                self._remember(r, "draw line {:s} {:s} {:s}".format(c1s, c3s, style))
                self._remember(r, "draw line {:s} {:s} {:s}".format(c2s, c3s, style))
                r.write("\n")

        cmd = r.getvalue()
        vmd = getCurrentViewer()
        vmd.send(cmd)


if __name__ == "__main__":
    from htmd import *
    from htmd import vmdgraphics

    m = Molecule("3PTB")
    m.view()
    mf = m.copy()
    mf.filter("protein ")
    v=htmd.vmdgraphics.VMDConvexHull(mf)

    import doctest
    doctest.testmod()
