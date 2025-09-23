"""Functions for dealing with matplotlib objects."""

import copy

import matplotlib as mpl
from matplotlib.path import Path as mpl_Path
from matplotlib.patches import PathPatch
import numpy as np


def assign_collection(axs, collection):
    """Add an artist collection to a (set of) matplotlib axes.

    Parameters
    ----------
    axs : matplotlib.axes.Axes instance, or iterable of such
        The axes(s) to add collection to.

    collection : matplotlib.collections.Collection, or instance
                 of a subclass of such (e.g., PatchCollection).
        The collection to be added to (each axes in) axs.

    """

    # matplotlib does not support adding the same artists
    # to multiple containers, so they must be (deep) copied
    # before assigning to each axes instance:
    for axes in np.array(axs).flatten():
        axes.add_collection(copy.deepcopy(collection))


def polygon_with_interiors(polys):
    """Construct a matplotlib patch representing a polygon with interior rings
    that will be drawn as holes in the exterior.

    This function has been adapted from GitHub user 'yohai's code:

    https://gist.github.com/yohai/81c5854eaa4f8eb5ad2256acd17433c8
    [Last accessed: 20 September 2025]

    which in turn was inspired by:

    https://sgillies.net/2010/04/06/painting-punctured-polygons-with-matplotlib.html
    [Last accessed: 20 September 2025]


    Parameters
    ----------
    polys: iterable of arrays, each of shape (*, 2)
        The coordinates of the vertices of the rings defining the polygon. The
        first array, polys[0], is taken as the exterior, with the rest being
        interior rings. Ring p has Np coordinates, and the horizontal/vertical
        coordinates are taken as polys[p][n] = [x, y] for n = 0..(Np-1).


    Returns
    -------
    matplotlib PathPatch instance
        This may then be added to a plot; e.g., using ax.add_patch().

    """

    def _reorder(poly, cw=True):
        """Reorders the polygon to run clockwise (cw) or counter-clockwise
        (ccw). It calculates whether a polygon is cw or ccw by summing
        (x2-x1)*(y2+y1) for all edges of the polygon [1].

        [1] https://stackoverflow.com/a/1165943/898213
            [Last accessed: 25 September 2025]
        """

        # Close polygon if not closed:
        if not np.allclose(poly[0,:], poly[-1,:]):
            poly = np.concatenate((poly, poly[[0],:]), axis=0)

        direction = (  (poly[:,0] - np.roll(poly[:,0], 1))
                     * (poly[:,1] + np.roll(poly[:,1], 1))).sum() < 0

        if direction == cw:
            return poly
        else:
            return poly[::-1,:]

    def _ring_coding(n):
        """Returns a length-n list of this format:
        [MOVETO, LINETO, LINETO, ..., LINETO, LINETO, CLOSEPOLY]
        """
        codes = [mpl_Path.LINETO] * n
        codes[0] = mpl_Path.MOVETO

        if n > 1:
            codes[-1] = mpl_Path.CLOSEPOLY

        return codes

    ccw = [True] + ([False] * (len(polys) - 1))
    polys = [_reorder(poly, c) for poly, c in zip(polys, ccw)]
    codes = np.concatenate([_ring_coding(p.shape[0]) for p in polys])
    vertices = np.concatenate(polys, axis=0)

    return PathPatch(mpl_Path(vertices, codes))

