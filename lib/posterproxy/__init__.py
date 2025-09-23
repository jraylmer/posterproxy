"""Implementation of the POlar STEReographic map PROjection in cartesian XY
coordinates for use with matplotlib plotting.

Points P on the surface of Earth in the northern hemisphere are projected
P(lon,lat) -> P'(x,y) onto the equatorial plane Eq, where P' is the
intersection of Eq with the line connecting the opposite pole (S) and P:

                                      N
                                  , - ~ - ,     P(lon,lat)
                              , '     |     ' ,*
                            ,         |     .' ^,
                           ,          |   .'  /  ,
                          ,           | .'   /    ,
                       <--------------+-----*--------> Eq
                          ,           |    / P'(x,y)
                           ,          |   /      ,
                            ,         |  /      ,
                              ,       | /     ,
                                 ,    |/  , '
                                    ' ~ '
                                      S

Currently, only the north polar stereographic projection (NPSP) is explicitly
implemented (though the south PSP is possible by manually adapting the inputs).

The package also includes functions for adding grid lines and land/coast features
using, for the latter, data from Natural Earth.


Basic usage
-----------
Given a set of longitude and latitude coordinates in 2D arrays called
lon and lat, with an array data located at these coodinates:

    >>> import posterproxy as psp
    >>> x, y = psp.lonlat_to_xy_npsp(lon, lat)

This generates cartesian coordinates (x,y) corresponding to the NPSP, that can
be used in regular plotting functions with matplotlib, such as:

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> psp.prepare_axes(ax)
    >>> ax.contourf(x, y, data)
    >>> psp.xy_land_overlay(ax)
    >>> psp.xy_gridlines(ax)

where the intermediate step psp.prepare_axes(ax) is not strictly necessary but
sets the appropriate axes display (equal aspect, no tick labels) and limits.

"""

__author__  = "Jake Aylmer"
__version__ = "0.1.0"

import numpy as np

from . import config

config._parse_config()

from .transform import lonlat_to_xy_npsp
from .transform import rotate_vectors_npsp
from .transform import spsp_from_global_data

from .decor.naturalearth import xy_land_overlay

from .decor.gridlines import xy_longitude_gridlines, xy_full_latitude_gridlines
from .decor.gridlines import xy_latitude_gridlines , xy_gridlines


def prepare_axes(axs, xy_extent_0=config.defaults["xy_extent_0"],
                      xy_extent_1=config.defaults["xy_extent_1"]):
    """Prepare one or more sets of matplotlib Axes instances for plotting with
    polar stereographic projection (x,y) coordinates.

    Specifically, this fixes the axes limits, sets the axes aspect ratio to be
    'equal', removes all axes ticks and tick labels, sets all axes spines to be
    visible and on top of plot elements, and deactivates the (cartesian) grid.


    Parameters
    ----------
    axs : matplotlib.axes.Axes instance, or array of such
        The axes(s) to set up.


    Optional parameters
    -------------------
    xy_extent_0, xy_extent_1 : length-2 iterables of float
        The xy coordinates at the lower-left and upper-right corners,
        respectively. Defaults are taken from the configuration. Note that
        xy_extent_0 = -xy_extent_1 is required. There should generally
        never be any need to change the default values (the actual x,y
        coordinates are arbitrary).

    """

    if not (    xy_extent_0[0] == -xy_extent_1[0]
            and xy_extent_0[1] == -xy_extent_1[1]):
        raise ValueError("xy_extent_0 must equal minus xy_extent_1")

    for ax in np.array(axs).flatten():

        # Remove axes ticks and tick labels:
        ax.tick_params(which="both", axis="both", left=False, right=False,
                       top=False, bottom=False, labelleft=False,
                       labelright=False, labeltop=False, labelbottom=False)

        # Cartesian grid lines off and 'equal' aspect ratio by default:
        ax.grid(False)
        ax.set_aspect("equal")

        # Draw axes spines (outline/border) above all plot elements
        # (set the z order to very high value):
        for spine in ax.spines:
            ax.spines[spine].set_visible(True)
            ax.spines[spine].set_zorder(1E10)

        ax.set_xlim(xy_extent_0[0], xy_extent_1[0])
        ax.set_ylim(xy_extent_0[1], xy_extent_1[1])

