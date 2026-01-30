"""Provide functions for adding grid lines (of constant longitude/latitude) to
a plot in the polar stereographic projection in cartesian (x,y) coordinates.
"""

import matplotlib as mpl
import numpy as np

from .. import config, transform
from ..utils import gen as gen_utils


def xy_longitude_gridlines(axs, lons=None, lat_lims=None, line_kw=None,
                           psp_xy_kw=None):
    """Add longitude grid lines (meridians) to one or more set of axes in (x,y)
    coordinates. In the (x,y) coordinate system, these are straight lines.


    Parameters
    ----------
    axs : matplotlib.axes.Axes instance or array of such
        The axes on which to draw grid lines.


    Optional parameters
    -------------------
    If any of the following are None, default values are taken from the
    configuration.

    lons : iterable of float, length M
        Longitude values at which to draw grid lines, in degrees_east (0-360).

    lat_lims : length-2 iterable or length-M list of length-2 lists or array (M, 2)
        Latitude limits (min and max) to draw meridians between. If a single
        set of limits is provided, this is applied to all meridians. Otherwise,
        individual limits (lat_lims[m]) can be provided for each meridian.

    line_kw : dictionary
        Appearance of grid lines, as keyword arguments passed to
        matplotlib.collections.LineCollection.

    psp_xy_kw : dictionary
        Properties of the polar stereographic projection in (x,y) coordinates
        as keyword arguments passed to function lonlat_to_xy_npsp().

    """

    if lons is None:
        lons = config.defaults["grid_longitudes"]

    if lat_lims is None:
        lat_lims = config.defaults["grid_longitudes_lat_lims"]

    line_kw = gen_utils.set_default_kw(line_kw, config.defaults["grid_line_kw"])

    if psp_xy_kw is None:
        psp_xy_kw = config.defaults["psp_xy_kw"]

    lat_lims = np.sort(lat_lims)

    if np.ndim(lat_lims) == 1:
        lat_lims = np.array([lat_lims for m in range(len(lons))])

    lines = []
    # Loop over requested meridians:
    for m in range(len(lons)):

        # Transform the pair of coordinates (lon_0, lat_0) and (lon_1, lat_1)
        # specifying the end points of meridian m to (x_0, y_0) and (x_1, y_1):
        xm, ym = transform.lonlat_to_xy_npsp(np.array([lons[m]]*2),
                                             lat_lims[m], **psp_xy_kw)

        # Re-shape into format required by matplotlib.collections.LineCollection:
        lines.append(np.vstack((xm, ym)).T)

    for ax in np.array(axs).flatten():
        ax.add_collection(mpl.collections.LineCollection(lines, **line_kw))


def xy_full_latitude_gridlines(axs, lats=None, line_kw=None, psp_xy_kw=None):
    """Add full (0-360 degrees_east) latitude grid lines to one or more sets of
    axes in (x,y) coordinates. In the (x,y) coordinate system, these are circles.


    Parameters
    ----------
    axs : matplotlib.axes.Axes instance or array of such
        The axes on which to draw grid lines.


    Optional parameters
    -------------------
    If any of the following are None, default values are taken from the
    configuration.

    lats : iterable of float, length L
        Latitude values at which to draw grid lines, in degrees_north.

    line_kw : dictionary
        Appearance of grid lines, as keyword arguments passed to
        matplotlib.collections.PatchCollection.

    psp_xy_kw : dictionary
        Properties of the polar stereographic projection in (x,y) coordinates
        as keyword arguments passed to function lonlat_to_xy_npsp().

    """

    if lats is None:
        lats = config.defaults["grid_latitudes"]

    line_kw = gen_utils.set_default_kw(line_kw, config.defaults["grid_line_kw"])

    if psp_xy_kw is None:
        psp_xy_kw = config.defaults["psp_xy_kw"]

    # Create copy of psp_xy_kw with xyp = (0,0), because this is needed to
    # properly calculated circle radii, and use the provided or default value
    # of xyp to determine the center:
    if "xy_offset" in psp_xy_kw.keys():
        circle_center = psp_xy_kw["xy_offset"]
    else:
        circle_center = config.defaults["xy_offset"]

    psp_xy_kw_copy = psp_xy_kw.copy()
    psp_xy_kw_copy["xy_offset"] = (0., 0.)

    # Unless specified these should not be filled:
    if "facecolor" not in line_kw.keys():
        line_kw["facecolor"] = "none"
    if "fill" in line_kw.keys():
        del line_kw["fill"]

    # Loop over each latitude:
    circles = []
    for j in range(len(lats)):

        # Only the circle radius in xy-space is required: transform one set of
        # coordinates with arbitrary longitude and specified latitude:
        xj, yj = transform.lonlat_to_xy_npsp(np.array([90]),
                                             np.array([lats[j]]),
                                             **psp_xy_kw_copy)

        circles.append(mpl.patches.Circle(circle_center,
                                          radius=np.sqrt(xj[0]**2 + yj[0]**2)))

    for ax in np.array(axs).flatten():
        ax.add_collection(mpl.collections.PatchCollection(circles, **line_kw))


def xy_latitude_gridlines(axs, lon_lims, lats=None, line_kw=None,
                          psp_xy_kw=None):
    """Add latitude grid lines between specified longitude limits to one or
    more sets of axes in (x,y) coordinates. In the (x,y) coordinate system,
    these are, in general, arcs.

    Note that this adds matplotlib.patches.Arc instances to each axes. If full
    latitude circles (0 - 360 degrees_east) are required, use the function
    xy_full_latitude_gridlines() which adds matplotlib.patches.Circle patches
    which are more efficiently drawn and are precisely circular.


    Parameters
    ----------
    axs : matplotlib.axes.Axes instance or array of such
        The axes on which to draw grid lines.

    lon_lims : length-2 iterable or length-L list of length-2 lists or array (L, 2)
        Longitude limits (min and max) to draw latitude arcs between. If a
        single set of limits is provided, this is applied to all latitude circles.
        Otherwise, individual limits (lon_lims[l]) can be provided for each
        latitude.


    Optional parameters
    -------------------
    If any of the following are None, default values are taken from the
    configuration.

    lats : iterable of float, length L
        Latitude values at which to draw grid lines, in degrees_north.

    line_kw : dictionary
        Appearance of grid lines, as keyword arguments passed to
        matplotlib.collections.LineCollection.

    psp_xy_kw : dictionary
        Properties of the polar stereographic projection in (x,y) coordinates
        as keyword arguments passed to function lonlat_to_xy_npsp().

    """

    if lats is None:
        lats = config.defaults["grid_latitudes"]

    line_kw = gen_utils.set_default_kw(line_kw, config.defaults["grid_line_kw"])

    if psp_xy_kw is None:
        psp_xy_kw = config.defaults["psp_xy_kw"]

    # matplotlib.patches.Arc angle parameter is measured with 0 corresponding
    # to positive x-axis; account for longitude offset and central longitude
    # (rotation) of plot:
    arc_angle = -90.
    if "central_longitude" in psp_xy_kw.keys():
        arc_angle += psp_xy_kw["central_longitude"]
    else:
        arc_angle += config.defaults["central_longitude"]

    # Create copy of psp_xy_kw with xy_offset = (0, 0), because this is needed
    # to properly calculated arc radii, and use the provided or default value of
    # xy_offset to determine the center:
    if "xy_offset" in psp_xy_kw.keys():
        arc_center = psp_xy_kw["xy_offset"]
    else:
        arc_center = config.defaults["xy_offset"]

    psp_xy_kw_copy = psp_xy_kw.copy()
    psp_xy_kw_copy["xy_offset"] = (0., 0.)

    # Arcs cannot be filled and will raise an exceptionif 'fill' is passed:
    if "fill" in line_kw.keys():
        del line_kw["fill"]
    if not "facecolor" in line_kw.keys():
        line_kw["facecolor"] = "none"

    if "color" in line_kw.keys():
        if not "edgecolor" in line_kw.keys():
            line_kw["edgecolor"] = line_kw["color"]
        del line_kw["color"]

    if np.ndim(lon_lims) == 1:
        lon_lims = np.array([lon_lims for j in range(len(lats))])

    # Loop over latitudes:
    for j in range(len(lats)):
        
        # Only the arc radius in xy-space is required: transform one set of
        # coordinates with arbitrary longitude and specified latitude:
        xj, yj = transform.lonlat_to_xy_npsp(np.array([90]),
                                             np.array([lats[j]]),
                                             **psp_xy_kw_copy)

        radius = np.sqrt(xj[0]**2 + yj[0]**2)

        # matplotlib Arcs cannot be added as a patch collection; see:
        #     https://github.com/matplotlib/matplotlib/issues/11266
        #     [Accessed: 10 Aug 2022]
        #
        # **UPDATE 18 Sep 2025**
        # Potentially fixed in matplotlib (TODO: use patch collection here):
        #     https://github.com/matplotlib/matplotlib/pull/23340
        #     [Accessed 18 Sep 2025]
        #
        for ax in np.array(axs).flatten():
            ax.add_patch(mpl.patches.Arc(arc_center, 2.*radius, 2.*radius,
                                         angle=arc_angle, theta1=lon_lims[j][0],
                                         theta2=lon_lims[j][1], **line_kw))


def xy_gridlines(axs, lons=None, lat_lims=None, lats=None, lon_lims=None,
                 line_kw=None, psp_xy_kw=None):
    """Add longitude and latitude grid lines to one or more plot axes in polar
    stereographic projection (x,y) coordinates.


    Parameters
    ----------
    axs : matplotlib.axes.Axes instance or array of such
        The axes on which to draw grid lines.


    Optional parameters
    -------------------
    If any of the following are None, default values are taken from the
    configuration.

    lons : iterable of float, length M
        Longitude values at which to draw grid lines, in degrees_east (0-360).

    lat_lims : length-2 iterable or length-M list of length-2 lists or array (M, 2)
        Latitude limits (min and max) to draw meridians between. If a single
        set of limits is provided, this is applied to all meridians. Otherwise,
        individual limits (lat_lims[m]) can be provided for each meridian.

    lats : iterable of float, length L
        Latitude values at which to draw grid lines, in degrees_north.

    lon_lims : length-2 iterable or length-L list of length-2 lists or array (L, 2)
        Longitude limits (min and max) to draw latitude arcs between. If a
        single set of limits is provided, this is applied to all latitude circles.
        Otherwise, individual limits (lon_lims[l]) can be provided for each
        latitude.

    line_kw : dictionary
        Appearance of grid lines, as keyword arguments passed to
        matplotlib.collections.LineCollection.

    psp_xy_kw : dictionary
        Properties of the polar stereographic projection in (x,y) coordinates
        as keyword arguments passed to function lonlat_to_xy_npsp().

    """

    if lons is None:
        lons = config.defaults["grid_longitudes"]

    if lat_lims is None:
        lat_lims = config.defaults["grid_longitudes_lat_lims"]

    if lats is None:
        lats = config.defaults["grid_latitudes"]

    if lon_lims is None:
        lon_lims = config.defaults["grid_latitudes_lon_lims"]

    line_kw = gen_utils.set_default_kw(line_kw, config.defaults["grid_line_kw"])

    if psp_xy_kw is None:
        psp_xy_kw = config.defaults["psp_xy_kw"]

    # Gridline style and projection specification are common to all gridlines
    # and use the same keyword (kw) arguments:
    common_kw = {"line_kw": line_kw, "psp_xy_kw": psp_xy_kw}

    if int(abs(lon_lims[1] - lon_lims[0])) == 360:

        # Plot latitude lines as Circle patches:
        xy_full_latitude_gridlines(axs, lats=lats, **common_kw)

        # Create longitude gridlines (meridians)
        xy_longitude_gridlines(axs, lons=lons, lat_lims=lat_lims, **common_kw)

    else:

        # Plot latitude lines as Arc patches between longitude limits:
        xy_latitude_gridlines(axs, lon_lims, lats=lats, **common_kw)

        # Create longitude gridlines (meridians) but only where latitude lines
        # are drawn (i.e., limit longitude range):
        lons = np.array(lons)[np.where(  (np.array(lons) >= min(lon_lims))
                                       & (np.array(lons) <= max(lon_lims)))]

        xy_longitude_gridlines(axs, lons=lons, lat_lims=lat_lims, **common_kw)

