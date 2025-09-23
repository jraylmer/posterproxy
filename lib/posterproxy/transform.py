"""Implementation of the polar stereographic projection."""

import numpy as np
from . import config


def _stereographic_projection_polar_coordinates(lon, lat, central_longitude=0.):
    """Transforms the coordinates of a given point P (or many points at once) on
    the surface of Earth onto a plane perpendicular to Earth's axis of rotation
    and passing through the equator (i.e., the polar stereographic projection).

    Note: this transformation works in units where the Earth radius is 1 and
    assumes the input coordinates are in degrees.

    The transformed coordinates are given in the polar coordinate system
    (r, theta) on the aforementioned plane, in which r is the distance from the
    origin (O; the pole) to the projected point, P', and theta is the angle
    between the axis oriented eastwards and the vector OP'.


    Parameters
    ----------
    lon, lat : float or array
        The longitude(s) and latitude(s) respectively of the point(s) P to be
        transformed, in degrees_east and degrees_north, respectively.


    Optional parameters
    -------------------
    central_longitude : float, default = 0.
        This parameter is analagous to that of the same name in
        cartopy.crs.NorthPolarStereo(), and defines the orientation of the data
        in P' space.

        By default, coordinates with longitude 0. have theta = 270. (i.e., on a
        map would appear at the bottom centre). By increasing this parameter to,
        e.g., 180., a rotation transformation is applied to the (r, theta)
        coordinates such that points with longitude 0. would have theta = 90.
        (i.e., on a map would now appear at the top center).


    Returns
    -------
    r, theta : float or array
        The transformed polar coordinates as described above.


    Notes
    -----
    The transformation is not valid for the opposite pole so a division-by-zero
    warning will be raised for lat = -90 (for north polar projection;
    effectively,such points are projected onto the plane at infinite radius).
    This may cause plotting artefacts (the NPSP should not be used for plotting
    data at latitudes less than around 50 degrees from the pole, anyway, but note
    that global datasets should thus be sub-setted to the plot region before
    plotting them).

    """

    r     = np.cos(lat * np.pi/180.) / (1 + np.sin(lat * np.pi/180.))
    theta = ((lon % 360.) - 90. - central_longitude) % 360.

    return r, theta


def _stereographic_projection_cartesian_coordinates(lon, lat, central_longitude=0.):
    """Transform the coordinates of a given point P (or many points at once) on
    the surface of Earth onto a plane perpendicular to Earth's axis of rotation
    and passing through the equator (i.e., the polar stereographic projection).

    The transformed coordinates are given in a cartesian coordinate system
    (x, y) on the aforementioned plane, in which x (y) is the horizontal
    (vertical) distance from the origin (O; the pole) to the projected point P'.

    Note that the actual projection transforms P to the plane in polar
    coordinates before transforming to cartesian coordinates. Units are used
    such that the Earth radius is 1, and the input coordinates are assumed to
    be in degrees.

    Use the wrapper function lonlat_to_xy_npsp() to transform longitudes and
    latitudes to NPSP with specified scale factors.


    Parameters
    ----------
    lon, lat : float or array
        The longitude(s) and latitude(s) respectively of the point(s) P to be
        transformed, in degrees_east and degrees_north, respectively.


    Optional parameters
    ------------------
    central_longitude : float, default = 0.0
        This parameter is analgous to that of the same name in
        cartopy.crs.NorthPolarStereo(), and defines the orientation of the data
        in P' space.

        By default, coordinates with longitude 0. have theta = 270. (i.e., on a
        map would appear at the bottom centre). By increasing this parameter to,
        e.g., 180., a rotation transformation is applied to the (r, theta)
        coordinates such that points with longitude 0. would have theta = 90.
        (i.e., on a map would now appear at the top center).


    Returns
    -------
    x, y : float or array
        The transformed polar coordinates as described above.

    """

    # Carry out the projection --> result is in polar coordinates:
    r, theta = _stereographic_projection_polar_coordinates(lon, lat,
        central_longitude=central_longitude)

    # Transform to (x,y) system and return:
    return r * np.cos(theta * np.pi/180.), r*np.sin(theta * np.pi/180.)


def lonlat_to_xy_npsp(lon, lat,
                      extent_latitude   = config.defaults["extent_latitude"],
                      central_longitude = config.defaults["central_longitude"],
                      xy_extent_0       = config.defaults["xy_extent_0"],
                      xy_extent_1       = config.defaults["xy_extent_1"],
                      xy_offset         = config.defaults["xy_offset"]):
    """Transform (lon, lat) coordinates into a north-polar stereographic
    projection (NPSP) defined in cartesian coordinates (x,y).


    Parameters
    ----------
    lon, lat : float or arrays of float
        Longitude and latitude coordinates in degrees_east and degrees_north,
        respectively. If arrays, they must be the same shape.


    Optional parameters
    -------------------
    Default values for the following parameters are set in the configuration:

    extent_latitude : float
        The latitude in degrees_north corresponding to the outer limit of the
        plot boundary. This sets the scale of the (x,y) coordinate system.

    central_longitude : float, default = 0.
        Sets the orientation of the plot, analogous to
        cartopy.crs.NorthPolarStereo(central_longitude=0.).

    The following parameters determine the scale and offset of the (x,y)
    coordinate system. The xy_extent_* parameters should rarely need to be
    changed since the scales of x and y are arbitrary, but for plotting it is
    essential that the coordinate axes match these. Also it is currently
    required that xy_extent_0 = -xy_extent_1 [so that the origin is (0,0)].
    A ValueError is raised if this is not the case.

    xy_extent_0 : 2-tuple of float
        The (x,y) coordinates of the lower-left corner of the plot. The
        standard value is (-1, -1).

    xy_extent_1 : 2-tuple of float
        The (x,y) coordinates of the upper-right corner of the plot. The
        standard value is (1, 1).

    xy_offset : 2-tuple of float
        The (x,y) coordinate shift by which to offset the overall plot area.
        Default is (0, 0), i.e., no offset.

        This parameter may be used to offset the centre of the plot, but note
        that this has no effect on the coordinate projection: it simply
        translates coordinates so that pole appears at coordinates xy_offset
        (an analogy would be making a standard stereographic projection with the
        pole in the centre of the plot, then using the matplotlib interactive
        toolbar to manually click and drag the data around).


    Returns
    -------
    x, y : float or arrays of the same shape as lon and lat
        The coordinate(s) in the transformed [and translated if if applicable]
        coordinate system.

    """

    if not (    xy_extent_0[0] == -xy_extent_1[0]
            and xy_extent_0[1] == -xy_extent_1[1]):
        raise ValueError("xy_extent_0 must equal -xy_extent_1")

    lon_360 = lon.copy() % 360.

    # Project (lon, lat) --> (r, theta) (2D polar-coordinate system):
    r, theta = _stereographic_projection_polar_coordinates(lon_360, lat,
        central_longitude=central_longitude)

    # Scale data such that the outer data limit matches the specified
    # extent_latitude (i.e., stretch transformation):
    r_lim, _ = _stereographic_projection_polar_coordinates(lon_360,
        extent_latitude, central_longitude=central_longitude)

    r /= r_lim

    # Transform to (x,y) coordinates (pole in the center):
    #
    # Note: this is why xy_extent_0 must equal -xy_extent_1: not clear how to
    # ----- handle cases where this is not so. Maybe a TODO item but not hugely
    #       important since (x,y) are arbitrary (any use cases?)
    #
    x = abs(xy_extent_1[0]) * r * np.cos( theta * np.pi/180. )
    y = abs(xy_extent_1[1]) * r * np.sin( theta * np.pi/180. )

    # Account for offset [simply translate coordinates in (x,y) space]:
    x += xy_offset[0] - xy_extent_0[0] - .5 * (xy_extent_1[0] - xy_extent_0[0])
    y += xy_offset[1] - xy_extent_0[1] - .5 * (xy_extent_1[1] - xy_extent_0[1])

    return x, y


def rotate_vectors_npsp(lon, u_lonlat, v_lonlat,
                        central_longitude=config.defaults["central_longitude"]):
    """Transform 2D vector components (u,v) defined with respect to a regular
    longitude/latitude grid basis onto a polar stereographic grid (u',p').


    Parameters
    ----------
    lon : float or array of float
        Longitude coordinates in degrees_east.

    u_lonlat, v_lonlat : flot or array of float, same shape as lon
        Vector components defined on the longitude/latitude grid.


    Optional parameters
    -------------------
    central_longitude : float
        Analogous to the parameter of the same name in
        cartopy.crs.NorthPolarStereo(). Default set in configuration.


    Returns
    -------
    u_xy, v_xy : arrays of float, same shape as lon
        Vector components on the polar stereographic grid.

    """

    rot_angle = (lon % 360.) - central_longitude

    cos_rot_angle = np.cos( rot_angle * np.pi/180. )
    sin_rot_angle = np.sin( rot_angle * np.pi/180. )

    u_xy = u_lonlat * cos_rot_angle - v_lonlat * sin_rot_angle
    v_xy = u_lonlat * sin_rot_angle + v_lonlat * cos_rot_angle

    return u_xy, v_xy


def spsp_from_global_data(lon, lat, data,
                          extent_latitude   = 65.,
                          central_longitude = 180.0,
                          xy_extent_0       = config.defaults["xy_extent_0"],
                          xy_extent_1       = config.defaults["xy_extent_1"],
                          xy_offset         = config.defaults["xy_offset"]):
    """Transform (lon, lat) coordinates into a south-polar stereographic
    projection (SPSP) defined in cartesian coordinates (x,y). This function
    also requires the data to be plotted (array indices are reordered).

    This function should be used for plotting the south polar region from
    global data. For data that is limited to the south pole already, just use
    the north polar function, as this function is under development
    (i.e., it does weird things).


    Parameters
    ----------
    lon, lat : arrays (nj,ni) of float
        Longitude and latitude coordinates in degrees_east and degrees_north,
        respectively. Must be the same shape.

    data : arrays (nt, nj, ni) or (nj, ni) of float
        Data that is required to be plotted. Can have an extra first dimension
        with respect to lon/lat (e.g., for a 2D time series).


    Optional parameters
    -------------------
    Default values for the following parameters are set in the configuration:

    extent_latitude : float
        The latitude in degrees_south corresponding to the outer limit of the
        plot boundary. This sets the scale of the (x,y) coordinate system.

    central_longitude : float
        Analogous to the parameter of the same name in
        cartopy.crs.SouthPolarStereo().

    The following parameters determine the scale and offset of the (x,y)
    coordinate system. The xy_extent_* parameters should rarely need to be
    changed since the scales of x and y are arbitrary, but for plotting it is
    essential that the coordinate axes match these. Also it is currently
    required that xy_extent_0 = -xy_extent_1 [so that the origin is (0,0)].
    A ValueError is raised if this is not the case.

    xy_extent_0 : 2-tuple of float
        The (x,y) coordinates of the lower-left corner of the plot. The
        standard value is (-1, -1).

    xy_extent_1 : 2-tuple of float
        The (x,y) coordinates of the upper-right corner of the plot. The
        standard value is (1, 1).

    xy_offset : 2-tuple of float
        The (x,y) coordinate shift by which to offset the overall plot area.
        Default is (0, 0), i.e., no offset.

        This parameter may be used to offset the centre of the plot, but note
        that this has no effect on the coordinate projection: it simply
        translates coordinates so that pole appears at coordinates xy_offset
        (an analogy would be making a standard stereographic projection with the
        pole in the centre of the plot, then using the matplotlib interactive
        toolbar to manually click and drag the data around).


    Returns
    -------
    x, y : arrays of float (nj, ni)
        The coordinates in the transformed [and translated if xyp is not at the
        center of xy0 and xy1] coordinate system.

    data_transformed : array of float, same shape as data
        The data which may now be plotted with x, y coordinates.

    """

    if not (    xy_extent_0[0] == -xy_extent_1[0]
            and xy_extent_0[1] == -xy_extent_1[1]):
        raise ValueError("xy_extent_0 must equal -xy_extent_1")

    # Note: I discovered this works by trial-and-error and, as of
    # ----- writing, did not work out exactly how/why it works

    x_n, y_n = lonlat_to_xy_npsp(lon, -lat,
                                 extent_latitude=abs(extent_latitude),
                                 central_longitude=central_longitude + 90.,
                                 xy_extent_0=xy_extent_0,
                                 xy_extent_1=xy_extent_1,
                                 xy_offset=xy_offset)
    # ^ reason for + 90 in central_longitude?

    # I particularly have no idea why the following works (well, it reverses
    # the mirroring of the array reversing, but I don't know why any of this is
    # necessary...
    x = y_n[::-1,:]  # << swapping x and y is necessary (apparently)
    y = x_n[::-1,:]

    if np.ndim(data) == 2:
        data_transformed = data[::-1,:]
    elif np.ndim(data) == 3:
        data_transformed = data[:,::-1,:]
    else:
        raise ValueError(f"Expected 2 or 3 dims (got {np.ndim(data)})")

    return x, y, data_transformed

