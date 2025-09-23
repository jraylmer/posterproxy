"""Provide functions for loading and preparing Natural Earth map data and for
plotting them in maps produced in the polar stereographic projection.
"""

from pathlib import Path
import warnings

import matplotlib as mpl
import numpy as np
import shapefile

from .. import config, transform
from ..utils import gen as gen_utils, mpl as mpl_utils


def _validate_ne_options(scale, bbox, reject_method):
    """Validates input options to Natural Earth data plotting interface (common
    options for loading land/coastlines/etc. If any input is invalid they are
    reverted to a default and a warning is raised.
    """

    valid_scales = ["110m", "50m"]

    valid_reject_methods = ["individual_points", "whole_parts",
                            "whole_geometries"]

    # Allowed bbox min/max:
    bbox_lon_min = 0.
    bbox_lon_max = 360.
    bbox_lat_min = -90.
    bbox_lat_max = 90.

    if scale not in valid_scales:
        if scale == "10m":
            warnings.warn(f"scale '10m' not currently supported; "
                          + f"setting scale = '50m'")
            scale = "50m"
        else:
            warnings.warn(f"scale '{scale}' not valid; "
                          + f"setting scale = '{valid_scales[0]}'")
        scale = valid_scales[0]

    # Ensure bbox is the correct shape:
    if np.shape(bbox) == 4:

        bbox = [bbox[:2], bbox[2:]]

        warnings.warn(f"Received wrong shape bbox: using bbox = [["
                      + f"{bbox[0][0]:.2f}, {bbox[0][1]:.2f}], ["
                      + f"{bbox[1][0]:.2f}, {bbox[1][1]:.2f}]]")

    elif np.shape(bbox) != (2,2):
        bbox = [[bbox_lon_min, bbox_lat_min],
                [bbox_lon_max, bbox_lat_max]]

        warnings.warn(f"Received wrong shape bbox: using default bbox = [["
                      + f"{bbox[0][0]:.2f}, {bbox[0][1]:.2f}], ["
                      + f"{bbox[1][0]:.2f}, {bbox[1][1]:.2f}]]")

    # Ensure input lon/lat bbox values are in the valid range (no need
    # to warn here as anything outside of allowed range is unphysical):
    bbox[0][0] = max(bbox_lon_min, bbox[0][0] % 360.)
    bbox[1][0] = min(bbox_lon_max, bbox[1][0])
    bbox[0][1] = max(bbox_lat_min, bbox[0][1])
    bbox[1][1] = min(bbox_lat_max, bbox[1][1])

    if reject_method not in valid_reject_methods:
        warnings.warn(f"reject_method '{reject_method}' not valid; setting "
                      + f"reject_method = {valid_reject_methods[0]}'")
        reject_method = valid_reject_methods[0]

    return scale, bbox, reject_method


def _ne_data_file(feature, scale):
    """Return the absolute path to data (*.shp and auxilliary files, *.cpg,
    *.dbf, *.prj, *.shx, in the same directory).
    """
    return Path(config.ne_data_dir, f"ne_{scale}_{feature}",
                f"ne_{scale}_{feature}.shp").resolve()


def _validate_lonlats_(x):
    """Used for ensuring data coordinates are valid (e.g., 110m land data has some
    latitudes < -90, for some reason).
    """
    x = np.maximum(x, np.array([[-180., -90.]]))
    x = np.minimum(x, np.array([[180., 90.]]))
    return x


def _get_ne_land_polygon_data(scale, bbox=[(0., -90.), (360., 90.)],
                              reject_method="individual_points"):
    """Load and prepare land polygon data from shapefile, for specified Natural
    Earth scale and bounding box.


    Parameters
    ----------
    scale : str {'110m', '50m', '10m'}
        Scale of data to load, where '110m' means 1:110 million, etc.; i.e.,
        '110m' is coarsest and '10m' is the most detailed.


    Optional parameters
    -------------------
    bbox : length-2 iterable of length-2 iterable of float,
           default = [ (0, -90), (360, 90) ]
        Rectangular (in lon-lat space) bounding box to filter data: the first
        iterable contains the longitude and latitude coordinates, respectively,
        of the south-west corner of the bbox, and the second contains that of
        the north-east corner. Longitudes are in degrees_east (0-360) and
        latitudes are in degrees_north. Default bbox is global.

        Note that all data is loaded initially as it is not possible to know
        in advance which specific geometries are within bbox.

    reject_method : str {'individual_points', 'whole_geometries', 'whole_parts'}
        Method for filtering geometries outside of bbox.

        'individual_points' (default option)
            Each node/vertex is considered individually for all geometries, and
            if it lies outside of bbox, it is removed from the corresponding
            geometry. This means that some geometries are returned incomplete.
            This is usually the desired behaviour for non-global plots. Note
            that polygons are always closed so that any gaps created by this
            method will be connected by a straight line.

        'whole_parts'
            For each part of a given geometry, if any node/vertex of any of its
            parts are outside of bbox, the whole part is removed from that
            geometry. If all parts of a geometry are removed in this way, the
            whole geometry record is removed.

        The following method is not currently implemented and behaves the same as
        to 'whole_parts':

        'whole_geometries'
            For each geometry, if any node/vertex of any of its parts are
            outside of bbox, the whole geometry is removed from the final
            collection of geometries.


    Returns
    -------
    geometries_list : list of list of list of array
        The land coordinate data.

        The data structure is not straightforward/intuitive, but this is a
        'private' function and the user should not normally have to deal with
        it in this form. There is one extra level of list nesting at the outer
        level compared to that for coastline data because polygons can have
        interior as well as exterior paths.

        The outer most list is that of all geometries, g, such that:

        geometries_list[g]       is a list of all shapes, s, belonging to g
        geometries_list[g][s]    is a list of all parts, p, belonging to shape s
        geometries_list[g][s][p] is an array of shape (Ngsp, 2)

        where Ngsp is the number of nodes/vertices describing p, so that
        geometries_list[g][s][p][:,2] are the corresponding longitudes and
        latitudes in degrees_east and degrees_north, respectively.

    """

    # Validate input options (setting to defaults if invalid):
    scale, bbox, reject_method = _validate_ne_options(scale, bbox, reject_method)

    shape_file = _ne_data_file("land", scale)

    geometries_list = []

    with shapefile.Reader(shape_file) as sf:

        # Still don't understand the terminology and it seems inconsistent: we
        # load 'shapes' from the shapefile (makes sense!) but then recast these
        # as 'geometries' using the geo_interface (also makes sense!). However,
        # some 'geometries' consist of multiple 'parts' (shapefiles are supposed
        # to contain only one type of geometry; alas, these natural earth data
        # files contain mixtures...). So I'm going to call the output of the
        # geo_interface 'geometries', each of which consist of one or more
        # 'shapes', each of which contain one or more 'parts'.
        sf_geometries = sf.shapes().__geo_interface__["geometries"]

        # Loop over all geometries, g:
        for g in range(len(sf_geometries)):
            
            # Get lists of lists of coordinate arrays (yes, just go with it)
            # of shape (Ngsp, 2).
            #
            # For all coordinates, also check that they are valid (i.e.,
            # latitudes not less than -90, which, from testing, is present in
            # the datasets presumably as some sort of rounding error).
            #
            # While doing this it is necessary to determine whether a 'Polygon'
            # or 'MultiPolygon' is being loaded, because the datatype needs to
            # be made consistent (i.e., convert all 'Polygon' to
            # 'MultiPolygon' of one shape).
            #
            if sf_geometries[g]["type"] == "Polygon":

                # Wrap an extra list around the data to make it a MultiPolygon:
                coords_g = [[_validate_lonlats_(np.array(dat))
                             for dat in sf_geometries[g]["coordinates"]]]
                # (where dat is the array of coordinates)

            elif sf_geometries[g]["type"] == "MultiPolygon":

                coords_g = [[_validate_lonlats_(np.array(dat)) for dat in part]
                            for part in sf_geometries[s]["coordinates"]]

            # (else do nothing: if there is something *else* in the data then
            # frankly I don't know what to do with it anyway so just ignore it)

            # Before saving coords_g to geometries_list, apply bbox filtering
            #
            # Each shape, s, of g must be checked, and within geometry each
            # part, p, (the exteriors and any intereriors) must also be checked.
            #
            # Loop over s in reverse so that, if required, the corresponding
            # list elements can be deleted without affecting later iterations:
            for s in range(len(coords_g)-1, -1, -1):

                # Loop over s's parts (p; also in reverse):
                for p in range(len(coords_g[s])-1, -1, -1):

                    outside_bbox_gsp = np.logical_or.reduce(
                        (coords_g[s][p][:,0] % 360 < bbox[0][0],
                         coords_g[s][p][:,0] % 360 > bbox[1][0],
                         coords_g[s][p][:,1]       < bbox[0][1],
                         coords_g[s][p][:,1]       > bbox[1][1]),
                        axis=0)
                    # the array outside_bbox_gsp is True where coordinates are
                    # OUTSIDE the bbox (i.e., reject if True)

                    if reject_method == "individual_points":

                        # If a coordinate is outside the bbox, reject that point
                        # alone (leave other vertices in place).
                        #
                        # Note that if this opens up a polygon, later code will
                        # always re-close it (all polygons drawn with matplotlib
                        # are closed by default).
                        #
                        # Here we are indexing the array, so can use boolean
                        # slicing:
                        if any(outside_bbox_gsp):
                            coords_g[s][p] = coords_g[s][p][~outside_bbox_gsp,:]

                    else:
                        # reject_method == 'whole_geometries'
                        # reject_method == 'whole_parts'
                        
                        # If any coordinate is outside bbox, reject whole part:
                        if any(outside_bbox_gsp):
                            del coords_g[s][p]

                    # Remove empty parts lists:
                    if len(coords_g[s][p]) == 0:
                        del coords_g[s][p]

                # Remove empty geometry lists:
                if not coords_g[s]:
                    del coords_g[s]

            # Finished processing this geometry (main outer loop):
            geometries_list.append(coords_g)

    return geometries_list


def _get_patches(coordinates, PatchCollection_kw=None):
    """Generate matplotlib patch collection for plotting land polygons.


    Parameters
    ----------
    coordinates : list of list of list of arrays of shape (Ngsp, 2)
        The coordinate list returned by function _get_ne_land_polygon_data().


    Optional parameters
    -------------------
    PatchCollection_kw : dict or None
        Keyword arguments passed to matplotlib.collection.PatchCollection(),
        e.g. to set facecolor of land polygons. If None, passes an empty
        dictionary.


    Returns
    -------
    matplotlib.collections.PatchCollection instance

    """

    patches = []  # list to append each patch to.

    # Patches are either simple polygons or polygons with holes if the
    # corresponding geometry's shape has multiple parts (in which case the
    # first part is interpreted as the exterior path, and the remaining parts
    # are considered to be the interior paths).

    # Loop over each geometry:
    for g in range(len(coordinates)):
        # Loop over each shape:
        for s in range(len(coordinates[g])):

            if len(coordinates[g][s]) > 1:
                # If this polygon has interiors (i.e., MultiPolygon with
                # multiple parts), need special treatment to set interiors:
                patches.append(mpl_utils.polygon_with_interiors(coordinates[g][s]))

            else:
                # Simple closed polygon without interiors (i.e., only one part
                # which is the exterior path):
                patches.append(mpl.patches.Polygon(coordinates[g][s][0],
                                                   closed=True))

    if PatchCollection_kw is None:
        PatchCollection_kw = {}

    return mpl.collections.PatchCollection(patches, **PatchCollection_kw)


def _get_ne_coastline_data(scale, bbox=[(0., -90.), (360., 90.)],
                           reject_method="individual_points"):
    """Load and prepare coastline path data from shapefile, for specified
    Natural Earth scale and bounding box.


    Parameters
    ----------
    scale : str {'110m', '50m', '10m'}
        Scale of data to load, where '110m' means 1:110 million, etc.; i.e.,
        '110m' is coarsest and '10m' is the most detailed.


    Optional parameters
    -------------------
    bbox : length-2 iterable of length-2 iterable of float,
           default = [ (0, -90), (360, 90) ]
        Rectangular (in lon-lat space) bounding box to filter data: the first
        iterable contains the longitude and latitude coordinates, respectively,
        of the south-west corner of the bbox, and the second contains that of
        the north-east corner. Longitudes are in degrees_east (0-360) and
        latitudes are in degrees_north. Default bbox is global.

        Note that all data is loaded initially as it is not possible to know
        in advance which specific geometries are within bbox.

    reject_method : str {'individual_points', 'whole_geometries', 'whole_parts'}
        Method for filtering geometries outside of bbox.

        'individual_points' (default option)
            Each node/vertex is considered individually for all geometries, and
            if it lies outside of bbox, it is removed from the corresponding
            geometry. This means that some geometries are returned incomplete.
            This is usually the desired behaviour for non-global plots. Note
            that polygons are always closed so that any gaps created by this
            method will be connected by a straight line.

        'whole_parts'
            For each part of a given geometry, if any node/vertex of any of its
            parts are outside of bbox, the whole part is removed from that
            geometry. If all parts of a geometry are removed in this way, the
            whole geometry record is removed.

        The following method is not currently implemented and behaves the same as
        to 'whole_parts':

        'whole_geometries'
            For each geometry, if any node/vertex of any of its parts are
            outside of bbox, the whole geometry is removed from the final
            collection of geometries.


    Returns
    -------
    geometries_list : list of list of array
        The coastline coordinate data.

        The data structure is not straightforward/intuitive, but this is a
        'private' function and the user should not normally have to deal with
        it in this form. There is one less level of list nesting compared to
        that for the land data because path data does not have interiors.

        The outer most list is that of all geometries, g, such that:

        geometries_list[g]    is a list of all paths, p, belonging to g
        geometries_list[g][p] is an array of shape (Ngp, 2)

        where Ngp is the number of nodes/vertices describing p, so that
        geometries_list[g][p][:,2] are the corresponding longitudes and
        latitudes in degrees_east and degrees_north respectively.

    """

    # Validate input options (setting to defaults if invalid):
    scale, bbox, reject_method = _validate_ne_options(scale, bbox, reject_method)

    shape_file_path = _ne_data_file("coastline", scale)

    geometries_list = []

    with shapefile.Reader(shape_file_path) as sf:
        
        # Still don't understand the terminology and it seems inconsistent: we
        # load 'shapes' from the shapefile (makes sense!) but then recast these
        # as 'geometries' using the geo_interface (also makes sense!). However,
        # some 'geometries' consist of multiple 'parts' (shapefiles are supposed
        # to contain only one type of geometry; alas, these natural earth data
        # files contain mixtures...). So I'm going to call the output of the
        # geo_interface 'geometries', each of which consist of (in the case of
        # coastlines) one or more "paths".
        sf_geometries = sf.shapes().__geo_interface__["geometries"]

        # Loop over all geometries, g:
        for g in range(len(sf_geometries)):

            # Get lists of coordinate arrays of shape (Ngp, 2).
            #
            # For all coordinates, also check that they are valid (i.e.,
            # latitudes not less than -90, which, from testing, is present
            # in the datasets presumably as some sort of rounding error).
            #
            # While doing this, it is necessary to determine whether a
            # 'LineString' or 'MultiLineString' is being loaded, because the
            # datatype needs to be made consistent (i.e., convert all
            # 'LineString' to 'MultiLineString' of one path).
            #
            if sf_geometries[g]["type"] == "LineString":

                # Put an extra list around to make this a 'MultiLineString'
                # with one path/part (i.e., same format for all sets of
                # coordinates):
                coords_g = [_validate_lonlats_(
                    np.array(sf_geometries[g]["coordinates"]))]

            elif sf_geometries[g]["type"] == "MultiLineString":
                coords_g = [_validate_lonlats_(np.array(dat))
                            for dat in sf_geometries[g]["coordinates"]]

            # (else do nothing: if there is something *else* in the data then
            # frankly I don't know what to do with it anyway so just ignore it)
            
            # Check coordinates are in the specified bbox.
            #
            # Each path/part, p, of geometry g must be checked. Loop over p in
            # reverse so that, if required, the corresponding list elements can
            # be deleted without consequence for later iterations:
            for p in range(len(coords_g)-1, -1, -1):

                outside_bbox_gp = np.logical_or.reduce(
                    (coords_g[p][:,0] % 360 < bbox[0][0],
                     coords_g[p][:,0] % 360 > bbox[1][0],
                     coords_g[p][:,1]       < bbox[0][1],
                     coords_g[p][:,1]       > bbox[1][1]),
                    axis=0)
                # the array outside_bbox_gp is True where coordinates are
                # OUTSIDE the bbox (i.e., reject if True)

                if reject_method == "individual_points":

                    # If a coordinate is outside the bbox, reject that point
                    # alone (leave other vertices in place).
                    #
                    # Here we are indexing the array, so can use boolean
                    # slicing:
                    if any(outside_bbox_gp):
                        coords_g[p] = coords_g[p][~outside_bbox_gp,:]

                else:
                    # reject_method == 'whole_geometries'
                    # reject_method == 'whole_parts'
                    
                    # If any coordinate is outside bbox, reject whole part:
                    if any(outside_bbox_gp):
                        del coords_g[p]

                # Remove empty parts lists:
                if len(coords_g[p]) == 0:
                    del coords_g[p]

            # Finished processing this geometry (main outer loop):
            geometries_list.append(coords_g)

    return geometries_list


def _get_lines(coordinates, PathCollection_kw=None):
    """Generate matplotlib path collection for plotting coastlines.


    Parameters
    ----------
    coordinates : list of list of arrays of shape (Ngp, 2)
        The coordinate list returned by function _get_ne_coastline_data().


    Optional parameters
    -------------------
    PathCollection_kw : dict or None
        Keyword arguments passed to matplotlib.collection.PathCollection(),
        e.g., to set color of coastlines. If None, passes an empty dictionary.


    Returns
    -------
    matplotlib.collections.PathCollection instance

    """

    paths_list = []

    for g in range(len(coordinates)):  # ------ for each geometry
        for p in range(len(coordinates[g])):  # for each path
            paths_list.append(mpl.patches.Path(coordinates[g][p]))

    if PathCollection_kw is None:
        PathCollection_kw = {}

    return mpl.collections.PathCollection(paths_list, **PathCollection_kw)


def land_overlay(axs, land=True, coastlines=True,
                 mpl_transform=None,
                 south_polar   = False,
                 psp_xy_kw     = config.defaults["psp_xy_kw"],
                 scale         = config.defaults["ne_scale"],
                 bbox          = config.defaults["ne_bbox"],
                 reject_method = config.defaults["ne_reject_method"],
                 coast_line_kw = config.defaults["coast_line_kw"],
                 land_patch_kw = config.defaults["land_patch_kw"]):
    """Apply a Natural Earth data land and/or coastline overlay
    to one or more matplotlib axes.


    Parameters
    ----------
    axs : matplotlib axes.Axes instance, or array of such
        The axes(s) to apply the land/coast overlay to.


    Optional parameters
    -------------------
    land, coastlines : bool (default: True for both)
        Whether to add filled land polygons and/or coastlines respectively.
        Note that coastlines are separate features to land (i.e., they are not
        the outline of land polygons).

    mpl_transform : matplotlib transformation instance or None (default)
        Coordinate transformation passed to both land and coastline patches for
        plotting. This parameter provides a means to load Natural Earth data
        using this interface but still be able to plot it on any map projection.
        By default, this is set to None and not used, instead applying the
        internal polar stereographic projection (via psp_xy_kw).

    south_polar : bool (default:  False)
        If True, assumes that the south polar region is needed (i.e., that
        function spsp_from_global() is being used).

    Default values for the following parameters are set in the configuration:

    psp_xy_kw : dictionary
        Properties of the polar stereographic projection in xy coordinates as
        keyword arguments passed to function psp.lonlat_to_xy_npsp().

    scale : str {'110m', '50m', '10m'}
        Scale of data to load, where '110m' means 1:110 million, etc.; i.e.,
        '110m' is coarsest and '10m' is the most detailed.

    bbox : length-2 iterable of length-2 iterable of float,
        Rectangular (in lon-lat space) bounding box to filter data: the first
        iterable contains the longitude and latitude coordinates, respectively,
        of the south-west corner of the bbox, and the second contains that of
        the north-east corner. Longitudes are in degrees_east (0-360) and
        latitudes are in degrees_north.

        Note that all data is loaded initially as it is not possible to know
        in advance which specific geometries are within bbox.

    reject_method : str {'individual_points', 'whole_geometries', 'whole_parts'}
        Method for filtering geometries outside of bbox.

        'individual_points' (default option)
            Each node/vertex is considered individually for all geometries, and
            if it lies outside of bbox, it is removed from the corresponding
            geometry. This means that some geometries are returned incomplete.
            This is usually the desired behaviour for non-global plots. Note
            that polygons are always closed so that any gaps created by this
            method will be connected by a straight line.

        'whole_parts'
            For each part of a given geometry, if any node/vertex of any of its
            parts are outside of bbox, the whole part is removed from that
            geometry. If all parts of a geometry are removed in this way, the
            whole geometry record is removed.

        The following method is not currently implemented and behaves the same as
        to 'whole_parts':

        'whole_geometries'
            For each geometry, if any node/vertex of any of its parts are
            outside of bbox, the whole geometry is removed from the final
            collection of geometries.

    land_patch_kw = dict
        Properties (color, etc.) of land patches as keyword arguments passed to
        matplotlib.collections.PatchCollection.

    coast_line_kw = dict
        Properties (color, etc.) of coastline lines as keyword arguments passed
        to matplotlib.collections.PathCollection.

    """

    if mpl_transform is None:
        # Prepare for explicit polar stereographic projection
        mod_psp_xy_kw = psp_xy_kw.copy()  # need to modify some

        if south_polar:
            # TODO: de-hardcode this property:
            bbox = [(0., -90.), (360., 45.)]

            # TODO: abstract this into psp module:
            if "central_longitude" in mod_psp_xy_kw.keys():
                mod_psp_xy_kw["central_longitude"] += 90.

        else:
            # For north polar stereo, need a minimum latitude so that 90 S is
            # not included (south pole is projected to a radius of infinity)
            # Convert to list of lists (can't update tuple values):
            bbox = [[bbox[0][0], bbox[0][1]], [bbox[1][0], bbox[1][1]]]
            bbox[0][1] = max(config.defaults["ne_bbox_min_lat"], bbox[0][1])
            bbox[1][1] = max(config.defaults["ne_bbox_min_lat"], bbox[1][1])

    # Common arguments passed to functions loading Natural Earth coordinates:
    nef_kw = {"scale": scale, "bbox": bbox, "reject_method": reject_method}

    if land:

        land_coords = _get_ne_land_polygon_data(**nef_kw)

        if mpl_transform is None:
            
            # Explicitly transform coordinates to polar stereographic:
            if south_polar:
                for g in range(len(land_coords)):
                    for s in range(len(land_coords[g])):
                        for p in range(len(land_coords[g][s])):
                            x_gsp, y_gsp = transform.lonlat_to_xy_npsp(
                                 land_coords[g][s][p][:,0],  # longitudes
                                -land_coords[g][s][p][:,1],  # latitudes
                                **mod_psp_xy_kw)

                            # Re-shape into required format (Npoints, 2):
                            land_coords[g][s][p] = np.vstack((y_gsp, x_gsp)).T

            else:  # north polar
                for g in range(len(land_coords)):
                    for s in range(len(land_coords[g])):
                        for p in range(len(land_coords[g][s])):
                            x_gsp, y_gsp = transform.lonlat_to_xy_npsp(
                                land_coords[g][s][p][:,0],  # longitudes
                                land_coords[g][s][p][:,1],  # latitudes
                                **mod_psp_xy_kw)

                            # Re-shape into required format (Npoints, 2):
                            land_coords[g][s][p] = np.vstack((x_gsp, y_gsp)).T

        # Create patches:
        patch_kw = gen_utils.set_default_kw(land_patch_kw,
                                            config.defaults["land_patch_kw"])
        land_patches = _get_patches(land_coords, PatchCollection_kw=patch_kw)

        if mpl_transform is not None:
            # Using implicit transform, so apply to patches:
            land_patches.set_transform(mpl_transform)

        mpl_utils.assign_collection(axs, land_patches)

    if coastlines:

        coastline_coords = _get_ne_coastline_data(**nef_kw)

        if mpl_transform is None:
            
            # Explicitly transform coordinates to polar stereographic:
            if south_polar:
                for g in range(len(coastline_coords)):
                    for p in range(len(coastline_coords[g])):
                        x_gp, y_gp = transform.lonlat_to_xy_npsp(
                             coastline_coords[g][p][:,0],  # longitudes
                            -coastline_coords[g][p][:,1],  # latitudes
                            **mod_psp_xy_kw)

                        # Re-shape into required format (Npoints, 2):
                        coastline_coords[g][p] = np.vstack((y_gp, x_gp)).T

            else:  # north polar

                # Transform to north polar stereographic projection:
                for g in range(len(coastline_coords)):
                    for p in range(len(coastline_coords[g])):
                        x_gp, y_gp = transform.lonlat_to_xy_npsp(
                            coastline_coords[g][p][:,0],  # longitudes
                            coastline_coords[g][p][:,1],  # latitudes
                            **mod_psp_xy_kw)

                        # Re-shape into required format (Npoints, 2):
                        coastline_coords[g][p] = np.vstack((x_gp, y_gp)).T

        # Force coastlines to never be filled:
        coast_line_kw["facecolors"] = "none"

        # Create Paths:
        line_kw = gen_utils.set_default_kw(coast_line_kw,
                                           config.defaults["coast_line_kw"])
        coastline_paths = _get_lines(coastline_coords, PathCollection_kw=line_kw)

        if mpl_transform is not None:
            # Using implicit transform, so apply to paths:
            coastline_paths.set_transform(mpl_transform)

        mpl_utils.assign_collection(axs, coastline_paths)


def xy_land_overlay(axs, land=True, coastlines=True, south_polar=False,
                    psp_xy_kw     = config.defaults["psp_xy_kw"],
                    scale         = config.defaults["ne_scale"],
                    bbox          = config.defaults["ne_bbox"],
                    reject_method = config.defaults["ne_reject_method"],
                    coast_line_kw = config.defaults["coast_line_kw"],
                    land_patch_kw = config.defaults["land_patch_kw"]):
    """"Apply a Natural Earth data land and/or coastline overlay
    in polar stereographic projection in cartesian (x,y) coordinates, to
    one or more matplotlib axes.


    Parameters
    ----------
    axs : matplotlib axes.Axes instance, or array of such
        The axes(s) to apply the land/coast overlay to (which should have,
        or have had already, a call to posterproxy.prepare_axes() before it
        makes sense to add land/coastlines with this function.


    Optional parameters
    -------------------
    land, coastlines : bool (default: True for both)
        Whether to add filled land polygons and/or coastlines respectively.
        Note that coastlines are separate features to land (i.e., they are not
        the outline of land polygons).

    south_polar : bool (default:  False)
        If True, assumes that the south polar region is needed (i.e., that
        function spsp_from_global() is being used).

    Default values for the following parameters are set in the configuration:

    psp_xy_kw : dictionary
        Properties of the polar stereographic projection in xy coordinates as
        keyword arguments passed to function psp.lonlat_to_xy_npsp().

    scale : str {'110m', '50m', '10m'}
        Scale of data to load, where '110m' means 1:110 million, etc.; i.e.,
        '110m' is coarsest and '10m' is the most detailed.

    bbox : length-2 iterable of length-2 iterable of float,
        Rectangular (in lon-lat space) bounding box to filter data: the first
        iterable contains the longitude and latitude coordinates, respectively,
        of the south-west corner of the bbox, and the second contains that of
        the north-east corner. Longitudes are in degrees_east (0-360) and
        latitudes are in degrees_north.

        Note that all data is loaded initially as it is not possible to know
        in advance which specific geometries are within bbox.

    reject_method : str {'individual_points', 'whole_geometries', 'whole_parts'}
        Method for filtering geometries outside of bbox.

        'individual_points' (default option)
            Each node/vertex is considered individually for all geometries, and
            if it lies outside of bbox, it is removed from the corresponding
            geometry. This means that some geometries are returned incomplete.
            This is usually the desired behaviour for non-global plots. Note
            that polygons are always closed so that any gaps created by this
            method will be connected by a straight line.

        'whole_parts'
            For each part of a given geometry, if any node/vertex of any of its
            parts are outside of bbox, the whole part is removed from that
            geometry. If all parts of a geometry are removed in this way, the
            whole geometry record is removed.

        The following method is not currently implemented and behaves the same as
        to 'whole_parts':

        'whole_geometries'
            For each geometry, if any node/vertex of any of its parts are
            outside of bbox, the whole geometry is removed from the final
            collection of geometries.

    land_patch_kw = dict
        Properties (color, etc.) of land patches as keyword arguments passed to
        matplotlib.collections.PatchCollection.

    coast_line_kw = dict
        Properties (color, etc.) of coastline lines as keyword arguments passed
        to matplotlib.collections.PathCollection.

    """
    land_overlay(axs, land=land, coastlines=coastlines, mpl_transform=None,
                 south_polar=south_polar, psp_xy_kw=psp_xy_kw, scale=scale,
                 bbox=bbox, reject_method=reject_method,
                 coast_line_kw=coast_line_kw, land_patch_kw=land_patch_kw)

