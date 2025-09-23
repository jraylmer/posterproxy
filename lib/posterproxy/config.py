"""Store default parameters and provide a function to parse configuration files
and update when the main package is imported.
"""

from configparser import ConfigParser
from pathlib import Path
import os

from matplotlib import rcParams as mplrc


# Define default parameters to be accessible from dictionary called 'defaults'
# Set initial values here, update from configuration files upon module import:
#
defaults = {
    "central_longitude"       : 0.,
    "extent_latitude"         : 60.,
    "xy_extent_0"             : [-1., -1.],
    "xy_extent_1"             : [ 1.,  1.],
    "xy_offset"               : [ 0.,  0.],
    "grid_longitudes"         : [float(j) for j in range( 0, 360, 60)],
    "grid_latitudes"          : [float(j) for j in range(50,  90, 10)],
    "grid_longitudes_lat_lims": [50., 89.75],
    "grid_latitudes_lon_lims" : [0., 360.],
    "ne_scale"                : "110m",
    "ne_bbox"                 : [[0., -60.], [360., 90.]],
    "ne_bbox_min_lat"         : -60.,
    "ne_reject_method"        : "individual_points",
    "ne_land"                 : True,
    "ne_coastlines"           : True
}

# Add default input parameters to main function psp.lonlat_to_xy_npsp():
defaults["psp_xy_kw"] = {
    "central_longitude"       : defaults["central_longitude"],
    "xy_extent_0"             : defaults["xy_extent_0"],
    "xy_extent_1"             : defaults["xy_extent_1"],
    "xy_offset"               : defaults["xy_offset"]}

# Add default keyword arguments to matplotlib line/patches. Note zorder set
# so that these elements are above contours (default zorder=1) and below lines
# (default zorder=2).
#
# Note also: land and coastlines are drawn as path/patch collections, so use
# 'facecolors', not 'facecolor' (etc.):
#
defaults["grid_line_kw"] = {
    "color"                   : mplrc["grid.color"],
    "linewidth"               : mplrc["grid.linewidth"],
    "linestyle"               : mplrc["grid.linestyle"],
    "zorder"                  : 1.7}

defaults["land_patch_kw"] = {
    "facecolors"              : "#EFEFDB",
    "edgecolors"              : "none",
    "zorder"                  : 1.5}

defaults["coast_line_kw"] = {
    "edgecolors"              : "k",
    "linewidth"               : mplrc["grid.linewidth"],
    "zorder"                  : 1.6}

ne_data_dir = Path(os.path.dirname(__file__), "..", "..", "dat").resolve()


def _parse_config():
    """Read configuration files and update default parameters. Also creates
    a local configuration file if one does not exist.
    """

    # Paths for global/default and local configuration INI files:
    dir_cfgs = Path(os.path.dirname(__file__), "..", "..", "cfg").resolve()
    cfg_file_global = Path(dir_cfgs, "CONFIG.ini")
    cfg_file_user  = Path(dir_cfgs, "CONFIG_USER.ini")

    # Check that the global configuration file exists:
    if not os.path.isfile(cfg_file_global):
        raise FileNotFoundError("Global configuration file not found (may need"
                                + " to pull or restore package repository)")

    # Create the user configuration file if it does not exist:
    if not os.path.isfile(cfg_file_user):
        with open(cfg_file_global, "r") as fg, open(cfg_file_user, "w") as fu:
            # Write new header and copy from line 6 in global configuration
            # file (removes the "Global configuration part):
            fu.writelines(["# User configuration\n"] + fg.readlines()[5:])

        if os.path.isfile(cfg_file_user):
            print(f"{__name__}: created user configuration file: "
                  + str(cfg_file_user))
        else:
            warnings.warn(f"{__name__}: unable to create user config. file")

    # Read the configuration files (any option in the same 'section' of the
    # local config takes precedence over that in the global config):
    cprsr = ConfigParser()
    cprsr.read([cfg_file_global, cfg_file_user])

    # Define a series of helper functions to parse options from cprsr,
    # carry out any required additional parsing, then assign to correct key in
    # the global defaults dictionary:
    #
    global defaults

    def _parse_string(section, option, defaults_key):
        """Parameters that are meant to be a string."""
        if cprsr.has_option(section, option):
            defaults[defaults_key] = cprsr.get(section, option)

    def _parse_float(section, option, defaults_key):
        """Parameters that are meant to be a single float value."""
        if cprsr.has_option(section, option):
            defaults[defaults_key] = cprsr.getfloat(section, option)

    def _parse_bool(section, option, defaults_key):
        """Parameters that are meant to be a single boolean value."""
        if cprsr.has_option(section, option):
            defaults[defaults_key] = cprsr.getboolean(section, option)

    def _parse_list(section, option, defaults_key, fmt_func=lambda x: float(x),
                    maxlist=None, delimiter=","):
        """For parameters that are a comma-separated list of values."""
        if cprsr.has_option(section, option):
            vals = [fmt_func(i)
                    for i in cprsr.get(section, option).split(delimiter)]
            if maxlist is not None:
                vals = vals[:min(len(vals),maxlist)]
            defaults[defaults_key] = vals

    def _parse_bbox(section, option, defaults_key, delimiter=","):
        """Parse a bounding box (bbox), which in the config file is a list of
        four values, but internally needs to be a list of two length-2 lists.
        """
        if cprsr.has_option(section, option):
            vals = [float(i)
                    for i in cprsr.get(section, option).split(delimiter)][:4]
            defaults[defaults_key] = [ [vals[0], vals[1]], [vals[2], vals[3]] ]

    def _parse_mpl_param(section, option, defaults_key, mplrc_key, fmt_func):
        """For a parameter which is intended as a matplotlib plot style parameter
        (color, linewidth, etc.). These may have the special value 'mplrc' in the
        config files, in which case it is set to a matplotlib rcParam default.
        """
        if cprsr.has_option(section, option):
            cfg_option = cprsr.get(section, option)
            if cfg_option == "mplrc" and mplrc_key is not None:
                defaults[defaults_key][option] = mplrc[f"{mplrc_key}"]
            else:
                defaults[defaults_key][option] = fmt_func(cfg_option)

    s = "MAP DEFAULTS"
    _parse_float(s, "central_longitude", "central_longitude")
    _parse_float(s, "extent_latitude"  , "extent_latitude")

    s = "XY DOMAIN PARAMETERS"
    _parse_list(s, "xy_extent_0", "xy_extent_0", lambda x: float(x), maxlist=2)
    _parse_list(s, "xy_extent_1", "xy_extent_1", lambda x: float(x), maxlist=2)
    _parse_list(s, "xy_offset"  , "xy_offset"  , lambda x: float(x), maxlist=2)

    s = "GRID LINES"
    _parse_list(s, "longitudes"         , "grid_longitudes")
    _parse_list(s, "latitudes"          , "grid_latitudes")
    _parse_list(s, "longitudes_lat_lims", "grid_longitudes_lat_lims")
    _parse_list(s, "latitudes_lon_lims" , "grid_latitudes_lon_lims")

    s = "GRID LINES"
    d = "grid_line_kw"
    _parse_mpl_param(s, "color"    , d, "grid.color"    , lambda x: x)
    _parse_mpl_param(s, "linestyle", d, "grid.linestyle", lambda x: x)
    _parse_mpl_param(s, "linewidth", d, "grid.linewidth", lambda x: float(x))
    _parse_mpl_param(s, "alpha"    , d, "grid.alpha"    , lambda x: float(x))
    _parse_mpl_param(s, "zorder"   , d, None            , lambda x: float(x))

    s = "NATURAL EARTH DATA"
    _parse_float( s, "bbox_min_lat" , "ne_bbox_min_lat")
    _parse_string(s, "reject_method", "ne_reject_method")
    _parse_string(s, "scale"        , "ne_scale")
    _parse_bool(  s, "land"         , "ne_land")
    _parse_bool(  s, "coastlines"   , "ne_coastlines")
    _parse_bbox(  s, "bbox"         , "ne_bbox")

    s = "LAND PATCH PROPS"
    d = "land_patch_kw"
    _parse_mpl_param(s, "edgecolors", d, "patch.edgecolor", lambda x: x)
    _parse_mpl_param(s, "facecolors", d, "patch.facecolor", lambda x: x)
    _parse_mpl_param(s, "zorder"    , d, None             , lambda x: float(x))

    s = "COASTLINE PROPS"
    d = "coast_line_kw"
    _parse_mpl_param(s, "edgecolors", d, "lines.color"   , lambda x: x)
    _parse_mpl_param(s, "linewidths", d, "grid.linewidth", lambda x: float(x))
    _parse_mpl_param(s, "zorder"    , d, None            , lambda x: float(x))

    # Special case: path to Natural Earth data (not part of defaults dictionary):
    if cprsr.has_option("NATURAL EARTH DATA", "directory"):
        global ne_data_dir
        c_ne_data_dir = cprsr.get("NATURAL EARTH DATA", "directory")

        if c_ne_data_dir.lower() != "default":
            ne_data_dir = Path(c_ne_data_dir).resolve()

    # Check main directory exists regardless of whether it was just updated:
    if not os.path.isdir(ne_data_dir):
        raise ValueError("Natural Earth data directory does not exist: "
                         + str(ne_data_dir))

