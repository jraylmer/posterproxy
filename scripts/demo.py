"""Minimal work example comparing use of this package against cartopy."""

import matplotlib as mpl
import matplotlib.pyplot as plt

import cartopy
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature, COLORS

import posterproxy as psp


if __name__ == "__main__":

    scale    = "110m"
    latitude = 60.

    # Create figure:
    fig = plt.figure(figsize=(6.2,3.4))

    # ----------------------------------------------------
    # Create a north polar stereographic map using cartopy
    # ----------------------------------------------------
    
    # Create geoAxes instance, positioned on the left-hand side of the figure:
    ax1 = fig.add_axes([.05, .05, .425, .9], projection=ccrs.NorthPolarStereo())
    
    # Set the map limits:
    ax1.set_extent((0., 360., latitude, 90.), crs=ccrs.PlateCarree())

    # Plotting data would go here, e.g.:
    #ax1.contourf(lon, lat, data, transform=ccrs.PlateCarree())

    # Add land and coastlines:
    ax1.add_feature(NaturalEarthFeature("physical", "land", scale,
                                        facecolor=COLORS["land"],
                                        edgecolor="none"))

    ax1.add_feature(NaturalEarthFeature("physical", "coastline", scale,
                                        linewidth=mpl.rcParams["grid.linewidth"],
                                        facecolor="none", edgecolor="k"))

    # Add gridlines (to configure and/or change appearance, we would keep a
    # reference to the return value below and use gl.set_xlocator(), etc.):
    ax1.gridlines()

    # Label this subplot with the module and version:
    ax1.set_title(f"{cartopy.__name__} v{cartopy.__version__}")


    # --------------------------------------------------------
    # Create a north polar stereographic map using posterproxy
    # --------------------------------------------------------
    
    # Specify the projection parameters:
    psp_xy_kw = {"extent_latitude": latitude}

    # Create regular/cartesian axes instance on the right-hand side of figure:
    ax2 = fig.add_axes([.525, .05, .425, .9])

    # Prepare the axes for appropriately displaying polar stereographic data
    # (this sets the axes limits, aspect ratio, and removes ticks/tick labels):
    psp.prepare_axes(ax2)

    # Ploting data would go here: generate the (lon,lat) coordinates in (x,y)
    # polar steoreographic system then use pcolormesh or other plotting
    # commands as usual (without 'transform' keyword argument):
    #
    #x, y = psp.lonlat_to_xy_npsp(lon_bnds, lat_bnds, **psp_xy_kw)
    #ax2.contourf(x, y, data)

    # Add land, coastlines, and grid lines:
    psp.xy_land_overlay(ax2, scale=scale, psp_xy_kw=psp_xy_kw)
    psp.xy_gridlines(ax2, psp_xy_kw=psp_xy_kw)

    # Label this subplot with the module and version:
    ax2.set_title(f"{psp.__name__} v{psp.__version__}")

    # --------------------------------------------------------

    # Finish figure:
    fig.text(.5, .5*ax1.get_position().y0,
             f"Generated with {mpl.__name__} v{mpl.__version__}",
             ha="center", va="center")

    plt.show()
