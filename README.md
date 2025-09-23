# posterproxy
An explicit implementation of the **po**lar **ster**eographic map **pro**jection in cartesian (**_x_**,&nbsp;**_y_**) coordinates, for plotting data in Python with matplotlib.

Instead of implicitly transforming coordinates via arguments to matplotlib functions, the transformed coordinates are calculated for direct use in plotting commands on regular cartesian axes. The package also includes an interface to add decorators such as grid lines and land/coastlines, the latter using [Natural Earth](https://www.naturalearthdata.com/) data.

This code was originally motivated by some issues encountered using implicit transformations, such as artefacts arising across the 0&#x00B0; meridian and at the pole. These do not occur in the explicit (_x_,&nbsp;_y_) coordinate system because the pole is just the origin (0,&nbsp;0) and there are no discontinuities.

## Required Python packages
* matplotlib
* NumPy
* [Python Shapefile Library (PyShp)](https://pypi.org/project/pyshp/)
* [requests](https://pypi.org/project/requests/) (only needed for the script downloading Natural Earth data)
 
## Limitations
These may be addressed in the future:
* Limited API support for the south polar stereographic projection
* Latitudes are assumed to be geocentric, not geodetic
* Land/coastline overlays only work with 1:110m and 1:50m scale [Natural Earth](https://www.naturalearthdata.com/) data
* Cursor coordinates in the interactive plot window display (_x_,&nbsp;_y_) rather than (_&#x03BB;_,&nbsp;_&#x03D5;_)

## Installation
1. Clone the repository
2. Make the package visible on the python PATH (e.g., on UNIX, `export PYTHONPATH=${PYTHONPATH}:$(pwd)/lib`, executed from the repository top-level directory)
3. Download the Natural Earth data (`scripts/ne_download.py`)

## Basic usage
Assuming data in an array `data` defined on coordinates `lon` and `lat`:

```python
import matplotlib.pyplot as plt
import posterproxy as psp

fig, ax = plt.subplots()
psp.prepare_axes(ax)

x, y = psp.lonlat_to_xy_npsp(lon, lat)
ax.contour(x, y, data)

psp.xy_land_overlay(ax)
psp.xy_gridlines(ax)
```

The equivalent using cartopy:

```python
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

ax = plt.axes(projection=ccrs.NorthPolarStereo())
ax.set_extent((0, 360, 60, 90), crs=ccrs.PlateCarree())

ax.contour(lon, lat, data, transform=ccrs.PlateCarree())

ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
ax.gridlines()
```
The default land, coastlines, and grid lines are set to resemble cartopy's defaults (but can also be changed via configuration files or on-the-fly using function arguments):

![Comparison of north polar stereographic map projection produced in cartopy (left) and posterproxy (right)](/misc/demo.png?raw=true)
