import xarray as xr

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy

def add_coords(var: xr.DataArray, rename: bool = False) -> xr.DataArray:
    """
    Add coordinate information to an xarray DataArray by averaging longitude and latitude over 
    specific dimensions and assigning these as coordinates. Optionally, rename the coordinates.

    Parameters:
    var (xr.DataArray): The input DataArray that contains 'XLONG' and 'XLAT' coordinates.
    rename (bool, optional): If True, rename 'west_east' to 'longitudes' and 'south_north' 
                             to 'latitudes'. Default is False.

    Returns:
    xr.DataArray: The DataArray with updated coordinates. If rename is True, coordinates are also renamed.

    Example:
    --------
    import xarray as xr
    import numpy as np

    # Create a sample DataArray with dummy dimensions and coordinates
    data = np.random.rand(10, 10)
    lon = np.linspace(-180, 180, 10)
    lat = np.linspace(-90, 90, 10)
    ds = xr.DataArray(data, coords=[('south_north', lat), ('west_east', lon)],
                      name='temperature', attrs={'XLONG': lon, 'XLAT': lat})

    # Add coordinates and rename
    ds = add_coords(ds, rename=True)
    print(ds)
    """
    lon_1d = var.XLONG.mean(dim="south_north").values
    lat_1d = var.XLAT.mean(dim="west_east").values
    updated_var = var.drop(["XLONG", "XLAT"]).assign_coords(
        south_north=lat_1d, west_east=lon_1d
    )
    
    if rename:
        updated_var = renamelatlon(updated_var)
    
    return updated_var

def renamelatlon(var: xr.DataArray) -> xr.DataArray:
    """
    Rename coordinates in an xarray DataArray from 'west_east' to 'longitudes' 
    and from 'south_north' to 'latitudes'.

    Parameters:
    var (xr.DataArray): The input DataArray with 'west_east' and 'south_north' coordinates.

    Returns:
    xr.DataArray: The DataArray with renamed coordinates.

    Example:
    --------
    import xarray as xr
    import numpy as np

    # Create a sample DataArray
    data = np.random.rand(10, 10)
    lon = np.linspace(-180, 180, 10)
    lat = np.linspace(-90, 90, 10)
    ds = xr.DataArray(data, coords=[('latitudes', lat), ('longitudes', lon)],
                      name='temperature')

    # Rename coordinates
    ds = renamelatlon(ds)
    print(ds)
    """
    return var.rename({"west_east": "longitudes", "south_north": "latitudes"})


def plot_coast(axes: cartopy.mpl.geoaxes.GeoAxes, color='black', linewidth=2, gridlines_alpha=0.5, states=False) -> None:
    """
    Add coastlines, country borders, and optional state/provincial borders to a Cartopy GeoAxes.

    Parameters:
    axes (cartopy.mpl.geoaxes.GeoAxes): The GeoAxes instance to plot on.
    color (str, optional): Color of the coastlines and borders. Default is 'black'.
    linewidth (int or float, optional): Width of the coastlines and borders. Default is 2.
    gridlines_alpha (float, optional): Transparency level of the gridlines. Default is 0.5.
    states (bool, optional): If True, include state/provincial borders. Default is False.

    Returns:
    gl (cartopy.mpl.gridliner.Gridliner): The gridliner instance with longitude and latitude formatting.

    Example:
    --------
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 5), subplot_kw={'projection': ccrs.PlateCarree()})
    plot_coast(ax, color='blue', linewidth=1.5, gridlines_alpha=0.7, states=True)
    plt.show()
    """

    countries = cfeature.NaturalEarthFeature(
        scale="10m", category="cultural", name="admin_0_countries", facecolor="none"
    )
    axes.add_feature(countries, edgecolor=color, linewidth=linewidth)

    if states:
        states = cfeature.NaturalEarthFeature(
            scale="10m",
            category="cultural",
            name="admin_1_states_provinces_lines",
            facecolor="none",
        )
        axes.add_feature(states, edgecolor=color, linewidth=linewidth)
    
    gl = axes.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=1,
        color="gray",
        alpha=gridlines_alpha,
        linestyle="--",
    )
    gl.top_labels = False
    gl.bottom_labels = True
    gl.left_labels = True
    gl.right_labels = False
    gl.xlines = True
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    return gl

