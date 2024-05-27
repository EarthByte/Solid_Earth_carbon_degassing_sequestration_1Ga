#!/usr/bin/env python3
import math
import sys
import time

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft
import gplately
from netCDF4 import Dataset
from scipy import ndimage

# use gplately conda env to run this script
# conda create -n gplately-env gplately


def compute_sigma_from_km(km, cell_width=0.1, lat=None, earth_radius=6378.0):
    """convert distance in kms to sigma

    :param km: distance in kms
    :param cell_width: the grid cell width in degrees
    :param lat: the latitude in degress. cos(lat)
    :param earth_radius: earth radius in kms

    :returns: the sigma to be used in ndimage.fourier_gaussian

    """
    # assume kernel width = 2*ceiling(sigma)+1

    cell_width_in_kms = math.radians(cell_width) * earth_radius
    if lat is not None:
        cell_width_in_kms = cell_width_in_kms * math.cos(math.radians(lat))

    if cell_width_in_kms > 0:
        return km / cell_width_in_kms / 2 - 1
    else:
        return None


def fft_gaussian_filter(data, distance_km=200, extent=[-180, 180, -90, 90]):
    """fft gaussian filter

    :param data: input 2D array numpy array
    :param distance_km: gaussian kernel width in kms
    :param extent: the extent of data

    :returns: None. the output data will be placed in data(in place)
    """
    lat_n = data.shape[0]
    lon_n = data.shape[1]
    lat_cell_width = (extent[3] - extent[2]) / lat_n
    lon_cell_width = (extent[1] - extent[0]) / lon_n

    # filter horizontally
    for idx, row in enumerate(data):
        nan_mask = np.isnan(row)
        if nan_mask.any():
            row[nan_mask] = np.interp(
                np.flatnonzero(nan_mask), np.flatnonzero(~nan_mask), row[~nan_mask]
            )
        input_ = numpy.fft.fft(row)
        sigma = compute_sigma_from_km(
            distance_km, cell_width=lon_cell_width, lat=lat_cell_width * idx + extent[2]
        )
        if sigma > lon_n / 2 or sigma is None:
            sigma = (lon_n - 1) / 2
        output_ = ndimage.fourier_gaussian(input_, sigma=sigma)
        output_ = numpy.fft.ifft(output_)
        data[idx][:] = output_.real

    # filter vertically
    sigma = compute_sigma_from_km(distance_km, cell_width=lat_cell_width)
    for i in range(lon_n):
        input_ = numpy.fft.fft(data[:, i])
        output_ = ndimage.fourier_gaussian(input_, sigma=sigma)
        output_ = numpy.fft.ifft(output_)
        data[:, i] = output_.real

    return data


def plot(data, result):
    fig, (ax1, ax2) = plt.subplots(
        nrows=1,
        ncols=2,
        subplot_kw={"projection": ccrs.PlateCarree(), "frameon": False},
        figsize=(32, 16),
        dpi=540,
    )

    plt.gray()  # show the filtered result in grayscale

    ax1.imshow(
        data,
        origin="lower",
        transform=ccrs.PlateCarree(),
        extent=[-180, 180, -90, 90],
    )
    ax2.imshow(
        result,
        origin="lower",
        transform=ccrs.PlateCarree(),
        extent=[-180, 180, -90, 90],
    )

    # save the figure without frame so that the image can be used to project onto a globe
    fig.savefig(
        f"test.pdf",
        # bbox_inches="tight",
        # pad_inches=0,
        dpi=120,
        transparent=True,
    )

    plt.close(fig)


def save_netcdf(data, out_filename, extent=[-180, 180, -90, 90], var_name="z"):
    nc_dataset = Dataset(out_filename, "w", "NETCDF4")
    nc_dataset.createDimension("lat", data.shape[0])
    nc_dataset.createDimension("lon", data.shape[1])

    lat_var = nc_dataset.createVariable("lat", "float32", ("lat"))
    lat_var[:] = np.linspace(extent[2], extent[3], num=data.shape[0])
    lon_var = nc_dataset.createVariable("lon", "float32", ("lon"))
    lon_var[:] = np.linspace(extent[0], extent[1], num=data.shape[1])
    dem_bathy_var = nc_dataset.createVariable(var_name, "float32", ("lat", "lon"))
    # dem_bathy_var[:] = np.transpose(dem_bathy)
    dem_bathy_var[:] = data
    nc_dataset.close()


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(sys.argv)
        print("usage: fft_gaussian_filter.py input_file distance_km output_file")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[3]
    distance_km = float(sys.argv[2])
    grid = Dataset(input_file)
    if "z" not in grid.variables.keys():
        print(grid.variables.keys())
        raise Exception(
            "variable name z not found. You need to find out what the name is and change it accordingly in next line."
        )
    data = np.asarray(grid.variables["z"])

    buf = np.copy(data)
    mask = np.isnan(data)

    #start_time = time.time()

    fft_gaussian_filter(buf, distance_km=distance_km)

    #end_time = time.time()

    #print(f"The time of fft_gaussian_filter is :{end_time-start_time}")

    # keep the NaNs
    buf[mask] = np.nan

    # plot for debug
    # plot(data, buf)

    # save filtered grid data to file
    gplately.grids.write_netcdf_grid(buf, output_file)
