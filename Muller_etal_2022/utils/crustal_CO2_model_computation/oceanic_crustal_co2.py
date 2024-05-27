# Dependencies:

# netCDF4
# numpy
# scipy

# Reference:

# Müller, R.D. and Dutkiewicz, A., 2018, Oceanic crustal carbon cycle drives 26 million-year atmospheric carbon dioxide periodicities, Science Advances, 4:eaaq0500, 1-7.

# Scotese, C.R., Song, H., Mills, B.J. and van der Meer, D.G., 2021,  Phanerozoic paleotemperatures: The earth’s changing climate during the last 540 million years. Earth-Science Reviews, 215, p.103503.

# Steps:
# 1. Load a text file with a time series of age versus global bottom water temperatures (bwt) (age_deep-ocean-temp_Scotese2021.txt)
# 2. Load an agegrid from a set of paleoage-grids
# 3. Interpolate file from (1) at 1my intervals using linear interpolation
# 4. The bottom water temperature (bwt) that needs to be used for any given age grid cell is the average water temp for the following 20 my (ie the 20 my after the formation of the crust in that cell), simply create a new bwt time series that contains age versus meanbwt20 = mean(bwt) between age and age-20
# 5. For each grid cell in the age grid, compute crustal CO2 by looking up Co2 as a function of meanbwt20 and age using age_bwt_co2_model_bilinear_log.nc as a lookup table.
# 6. Use the collection of computed CO2 values to create a new netcdf file(“upper_crustal_CO2_$age.nc”) in the same dimensions as the input age grid file
# 7. Repeat for all other netcdf age grids.
# 8. Make a summary plot: Reconstruction age versus CO2

import csv
import datetime
import os
import statistics
import pandas as pd
from multiprocessing import Pool
from pathlib import Path

import netCDF4
import numpy as np
from scipy.spatial import KDTree

import gplately


# *********************************************************************************

def time_averaged_bottom_water_temperature(interp_temps, average_time_range):
    """ Calculate the bottom water temperature averaged over a span of `average_time_range`
    million years.
    Looping through time, if the current time X < `average_time_range` Ma, the mean 
    bottom-water temperature between 0-`average_time_range` Ma (inclusive) is taken.
    """
    meanbwt = []
    for idx, x in enumerate(interp_temps):
        if idx < average_time_range:
            meanbwt.append(statistics.mean(interp_temps[0 : idx + 1]))
        else:
            meanbwt.append(statistics.mean(interp_temps[idx - (average_time_range-1) : idx + 1]))
    return meanbwt


def create_co2_grid(
    reconstruction_time,  
    age_temp_series_path, 
    agegrid_path, 
    co2_water_temperature_lookup_table_path, 
    output_path,
    mode,
    average_time_range = 20
    ):
    """Calculate crustal carbon (t/m^2) as a function of time (Myr) and mean bottom water 
    temperature (degree Celsius) over an `average_time_range` from a bilinear-log
    lookup table.

    The time-bottom water temperature series used for the table look-up is obtained from an 
    external file with the full path `age_temp_series_path`. 

    Carbon values are interpolated by age onto the original age grid's domain points, 
    forming a crustal carbon snapshot at `reconstruction_time`. The new grid is exported to 
    a netCDF format using gplately. 
    """

    # Create output path if it does not exist already
    Path(f"{output_path}/{mode}/{reconstruction_time}").mkdir(parents=True, exist_ok=True)

    # ====================== PREPARE AGE GRID ==========================================
    # Read the netCDF age grid, collect all valid times 
    age_grid_nc = netCDF4.Dataset(agegrid_path)

    # Identify netCDF variable names 
    age_grid_z = "z"  # the age grid netcdf variable name (age)

    age_grid_shape = age_grid_nc[age_grid_z].shape
    age_ravel = age_grid_nc[age_grid_z][:].ravel()
    ages_grid = age_ravel[age_ravel.mask == False]
    age_grid_index = np.where(age_ravel.mask == False)[0]
    ages_grid = np.clip(ages_grid, 0, 175)


    # ====================== PREPARE THE CARBON = F(TEMP, AGE) LOOKUP TABLE =========================
    age_temp_co2_nc = netCDF4.Dataset(co2_water_temperature_lookup_table_path)

    # netCDF4 variable names
    age_bwt_co2_grid_x = "x"  # age_bwt_co2_model_bilinear_log.nc variable name (age)
    age_bwt_co2_grid_y = "y"  # age_bwt_co2_model_bilinear_log.nc variable name (below water temperature)
    age_bwt_co2_grid_z = "z"  # age_bwt_co2_model_bilinear_log.nc variable name (co2)

    x, y = np.meshgrid(
    age_temp_co2_nc[age_bwt_co2_grid_x], age_temp_co2_nc[age_bwt_co2_grid_y])
    kd_tree = KDTree(np.c_[x.ravel(), y.ravel()])


    # ====================== PREPARE BOTTOM-OCEAN TEMPERATURES =========================
    # Read the temperature-time series textfile
    temperature_time_series = pd.read_csv(age_temp_series_path, delimiter="\t")
    ages = temperature_time_series.iloc[:,0].to_numpy()
    temps = temperature_time_series.iloc[:,1].to_numpy()

    # Perform linear interpolation on the time series so that ages are spaced in
    # 1my intervals
    interp_ages = range(int(ages[0]), int(ages[-1]) + 1, 1)
    interp_temps = np.interp(interp_ages, ages, temps)

    # At the current timestep, find the time-averaged bottom water 
    # temperature over an X Myr range (consider X-1 Myr before and 1 Myr in the future).
    # Do this for all linearly interpolated temperatures.
    meanbwt = time_averaged_bottom_water_temperature(interp_temps, average_time_range)

    # Interpolate these time-averaged temperatures and their corresponding
    # ages to the original age grid
    temps_grid = np.interp(ages_grid, interp_ages, meanbwt)

    # Use the age grid's valid age range and the temperature grid to query the tree
    nearest_neighbour_distances, nearest_neighbour_indices = kd_tree.query(np.c_[ages_grid, temps_grid], k=1)


    # ====================== CARBON LOOK-UP ============================================
    # Create a new co2 array and fill it with data queried from the lookup table
    co2_array = np.empty(age_grid_shape)
    co2_array[:] = np.nan

    #print("Commencing carbon lookup at {} Ma...".format(reconstruction_time))

    lookup_indices = np.unravel_index(
        nearest_neighbour_indices, 
        age_temp_co2_nc[age_bwt_co2_grid_z].shape
    )

    np.put(
        co2_array, 
        age_grid_index, 
        np.array(age_temp_co2_nc[age_bwt_co2_grid_z][:][lookup_indices[0], lookup_indices[1]])
    )

    # ===================== WRITE GRIDS TO NETCDF4 =====================================

    gplately.grids.write_netcdf_grid(
        f"{output_path}/{mode}/{reconstruction_time}/upper_crustal_CO2_{mode}_{reconstruction_time}.nc", co2_array
        )

    np.savetxt(
        f"{output_path}/{mode}/{reconstruction_time}/co2_{mode}_mean_{reconstruction_time}.txt",
        [statistics.mean(co2_array[~np.isnan(co2_array)])],
        fmt="%f",
        )
    np.savetxt(
        f"{output_path}/{mode}/{reconstruction_time}/co2_{mode}_std_{reconstruction_time}.txt",
        [statistics.stdev(co2_array[~np.isnan(co2_array)])],
        fmt="%f",
        )
    #print("Created crustal carbon grid for {} Ma!".format(reconstruction_time))
    

