<!-- #region -->
# Carbon review paper
Workflows from Muller et al. (in review) "Solid Earth carbon degassing and sequestration since 1 billion years ago".

## Dependencies

- [pyGPlates](https://www.gplates.org/docs/pygplates/pygplates_getting_started.html#installation)
- [PlateTectonicTools](https://github.com/EarthByte/PlateTectonicTools)
- [gplately](https://github.com/GPlates/gplately/tree/master)
- [Slab-Dip](https://github.com/brmather/Slab-Dip)
- [melt](https://github.com/brmather/melt)
- joblib *(these workflows use parallelisation with joblib and cannot be run on a single thread).*

## Input files
Run the following notebooks and workflows to prepare the input files required to run these notebooks. 

- **/utils/A-SeafloorGrids.ipynb**: This generates seafloor age grids and spreading rate grids based on an input plate model.  
- **/utils/Min-Mean-Max-Crustal-Carbon** This generates crustal carbon grids needed to run 01-Sources-of-Carbon.

- Generate **total sediment thickness grids** for 1000-0Ma from [EarthByte's predicting sediment thickness workflow](https://github.com/EarthByte/predicting-sediment-thickness).
- Generate **carbonate sediment thickness grids** for 170-0Ma from [EarthByte's CarbonateSedimentThickness workflow](https://github.com/EarthByte/CarbonateSedimentThickness).
- Generate **contoured continental masks** for 1000-0Ma from [EarthByte's continent contouring workflow](https://github.com/EarthByte/continent-contouring).

The notebooks will reference these grids using the following directory structure:

```
    # Change this: Directories to age grids and spreading rate grids
    grid_directory = "./Muller2022_InputGrids/"
    spreadrate_filename = grid_directory+"SpreadingRate/Muller2022_SPREADING_RATE_grid_{:.1f}Ma.nc"
    agegrid_filename = grid_directory+"SeafloorAge/Muller2022_SEAFLOOR_AGE_grid_{:.1f}Ma.nc"
```
This can be changed to suit the directory made for the input grids above, or the input grids can be saved to the following directory structure for ease of running the notebooks (i.e. they have been designed to run with the following directory structure):

1. CarbonateSediment - carbonate sediment thickness grids from [EarthByte's CarbonateSedimentThickness workflow](https://github.com/EarthByte/CarbonateSedimentThickness).
2. ContinentalMasks - passive margins, trenches and contoured continent polygon gpmls, as well as continent masks made from contoured/buffered continental polygons in [EarthByte's continent contouring workflow](https://github.com/EarthByte/continent-contouring). These use the Muller et al. (2022) ref frame. The current version is based on the plate model in webDAV: https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2022_SE/Muller_etal_2022_SE_Merdith21_PMAG-ref-frame-oceanic-crustal-agegrids_v1.2.1.zip

There is also ContinentalMasksNoArcs, which includes masks made from a manually-edited COB-continent combined file in /utils/COB_polygons_and_coastlines_combined_1000_0_Merdith_etal_no_arcs.gpml. We need this file to ignore island arcs when calculating continental platform degassing in notebook 4. 

3. CrustalCarbon - grids made using CO2_review_paper/Muller_etal_2022/utils/Min-Mean-Max-Crustal-Carbon.ipynb using the plate model in webDAV/
4. SeafloorAge - grids made with **/utils/A-SeafloorGrids.ipynb** from the continent masks in /ContinentalMasks/. Considers a full spreading rate in units of mm/yr.
5. SpreadingRate - as for SeafloorAge. Considers a full spreading rate in units of mm/yr.
6. TotalSediment - total sediment thickness grids from [EarthByte's predicting sediment thickness workflow](https://github.com/EarthByte/predicting-sediment-thickness).


## Notebooks: 
Once the input files above have been generated, run the following notebooks in order:


1. 01-Sources-of-Carbon
2. 02-Subducted-Carbon
3. 03-Carbon-Degassing

4. /H2O_review_paper/01-Sources-of-Water
5. /H2O_review_paper/02-Subducted-Water
6. /H2O_review_paper/06-Maps-of-water-storage.ipynb

7. 04-Carbonate-Platform-Degassing
8. 05-Atmospheric-Carbon
9. 06-PlateTectonicTools

Some auxiliary notebooks for animations and plots are in /utils/:
- 0A-Cumulative-Subducted-Carbon.ipynb
- 0B-Miscellaneous-Plots.ipynb
- 0C-Carbon-Panel-Plots.ipynb
- 0D-Carbon-Panel-Plot-Videos.ipynb


<!-- #endregion -->

```python

```
