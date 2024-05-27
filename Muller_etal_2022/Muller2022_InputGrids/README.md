After running the following workflows:

- **/utils/A-SeafloorGrids.ipynb**: This generates seafloor age grids and spreading rate grids based on an input plate model.  
- **/utils/Min-Mean-Max-Crustal-Carbon** This generates crustal carbon grids needed to run 01-Sources-of-Carbon.

- Generate **total sediment thickness grids** for 1000-0Ma from [EarthByte's predicting sediment thickness workflow](https://github.com/EarthByte/predicting-sediment-thickness).
- Generate **carbonate sediment thickness grids** for 170-0Ma from [EarthByte's CarbonateSedimentThickness workflow](https://github.com/EarthByte/CarbonateSedimentThickness).
- Generate **contoured continental masks** for 1000-0Ma from [EarthByte's continent contouring workflow](https://github.com/EarthByte/continent-contouring).


The notebooks should be stored in the folders in this directory for ease of running the notebooks (i.e. they have been designed to run with the following directory structure):

1. CarbonateSediment - carbonate sediment thickness grids from [EarthByte's CarbonateSedimentThickness workflow](https://github.com/EarthByte/CarbonateSedimentThickness).
2. ContinentalMasks - passive margins, trenches and contoured continent polygon gpmls, as well as continent masks made from contoured/buffered continental polygons in [EarthByte's continent contouring workflow](https://github.com/EarthByte/continent-contouring). These use the Muller et al. (2022) ref frame. The current version is based on the plate model in webDAV: https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2022_SE/Muller_etal_2022_SE_Merdith21_PMAG-ref-frame-oceanic-crustal-agegrids_v1.2.1.zip

There is also ContinentalMasksNoArcs, which includes masks made from a manually-edited COB-continent combined file in /utils/COB_polygons_and_coastlines_combined_1000_0_Merdith_etal_no_arcs.gpml. We need this file to ignore island arcs when calculating continental platform degassing in notebook 4. 

3. CrustalCarbon - grids made using CO2_review_paper/Muller_etal_2022/utils/Min-Mean-Max-Crustal-Carbon.ipynb using the plate model in webDAV/
4. SeafloorAge - grids made with **/utils/A-SeafloorGrids.ipynb** from the continent masks in /ContinentalMasks/. Considers a full spreading rate in units of mm/yr.
5. SpreadingRate - as for SeafloorAge. Considers a full spreading rate in units of mm/yr.
6. TotalSediment - total sediment thickness grids from [EarthByte's predicting sediment thickness workflow](https://github.com/EarthByte/predicting-sediment-thickness).