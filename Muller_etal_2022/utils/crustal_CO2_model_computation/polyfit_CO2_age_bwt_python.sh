#!/bin/zsh

# Dietmar Muller Sept 2022

# Determine the dependence of oceanic crustal CO2 content on age and botom water temperature
# based on data from Gillis and Coogan, EPSL, 2011.

# --- The resultant model is used to create predicted crustal CO2 grids through time -----

# set mean, min or max CO2
x=mean

# ---- Set input and output file names ----------

dataxyz=age_bwt_co2_weight_$x.txt
modelxym=age_bwt_co2_model_${x}_bilinear_log.xym
modelgrd=age_bwt_co2_model_${x}_bilinear_log.nc

psfile=age_bwt_co2_${x}_model_bilinear_log.ps


# ---- Compute model -----
# ---- fit polynomial N4 bilinear

python grid_trend2d.py -lx -ly -m 4 -xc 200 -yc 200 -xr 0.01 175 -yr 4 20 -w -- $dataxyz $modelxym

frame=0/175/4/20
# xyz2grd $modelxym1 -G$modelgrd1 -R$frame -I1 -V

gmt surface $modelxym -R$frame -Ll0.001 -Lu5 -Z1 -N1000 -G$modelgrd -I0.1 -T1 -V

# ---- Create colour palate ---
cpt=model.cpt
# makecpt -Chaxby -D -T0/5/0.5 -V >$cpt
gmt makecpt -Cno_green -D -T0/3/0.5 -V >$cpt
 
# ---- Plot grdfile of model -------


gmt grdview $modelgrd -Qs -Wc -Z2 -C$cpt -Z2 -Ba20f20:"Age [my]":/a5f1:"Bottom Water Temp [deg C]":WeSn -E-120/30 -Y6 -JZ4 -JX9 -R$frame -P -K >$psfile
# makecpt -Cno_green -D -T0/4/0.5 -V >$cpt
# psscale -L -C$cpt -D8/-1.5/18/.3h -B:"Age and bottom water temp dependence of crustal CO2 content": -O >> $psfile
gmt psscale -L -C$cpt -D13.5/4.5/9/.3 -O >> $psfile

gmt psconvert $psfile -E300 -Tg -A
open $psfile &
