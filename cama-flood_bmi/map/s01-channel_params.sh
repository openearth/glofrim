#!/bin/sh

# Estimating channel parameters (depth, width) by empirical equations of annual discharge
# Please execute this shell script in $CaMa-Flood/map/$MAPNAME directory
# % cd $MAPNAME
# % ../s01-channel_params.sh

################
# input runoff forcings

TYPE='bin'          # plain binary (Fortran direct access)
INTERP='inpmat'     # runoff interpolation using input matrix

DIMINFO='./diminfo_1deg.txt'                 # dimention info to specity CaMa-Flood and Input Data resolutions
CROFBIN='../data/runoff_1981-2000_day.bin'   # Climatology of Daily Runoff (365 day records in mm/day)

#######

## calculate annual discharge
##   (max of 30-day moving average of climatological discharge)

echo ""
echo "@@@ calc_outclm $TYPE $INTERP $DIMINFO $CROFBIN"

../calc_outclm $TYPE $INTERP $DIMINFO $CROFBIN


###  for netCDF runoff climatology  ###
# TYPE='cdf'
# CROFCDF='runoff.nc'  # netCDF runoff climatrogy file
# CROFVAR='RO'         # netCDF runoff variable name

# ../calc_outclm $TYPE $INTERP $DIMINFO $CROFCDF $CROFVAR

##############################################
# Channel bathymetry parameters

HC=0.14                                      ## HC:Coefficient, HP:Power, and HMIN:Minimum for channel depth 
HP=0.40                                      # (H=max(HMIN,HC*Qave**HP)
HMIN=2.00

WC=0.40                                      ## WC:Coefficient, WP:Power, and WMIN:Minimum for channel width
WP=0.75                                      #  (W=max(WMIN,WC*Qave**WP)
WMIN=10.0


## calculate channel parameters 
##   output width: rivwth.bin, depth: rivhgt.bin

echo "" 
echo "@@@ calc_rivwth $TYPE $DIMINFO $HC $HP $HMIN $WC $WP $WMIN"

../calc_rivwth $TYPE $DIMINFO $HC $HP $HMIN $WC $WP $WMIN


##############################################
# Channel width from Global Width Database for Large Rivers (GWD-LR)
#   Width of small channel (<300m) is given from rivwth.bin

## set channel width from GWD-LR
##   output GWD-LR width: rivwth_gwdlr.bin

echo "" 
echo "@@@ set_gwdlr $DIMINFO"

../set_gwdlr $DIMINFO
