#!/bin/sh
AREAS="reg"

IAREA=1
for AREA in $AREAS
do
  IAREA=`expr $IAREA + 1`
  WEST=`awk  'NR==3 {print $'$IAREA'}' ../hires/location.txt`
  NORTH=`awk 'NR==4 {print $'$IAREA'}' ../hires/location.txt`
  NX=`awk    'NR==5 {print $'$IAREA'}' ../hires/location.txt`
  NY=`awk    'NR==6 {print $'$IAREA'}' ../hires/location.txt`
  CSIZE=`awk 'NR==7 {print $'$IAREA'}' ../hires/location.txt`

  SOUTH=`echo "scle=7; $NORTH - $NY * $CSIZE" | bc`

######
  CTL="${AREA}.catmxy.ctl"

  echo "dset ^${AREA}.catmxy"            >  $CTL
  echo 'undef -9999'                     >> $CTL
  echo "title HydroSHEDS $CSIZE deg"     >> $CTL
  echo 'options yrev big_endian'         >> $CTL
  echo "xdef $NX linear $WEST  $CSIZE"   >> $CTL
  echo "ydef $NY linear $SOUTH $CSIZE"   >> $CTL
  echo 'tdef 1 linear 00Z01jan2000 1yr'  >> $CTL
  echo 'zdef 1 linear 1 1'               >> $CTL
  echo 'vars 2'                          >> $CTL
  echo 'catmx 1 -1,40,2,-1 ** catchment x y' >> $CTL
  echo 'catmy 1 -1,40,2,-1 ** 0: no corresponding catchment' >> $CTL
  echo 'ENDVARS'                         >> $CTL

  mv $CTL ../hires/

######
  CTL="${AREA}.flddif.ctl"

  echo "dset ^${AREA}.flddif"            >  $CTL
  echo 'undef -9999'                     >> $CTL
  echo "title HydroSHEDS $CSIZE deg"     >> $CTL
  echo 'options yrev little_endian'      >> $CTL
  echo "xdef $NX linear $WEST  $CSIZE"   >> $CTL
  echo "ydef $NY linear $SOUTH $CSIZE"   >> $CTL
  echo 'tdef 1 linear 00Z01jan2000 1yr'  >> $CTL
  echo 'zdef 1 linear 1 1'               >> $CTL
  echo 'vars 1'                          >> $CTL
  echo 'var 1 99 ** height above channel [m]' >> $CTL
  echo 'ENDVARS'                         >> $CTL

  mv $CTL ../hires/

############
done

