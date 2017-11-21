#!/bin/sh

# Delineate bifurcation channels by analysing hi-res database

# Please execute this shell script in $CaMa-Flood/map/$MAPNAME directory
# % cd $MAPNAME
# % ../s03-set_bifurcation.sh

#######################
# Specify original hires data dir (HydroSHEDS, GWD-LR, SRTM3)
# Option 1: make a link to FLOW directory
# Option 2: download from dvelopper's web

#SHEDSDIR="../sheds_0.005_140807"
#rm -f sheds
#ln -s $SHEDSDIR sheds


#######################
# set bifurcation channel
#   output bifurcation channel list: fldpth.txt

../set_fldpth

#######################
# modify channel depth
#   output modified depth: rivhgt_pth.bin

../set_rivhgt_pth


