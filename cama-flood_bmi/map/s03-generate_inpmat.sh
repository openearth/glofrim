#!/bin/sh

# Generate input matrix

# Please execute this shell script in $CaMa-Flood/map/$MAPNAME directory
# % cd $MAPNAME
# % ../s03-generate_inpmat.sh

#######################
# Specify input gridsize, input domain edge, input north-south order, diminfo & inpmat name

# Default setting
GRSIZEIN=0.5       # input grid size
WESTIN=-180.0      # input domain west east north south edge
EASTIN=180.0
NORTHIN=90.0
SOUTHIN=-90.0
OLAT="StoN"        # north-south order of input data
# OLAT="NtoS"

# For the test run of Elbe river with input from regional PCR input
#WESTIN=7.0      # input domain west east north south edge
#EASTIN=17.0
#NORTHIN=55.0
#SOUTHIN=48.0

DIMINFO="diminfo_30min.txt"
INPMATBIN="inpmat-30min.bin"
INPMATTXT="inpmat-30min.txt"

#######################
# generate inpmat
#   Fortran code is for linear Cartesian input
#   for other grid coordinate, code should be reqritten

../generate_inpmat $GRSIZEIN $WESTIN $EASTIN $NORTHIN $SOUTHIN $OLAT

mv diminfo_tmp.txt $DIMINFO
mv inpmat-tmp.bin  $INPMATBIN
mv inpmat-tmp.txt  $INPMATTXT

