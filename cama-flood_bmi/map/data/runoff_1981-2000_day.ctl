dset  ^runoff_1981-2000_day.bin
undef -9999
title 
options yrev little_endian
xdef  360 linear -179.5 1.
ydef  180 linear  -89.5 1.
tdef  365 linear 00Z01jan2000 1dy
zdef    1 linear 1 1
vars 1
var 1 99       ** daily runoff climatology [mm/day]
ENDVARS




