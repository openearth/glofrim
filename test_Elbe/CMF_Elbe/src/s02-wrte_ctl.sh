#!/bin/sh

./gradsinfo > tmp
NX=`awk    '{print $1}' tmp`
NY=`awk    '{print $2}' tmp`
GSIZE=`awk '{print $3}' tmp`
WEST=`awk  '{print $4}' tmp`
SOUTH=`awk '{print $5}' tmp`

rm -f tmp

VARS='elevtn fldhgt grarea nxtdst rivout rivlen rivhgt rivwth width basin bsncol rivseq uparea upgrid lsmask'
for VAR in $VARS
do
  TYPE=`awk  -v AVAR="$VAR" '$1==AVAR {print $2}' ./varname.txt`
  LAYER=`awk -v AVAR="$VAR" '$1==AVAR {print $3}' ./varname.txt`
  EXPL=`awk  -v AVAR="$VAR" '$1==AVAR {print $4}' ./varname.txt`
  UNIT=`awk  -v AVAR="$VAR" '$1==AVAR {print $5}' ./varname.txt`

  if [ "${TYPE}" = "Int" ];then
    UNDEF="-9999"
    TFORM="-1,40,4"
  else
    UNDEF="-9999"
    TFORM="99"
  fi

  CTLFILE="../${VAR}.ctl"
  rm -f $CTLFILE

  echo "dset    ^${VAR}.bin"                 >> $CTLFILE
  echo "undef   ${UNDEF}"                    >> $CTLFILE
  echo "title   Flow_River_Network"          >> $CTLFILE
  echo "options yrev little_endian"          >> $CTLFILE
  echo "xdef ${NX} linear ${WEST}   ${GSIZE}"  >> $CTLFILE
  echo "ydef ${NY} linear ${SOUTH}  ${GSIZE}" >> $CTLFILE
  echo "tdef 1 linear 00Z01jan2000 1dy"      >> $CTLFILE
  echo "zdef ${LAYER} linear 1 1"            >> $CTLFILE
  echo "vars 1"                              >> $CTLFILE
  echo "var ${LAYER} ${TFORM} ** ${EXPL} ${UNIT}" >> $CTLFILE
  echo "ENDVARS"                             >> $CTLFILE
done


##########################################################################

VARS='nextxy lonlat'
for VAR in $VARS
do
  TYPE=`awk  -v AVAR="$VAR" '$1==AVAR {print $2}' ./varname.txt`
  LAYER=`awk -v AVAR="$VAR" '$1==AVAR {print $3}' ./varname.txt`
  EXPL=`awk  -v AVAR="$VAR" '$1==AVAR {print $4}' ./varname.txt`
  UNIT=`awk  -v AVAR="$VAR" '$1==AVAR {print $5}' ./varname.txt`

  if [ "${TYPE}" = "Int" ];then
    UNDEF="-9999"
    TFORM="-1,40,4"
  else
    UNDEF="-9999"
    TFORM="99"
  fi

  if [ "${VAR}" = "nextxy" ];then
    VAR1="nextx"
    VAR2="nexty"
  else
    VAR1="dlon"
    VAR2="dlat"
  fi

  CTLFILE="../${VAR}.ctl"
  rm -f $CTLFILE

  echo "dset    ^${VAR}.bin"                 >> $CTLFILE
  echo "undef   ${UNDEF}"                    >> $CTLFILE
  echo "title   Flow_River_Network"          >> $CTLFILE
  echo "options yrev little_endian"          >> $CTLFILE
  echo "xdef ${NX} linear ${WEST}   ${GSIZE}"  >> $CTLFILE
  echo "ydef ${NY} linear ${SOUTH}  ${GSIZE}" >> $CTLFILE
  echo "tdef 1 linear 00Z01jan2000 1dy"      >> $CTLFILE
  echo "zdef 1 linear 1 1"                   >> $CTLFILE
  echo "vars 2"                              >> $CTLFILE
  echo "${VAR1} 1 ${TFORM} ** ${EXPL} ${UNIT}"   >> $CTLFILE
  echo "${VAR2} 1 ${TFORM} ** "                  >> $CTLFILE
  echo "ENDVARS"                             >> $CTLFILE
done


##########################################################################

VAR='inpmat-1deg'
INPN=`awk 'NR==6 {print $1}' ../diminfo_1deg.txt`

  CTLFILE="../${VAR}.ctl"
  rm -f $CTLFILE

  echo "dset    ^inpmat-1deg.bin"            >> $CTLFILE
  echo "undef   0"                           >> $CTLFILE
  echo "title   Flow_River_Network"          >> $CTLFILE
  echo "options yrev little_endian"          >> $CTLFILE
  echo "xdef ${NX} linear ${WEST}   ${GSIZE}" >> $CTLFILE
  echo "ydef ${NY} linear ${SOUTH}  ${GSIZE}" >> $CTLFILE
  echo "tdef 1 linear 00Z01jan2000 1dy"      >> $CTLFILE
  echo "zdef ${INPN} linear 1 1"             >> $CTLFILE
  echo "vars 3"                              >> $CTLFILE
  echo "inpx ${INPN} -1,40,4 ** input grid ixin"   >> $CTLFILE
  echo "inpy ${INPN} -1,40,4 ** input grid iyin"   >> $CTLFILE
  echo "inpa ${INPN} 99      ** input grid area"   >> $CTLFILE
  echo "ENDVARS"                             >> $CTLFILE
