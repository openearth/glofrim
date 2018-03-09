#!/bin/sh
# =========================================================
# Generate input_flood.nam file (simulation settings)
# , modified from original gosh script for GLOFRIM work
# =========================================================

# CaMa-Flood base directory
# CAMADIR="/Users/yamadai/work/CaMa-Flood"
#CAMADIR=`pwd`/..
#cd ${CAMADIR}/gosh

# OPTIONS #
# You can change detailed setting by editing the namelist below.
# for entire list of options, see $CaMa-Flood/mod/mod_input.F

##### Basic Settings ##################
#BASE=$CAMADIR                         #   base directory
#EXP="elbe_15min"

#RDIR=${BASE}/out/$EXP                 #   directory to run CaMa-Flood
#PROG=${BASE}/src/MAIN_day             #   main program
export OMP_NUM_THREADS=4              #   OpenMP cpu num
LFLDOUT=".TRUE."                      #   .TRUE. to activate floodplain discharge
LPTHOUT=".FALSE."                     #   .TRUE. to activate bifurcation flow, mainly for delta simulation
LSTOONLY=".FALSE."                    #   .TRUE. for restart only from storage (no previous-time-step discharge)

##### Time Step #######################   NOTE: spatial resolutions is set by dimension info file "diminfo.txt"
LADPSTP=".TRUE."                      #   .TRUE. for adaptive time step
DT=86400                              #   time step [sec]; set to 86400 for adaptive time step
DTIN=86400                            #   input runoff time step [sec]

##### Simulation Time #################
YSTART=1990                          #   start year (from YSTART / Jan / 1st )
YEND=1990

SPINUP=2                              #   1 for restart, 2 for spinup
NSP=1                                 #   spinup years
#CRESTSTO="${BASE}/"               #   restart file name

##### Map & Topography ################
#FMAP="${BASE}/map/elbe_15min" 
FMAP="."

CDIMINFO="${FMAP}/diminfo_1deg.txt"   #   dimention info (1deg, 0E->360E, 90N-90S)
#CDIMINFO="${FMAP}/diminfo_30min.txt"   #   dimention info (30min, 0E->360E, 90S->90N), generate a new matrix in map dir if needed

CNEXTXY=${FMAP}/nextxy.bin            #   downstream xy (river network map)
CGRAREA=${FMAP}/grarea.bin            #   unit-catchment area [m2]
CELEVTN=${FMAP}/elevtn.bin            #   base elevation      [m]
CNXTDST=${FMAP}/nxtdst.bin            #   downstream distance [m]
CRIVLEN=${FMAP}/rivlen.bin            #   channel length      [m]
CFLDHGT=${FMAP}/fldhgt.bin            #   floodplain elevation profile (height above 'elevtn') [m]

#CRIVWTH=${FMAP}/rivwth.bin            #   channel width       [m] (empirical power-low)
CRIVWTH=${FMAP}/rivwth_gwdlr.bin      #   channel width       [m] (GWD-LR + filled with empirical)
CRIVHGT=${FMAP}/rivhgt.bin            #   channel depth       [m] (empirical power-low)
#CRIVHGT=${FMAP}/rivhgt_pth.bin        #   channel depth       [m] (modified for bifurcation)
#CPTHOUT=${FMAP}/fldpth.txt            #   bifurcation channel list

##### Dike simulation
LDIKHGT=".FALSE."                     #    true: input dike data
#CDIKHGT=${FMAP}/dikhgt.bin            #    dike height [m]

##### Input Runoff Forcing ################# input runoff resolution should be consistent with "inpmat.bin"
LINTERP=".TRUE."                           #   true: runoff interpolation using inpmat, false: nearest point interpolation
#
CINPMAT=${FMAP}/inpmat-1deg.bin            #   runoff input matrix (1deg, 0E->360E, 90N-90S)     !! for sample bonary input
#CINPMAT=${FMAP}/inpmat-30min.bin           #   runoff input matrix (30min, for PCR-GLOBWB)
                                             #   generate a new matrix in map dir if needed

LBMIROF=".TRUE."                           #   true: read runoff via BMI

LINPCDF=".FALSE."                          #   true: netCDF input file
#CROFDIR="${BASE}/inp/ELSE_GPCC/Roff/"      #   runoff directory
#CROFDIR="../../sample_runoff/ELSE_GPCC/Roff/"      #   runoff directory
#CROFPRE="Roff____"                         #   runoff prefix/suffix  ( $(PREFIX)yyyymmdd$(SUFFIX) )
#CROFSUF=".one"

### Now the input runoff is updated in BMI ############
# LINPCDF=".TRUE."
# CROFDIR="${BASE}/inp/ELSE_GPCC/runoff_nc/" #   runoff directory
# CROFCDF="set-by-shell"                     #   netCDF runoff filename
# CROFVAR="runoff"                           #   netCDF runoff variable name
# SYEARIN="set-by-shell"                     #   netCDF runoff file, start date
# SMONIN="set-by-shell"
# SDAYIN="set-by-shell"

DROFUNIT=1.D-3                             #   runoff unit conversion (1.D-3 when input [mm] is converted to [m3/m2])

##### Mean sea level ##################
LMEANSL=".FALSE."                      # true for mean sea level
#CMEANSL="/work_n3/ikeuchi/GTSM/data/bangla_mean_sea_level.bin"
#CMEANSL="/work_n3/ikeuchi/GTSM/data/bangla_mean_sea_level_15min_oldmap.bin"

##### Input Boundary Sea Level ########
LBOUNDSL=".FALSE."                     # true for variable boundary sea level
#CBOUNDDIR="/work_n3/ikeuchi/GTSM/data/GTSM_Sidr_dry/"
#CBOUNDPRE="gtsm"
#CBOUNDSUF=".bin"

##### Output Settings #################
#TODO:THE OUTPUT DIRECTORY NAME SHOULD BE MODIFIED
LOUTCDF=".FALSE."                     # true for netCDF output, false for plain binary output
COUTDIR="./out/"                          # output directory 
#mkdir $COUTDIR
# output variables set "NONE" for no output
CRIVOUTDIR="NONE"                     #   river discharge         [m3/s]
CRIVSTODIR="NONE"                     #   river storage           [m3]
CRIVVELDIR="NONE"                     #   river flow velocity     [m/s]
CRIVDPHDIR="NONE"                 #   river water depth       [m]

CFLDOUTDIR="NONE"                     #   floodplain discharge    [m3/s]
CFLDSTODIR="NONE"                     #   floodplain storage      [m]
CFLDDPHDIR="NONE"                 #   floodplain water depth  [m]
CFLDAREDIR="NONE"                 #   flooded area            [m]
CFLDFRCDIR="NONE"                 #   flooded area fraction   [m2/m2]

CSFCELVDIR="NONE"                 #   water surface elevation [m]
COUTFLWDIR="$COUTDIR"                 #   total discharge (rivout+fldout)   [m3/s]
CSTORGEDIR="NONE"                 #   total storage (rivsto+fldsto)     [m3]

CPTHOUTDIR="NONE"                 #   net bifurcation flow (grid-based) [m3/s]
CPTHFLWDIR="NONE"                 #   bifurcation flow (channel based)  [m3/s]

##### Model Parameters ################
PMANRIV=0.03D0                        # manning coefficient river
PMANFLD=0.10D0                        # manning coefficient floodplain
PCADP=0.7                             # satety coefficient for CFL condition

##### Spatial Resolutions #############
# Set by "diminfo.txt".
# For manual setting, set CDIMINFO="NONE"
# NX="set-by-diminfo"                   #   number of grids in east-west
# NX="set-by-diminfo"                   #   number of grids in east-west
# NLFP="set-by-diminfo"                 #   floodplain layer

# NXIN="set-by-diminfo"                 #   number of input grids in east-west
# NYIN="set-by-diminfo"                 #   number of input grids in east-west
# INPN="set-by-diminfo"                 #   max number of input grids for one cama grid

# WEST="set-by-diminfo"                 #   domain west  edge
# EAST="set-by-diminfo"                 #   domain east  edge
# NORTH="set-by-diminfo"                #   domain north edge
# SOUTH="set-by-diminfo"                #   domain south edge

### End of Setting ####################


## create running dir 
#mkdir -p $RDIR
#cd $RDIR

## if new simulation, remove old files in running directory

#if [ $SPINUP -eq 2 ]; then
#  rm -rf ${RDIR}/????-sp*
#  rm -rf ${RDIR}/*.bin
#  rm -rf ${RDIR}/*.pth
#  rm -rf ${RDIR}/*.vec
#  rm -rf ${RDIR}/*.nc
#  rm -rf ${RDIR}/*.log
#  rm -rf ${RDIR}/*.txt
#  rm -rf ${RDIR}/restart*
#else
#  NSP=0
#fi

## loop 1-year simulation from $YSTART to $YEND

ISP=1         ## spinup count
IYR=$YSTART   ## curent year

#while [ $IYR -le $YEND ];
#do 

if [ $SPINUP -eq 2 ];then
  IRESTART=2                          ## from zero storage
  CRESTSTO=""
else
  IRESTART=1
  CRESTSTO="restart${IYR}0101.bin"    ## from restart file
fi
ISYEAR=$IYR
IEYEAR=`expr $ISYEAR + 1`

#ln -sf $PROG MAIN_day

CROFCDF="${CROFDIR}/runoff${ISYEAR}.nc"
SYEARIN=`expr $ISYEAR - 1`


########### create input namelist ##########
cat > input_flood.nam << EOF
&NRUNVER
IRESTART=$IRESTART                  ! 1=> restart;  2=>spinup
CRESTDIR="./"                       ! restart directory
CRESTSTO="$CRESTSTO"                ! restart file
LSTOONLY=$LSTOONLY                  ! true for restart only from storage
LRESTCDF=.FALSE.                    ! true for netCDF restart file
RESTFREQ=0                          ! 0: yearly restart file, 1: daily restart file
/
&NSIMTIME
ISYEAR=$ISYEAR                      ! start year
ISMON=1                             ! month 
ISDAY=1                             ! day        (assumed at 00UTC)
IEYEAR=$IEYEAR                      ! end year
IEMON=1                             ! end
IEDAY=1                             ! end        (assumed at 00UTC)
/
&NMAP
LMAPCDF=.false.                     ! true for netCDF map input
CDIMINFO="${CDIMINFO}"              ! dimention info
CNEXTXY="${CNEXTXY}"                ! downstream xy (river network map)
CGRAREA="${CGRAREA}"                ! unit-catchment area [m2]
CELEVTN="${CELEVTN}"                ! base elevation      [m]
CNXTDST="${CNXTDST}"                ! downstream distance [m]
CRIVWTH="${CRIVWTH}"                ! channel width       [m]
CRIVLEN="${CRIVLEN}"                ! channel length      [m]
CRIVHGT="${CRIVHGT}"                ! channel depth       [m]
CFLDHGT="${CFLDHGT}"                ! floodplain elevation profile [m]
CPTHOUT="${CPTHOUT}"                ! bifurcation channel list
CRIVCLINC="NONE"                    ! * netCDF river maps
LDIKHGT=${LDIKHGT}                  ! true for dike calculation
CDIKHGT="${CDIKHGT}"                ! dike height         [m]
CRIVPARNC="NONE"                    ! * netCDF river width & depth
/
&NINPUT 
LINTERP=${LINTERP}                  ! true for runoff interpolation using input matrix
LBMIROF=${LBMIROF}                  ! true for reading runoff via BMI
LINPCDF=${LINPCDF}                  ! true for netCDF input
CINPMAT="${CINPMAT}"                ! input matrix file name
CRUNOFFDIR="${CROFDIR}"             ! runoff input directory
CRUNOFFPRE="${CROFPRE}"             ! runoff input prefix
CRUNOFFSUF="${CROFSUF}"             ! runoff input suffix
CRUNOFFCDF="${CROFCDF}"             ! * netCDF input runoff file name
CROFCDFVAR="${CROFVAR}"             ! * netCDF input runoff variable name
SYEARIN=$SYEARIN                    ! * for netCDF input start date (start of the initial time step)
SMONIN=12
SDAYIN=31
LINTERPCDF=.FALSE.                  ! * true for netCDF input matrix
LMEANSL=${LMEANSL}                  ! true for mean sea level
CMEANSL="${CMEANSL}"                ! mean sea level
LBOUNDSL=${LBOUNDSL}                ! true for boundary condition for variable sea level
CBOUNDDIR="${CBOUNDDIR}"            ! boundary sea level directory
CBOUNDPRE="${CBOUNDPRE}"            ! boundary sea level prefix
CBOUNDSUF="${CBOUNDSUF}"            ! boundary sea level suffix
/
&NOUTPUT
LOUTCDF=${LOUTCDF}                  ! true for netCDF output
COUTDIR="${COUTDIR}"                ! output directory ("NONE" for no output)
CRIVOUTDIR="${CRIVOUTDIR}"          ! river discharge        [m3/s]
CRIVSTODIR="${CRIVSTODIR}"          ! river storage          [m3]
CRIVVELDIR="${CRIVVELDIR}"          ! river flow velocity    [m/s]
CRIVDPHDIR="${CRIVDPHDIR}"          ! river water depth      [m]
CFLDOUTDIR="${CFLDOUTDIR}"          ! floodplain discharge   [m3/s]
CFLDSTODIR="${CFLDSTODIR}"          ! floodplain storage     [m3]
CFLDDPHDIR="${CFLDDPHDIR}"          ! floodplain water depth [m]
CFLDFRCDIR="${CFLDFRCDIR}"          ! flooded area fraction  [m2/m2]
CFLDAREDIR="${CFLDAREDIR}"          ! flooded area           [m2]
CSFCELVDIR="${CSFCELVDIR}"          ! water surface elevation           [m]
COUTFLWDIR="${COUTFLWDIR}"          ! total discharge (rivout+fldout)   [m3/s]
CSTORGEDIR="${CSTORGEDIR}"          ! total storage   (rivsto+fldsto)   [m3]
CPTHOUTDIR="${CPTHOUTDIR}"          ! net bifurcation flow (grid-based) [m3/s]
CPTHFLWDIR="${CPTHFLWDIR}"          ! bifurcation flow (channel-based)  [m3/s]
COUTINSDIR="NONE"                   ! instantaneous discharge (no river routing, summation of upstream runoff)
LOUTVEC=.FALSE.                     ! for 1-D land-only output (small data size, post processing required)
/
&NCONF                              ! * NX, NY, NFLP, NXIN, NYIN, INPN, WEST, EAST, NORTH, SOUTH set by diminfo.txt
DT=$DT                              ! time step [sec]
DTIN=$DTIN                          ! input runoff time step [sec]
DROFUNIT=$DROFUNIT                  ! runoff unit conversion (1.D-3 when input [mm] is converted to [m3/m2]
LADPSTP=$LADPSTP                    ! true for adaptive time step
LFLDOUT=$LFLDOUT                    ! true to activate floodplain discharge
LPTHOUT=$LPTHOUT                    ! true to activate bifurcation channel flow
LFLD=.TRUE.                         ! true to activate floodplain inundation
LKINE=.FALSE.                        ! true for kinematic river routing
LMAPEND=.FALSE.                     ! true to convert map data endian
LINPEND=.FALSE.                     ! true to convert input data endian
LLEAPYR=.TRUE.                      ! true for leap year calculatuon, false: always 365days/year
/
&NPARAM
PMANRIV=$PMANRIV                    ! manning coefficient river
PMANFLD=$PMANFLD                    ! manning coefficient floodplain
PGRV=9.8D0                       ! accerelation due to gravity
PDSTMTH=30000.D0                    ! downstream distance at river mouth [m]
PCADP=$PCADP                      ! satety coefficient for CFL condition
PMINSLP=1.D-5                       ! * minimum slope (for kinematic wave)
/
EOF

exit 0

echo "start: ${ISYEAR}" `date` >> log.txt
time ./MAIN_day > run_${ISYEAR}.log 
echo "end:   ${ISYEAR}" `date` >> log.txt

###################

# if curent spinup time $ISP < required spinup time $NSP
#   copy the restart file restart$(IYR+1) to restart$(IYR)
#   copy the outputs to directory "${IYR}-sp1"

SPINUP=1
if [ $IYR -eq $YSTART ];
then
  if [ $ISP -le $NSP ];
    then
    IYR1=`expr $IYR + 1`
    cp restart${IYR1}0101.bin restart${IYR}0101.bin-sp${ISP}
    mv restart${IYR1}0101.bin restart${IYR}0101.bin

    if [ $LPTHOUT == ".TRUE." ];
    then
    cp restart${IYR1}0101.bin.pth restart${IYR}0101.bin.pth-sp${ISP}
    mv restart${IYR1}0101.bin.pth restart${IYR}0101.bin.pth
    fi

    mkdir ${IYR}-sp${ISP}
    mv ??????${IYR}.bin ${IYR}-sp${ISP}
    if [ $LPTHOUT == ".TRUE." ];
    then
    mv ??????${IYR}.pth ${IYR}-sp${ISP}
    fi
    mv ??????${IYR}.nc  ${IYR}-sp${ISP}
    mv *${IYR}.log      ${IYR}-sp${ISP}

    ISP=`expr $ISP + 1`
  else
    ISP=0
    IYR=`expr $IYR + 1`
  fi
else
  IYR=`expr $IYR + 1`
fi
##################

#done # loop to next year simulation

exit 0
