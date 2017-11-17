MODULE INIT_INPUTNAM_MOD
CONTAINS 
SUBROUTINE INIT_INPUTNAM
!==================================================
! set up of namelist inputs.
!==================================================
USE PARKIND1   ,ONLY: JPIM, JPRB
USE MOD_INPUT  ,ONLY: NSETFILE, NULNAM,  TMPNAM, &
                      NRUNVER,  IRESTART,CRESTDIR,CRESTSTO,LRESTCDF,LSTOONLY,RESTFREQ,&
                      NSIMTIME, ISYEAR,  ISMON,   ISDAY,   IEYEAR,  IEMON,   IEDAY,&
                      NMAP,     CDIMINFO,CNEXTXY, CGRAREA, CELEVTN, CNXTDST, CRIVWTH, CRIVLEN, CRIVHGT, CFLDHGT,&
                                CPTHOUT, CRIVCLINC,CRIVPARNC,LMAPCDF, &
                      NINPUT,   LINTERP, LINPCDF, CINPMAT,   CRUNOFFDIR,CRUNOFFPRE,CRUNOFFSUF,CRUNOFFCDF,CROFCDFVAR, &
                                LMEANSL,CMEANSL,LBOUNDSL,CBOUNDDIR,CBOUNDPRE,CBOUNDSUF, &
                      NOUTPUT,  COUTDIR, CRIVOUTDIR,CRIVSTODIR,CRIVVELDIR,CRIVDPHDIR,&
                                         CFLDOUTDIR,CFLDSTODIR,CFLDDPHDIR,CFLDFRCDIR,CFLDAREDIR,&
                                         CSFCELVDIR,COUTFLWDIR,CSTORGEDIR,CPTHOUTDIR,CPTHFLWDIR,&
                                         COUTINSDIR,LOUTCDF,   LOUTVEC,   LOUTYYYY, &
                      NCONF,    NX,NY,NLFP,DT,DT_DEF, NXIN,NYIN,DTIN,DROFUNIT,INPN,WEST,EAST,NORTH,SOUTH,&
                                SYEARIN, SMONIN,  SDAYIN, LFLD,LKINE,LADPSTP,LFLDOUT, LPTHOUT, LDAM, &
                                LMAPEND, LINPEND, LLEAPYR, &
                      NPARAM,   PMANRIV, PMANFLD, PGRV,   PDSTMTH,   PCADP,  PMINSLP, LINTERPCDF

USE MOD_MAP    ,ONLY: REGIONTHIS
USE MOD_OUTPUT ,ONLY: COUTVARDIR

IMPLICIT NONE

! ==============================
! *** 0. SET INPUT UNIT AND OPEN FILE 
OPEN(NSETFILE,FILE="input_flood.nam",STATUS="OLD")

IF(REGIONTHIS==1) then
  WRITE(NULNAM,*) " "
  WRITE(NULNAM,*) "INIT_INPUTNAM: input_flood.nam OPEN in unit:",NSETFILE 
ENDIF
! ==============================
! *** 1. basic simulation setting
! * defaults
IRESTART=1                       !! 1: start from "zero-storage", 2: restart from "restart file"
CRESTDIR="./"                    !! restart file directory
CRESTSTO="restartin"             !! restart file name
LRESTCDF=.FALSE.                 !! true for netCDF restart file
LSTOONLY=.FALSE.                 !! true for restart only from storage (no discharge restart, not strict)
RESTFREQ=0                       !! 0: restart file at end of simulation, 1: daily restart file

! * change
READ(NSETFILE,NRUNVER)
IF(REGIONTHIS==1) then
  WRITE(NULNAM,*) " "
  WRITE(NULNAM,*) "****INIT_INPUTNAM: 1. INIT_RUNVER****"
  WRITE(NULNAM,*) "IRESTART=",IRESTART
  WRITE(NULNAM,*) "CRESTDIR=",TRIM(CRESTDIR)
  WRITE(NULNAM,*) "CRESTSTO=",TRIM(CRESTSTO)
  WRITE(NULNAM,*) "LRESTCDF=",LRESTCDF
  WRITE(NULNAM,*) "LSTOONLY=",LSTOONLY
  WRITE(NULNAM,*) "RESTFREQ=",RESTFREQ
ENDIF


! ==============================
! *** 2. set simulation time 
! * defaults
ISYEAR=1990
ISMON=1
ISDAY=1
IEYEAR=1991
IEMON=1
IEDAY=1

! * change
REWIND(NSETFILE)
READ(NSETFILE,NSIMTIME)
IF(REGIONTHIS==1) then
  WRITE(NULNAM,*) " "
  WRITE(NULNAM,*) "****INIT_INPUTNAM: 2. SET_SIMTIME****"
  WRITE(NULNAM,*) "ISYEAR,ISMON,ISDAY:",ISYEAR,ISMON,ISDAY
  WRITE(NULNAM,*) "IEYEAR,IEMON,IEDAY:",IEYEAR,IEMON,IEDAY
ENDIF

! ==============================
! *** 3. set map
! * defaults
CDIMINFO="NONE"
CNEXTXY="./nextxy.bin "
CGRAREA="./grarea.bin"
CELEVTN="./elevtn.bin"
CNXTDST="./nxtdst.bin "
CRIVWTH="./rivwth.bin"
CRIVLEN="./rivlen.bin"
CRIVHGT="./rivhgt.bin"
CFLDHGT="./fldhgt.bin"
CPTHOUT="./fldpath.txt"
CRIVCLINC="NONE"
CRIVPARNC="NONE"
LMAPCDF=.FALSE.

! * change
REWIND(NSETFILE)
READ(NSETFILE,NMAP)
IF(REGIONTHIS==1) then
  WRITE(NULNAM,*) " "
  WRITE(NULNAM,*) "****INIT_INPUTNAM: 3. SET_SETMAP****"
  WRITE(NULNAM,*) "CDIMINFO:",TRIM(CDIMINFO)
  WRITE(NULNAM,*) "CNEXTXY:",TRIM(CNEXTXY)
  WRITE(NULNAM,*) "CGRAREA:",TRIM(CGRAREA)
  WRITE(NULNAM,*) "CELEVTN:",TRIM(CELEVTN)
  WRITE(NULNAM,*) "CNXTDST:",TRIM(CNXTDST)
  WRITE(NULNAM,*) "CRIVWTH:",TRIM(CRIVWTH)
  WRITE(NULNAM,*) "CRIVLEN:",TRIM(CRIVLEN)
  WRITE(NULNAM,*) "CRIVHGT:",TRIM(CRIVHGT)
  WRITE(NULNAM,*) "CFLDHGT:",TRIM(CFLDHGT)
  WRITE(NULNAM,*) "CPTHOUT:",TRIM(CPTHOUT)

  WRITE(NULNAM,*) "LMAPCDF:"  ,LMAPCDF
  WRITE(NULNAM,*) "CRIVCLINC:",TRIM(CRIVCLINC)
  WRITE(NULNAM,*) "CRIVPARNC:",TRIM(CRIVPARNC)
ENDIF

! ==============================
! *** 4. set INPUT
! * defaults
LINTERP=.FALSE.                 !! use nearest interpolation (.true. for interpolation with input matrix)
CINPMAT="NONE"

LINPCDF=.FALSE.                 !! true for netCDF runoff input

CRUNOFFDIR="./runoff/"
CRUNOFFPRE="Roff____"           !! defaults runoff file name Roff____YYYYMMDD.one
CRUNOFFSUF=".one"

CRUNOFFCDF="NONE"
CROFCDFVAR="runoff"
SYEARIN=0                       !! netCDF input file start date (set to 0 when not used)
SMONIN=0
SDAYIN=0
LINTERPCDF=.FALSE.              !! true for netCDF input matrix

LMEANSL=.FALSE.
CMEANSL="./meansl/"

LBOUNDSL=.FALSE.
CBOUNDDIR="./bound/"
CBOUNDPRE="bound"
CBOUNDSUF=".bin"

! * change
REWIND(NSETFILE)
READ(NSETFILE,NINPUT)
IF(REGIONTHIS==1) then
  WRITE(NULNAM,*) " "
  WRITE(NULNAM,*) "****INIT_INPUTNAM: 4. SET_INPUT****"
  WRITE(NULNAM,*) "LINTERP:"   ,LINTERP
  WRITE(NULNAM,*) "LINTERPCDF:"   ,LINTERPCDF
  WRITE(NULNAM,*) "CINPMAT:",   TRIM(CINPMAT)
  WRITE(NULNAM,*) "LINPCDF:"   ,LINPCDF
  WRITE(NULNAM,*) "CRUNOFFDIR:",TRIM(CRUNOFFDIR)
  WRITE(NULNAM,*) "CRUNOFFPRE:",TRIM(CRUNOFFPRE)
  WRITE(NULNAM,*) "CRUNOFFSUF:",TRIM(CRUNOFFSUF)
  WRITE(NULNAM,*) "CRUNOFFCDF:",TRIM(CRUNOFFCDF)
  WRITE(NULNAM,*) "CROFCDFVAR:",TRIM(CROFCDFVAR)
  WRITE(NULNAM,*) "SYEARIN,SMONIN,SDAYIN",SYEARIN,SMONIN,SDAYIN
  WRITE(NULNAM,*) "LMEANSL:",LMEANSL
  WRITE(NULNAM,*) "CMEANSL:",TRIM(CMEANSL)
  WRITE(NULNAM,*) "LBOUNDSL:",LBOUNDSL
  WRITE(NULNAM,*) "CBOUNDDIR:",TRIM(CBOUNDDIR)
  WRITE(NULNAM,*) "CBOUNDPRE:",TRIM(CBOUNDPRE)
  WRITE(NULNAM,*) "CBOUNDSUF:",TRIM(CBOUNDSUF)
ENDIF

! ==============================
! *** 5. set output
! * defaults
COUTDIR="./"
CRIVOUTDIR="NONE"
CRIVSTODIR="NONE"
CRIVVELDIR="NONE"
CRIVDPHDIR="./"

CFLDOUTDIR="NONE"
CFLDSTODIR="NONE"
CFLDDPHDIR="NONE"
CFLDFRCDIR="./"
CFLDAREDIR="./"

CSFCELVDIR="./"
COUTFLWDIR="./"
CSTORGEDIR="./"
COUTINSDIR="NONE"

CPTHOUTDIR="NONE"
CPTHFLWDIR="NONE"

LOUTVEC=.FALSE.
LOUTCDF=.FALSE.
LOUTYYYY=.TRUE.

! * change
REWIND(NSETFILE)
READ(NSETFILE,NOUTPUT)

IF(REGIONTHIS==1) then
  WRITE(NULNAM,*) " "
  WRITE(NULNAM,*) "****INIT_INPUTNAM: 5. SET_OUTPUT****"
  WRITE(NULNAM,*) "COUTDIR:",TRIM(COUTDIR)
  WRITE(NULNAM,*) "CRIVOUTDIR:",TRIM(CRIVOUTDIR)
  WRITE(NULNAM,*) "CRIVSTODIR:",TRIM(CRIVSTODIR)
  WRITE(NULNAM,*) "CRIVVELDIR:",TRIM(CRIVVELDIR)
  WRITE(NULNAM,*) "CRIVDPHDIR:",TRIM(CRIVDPHDIR)

  WRITE(NULNAM,*) "CFLDOUTDIR:",TRIM(CFLDOUTDIR)
  WRITE(NULNAM,*) "CFLDSTODIR:",TRIM(CFLDSTODIR)
  WRITE(NULNAM,*) "CFLDDPHDIR:",TRIM(CFLDDPHDIR)
  WRITE(NULNAM,*) "CFLDFRCDIR:",TRIM(CFLDFRCDIR)
  WRITE(NULNAM,*) "CFLDAREDIR:",TRIM(CFLDAREDIR)

  WRITE(NULNAM,*) "CSFCELVDIR:",TRIM(CSFCELVDIR)
  WRITE(NULNAM,*) "COUTFLWDIR:",TRIM(COUTFLWDIR)
  WRITE(NULNAM,*) "CSTORGEDIR:",TRIM(CSTORGEDIR)
  WRITE(NULNAM,*) "COUTINSDIR:",TRIM(COUTINSDIR)

  WRITE(NULNAM,*) "CPTHOUTDIR:",TRIM(CPTHOUTDIR)
  WRITE(NULNAM,*) "CPTHFLWDIR:",TRIM(CPTHFLWDIR)

  WRITE(NULNAM,*) "LOUTVEC:",LOUTVEC
  WRITE(NULNAM,*) "LOUTCDF:",LOUTCDF
  WRITE(NULNAM,*) "LOUTYYYY:",LOUTYYYY
ENDIF

IF( LOUTCDF ) LOUTVEC=.FALSE.

! ==============================
! *** 6. set CONF: GRID TSTEP CONFIGURATION
! * defaults
NX=1440                                     !! 15 minute resolution
NY=720
NLFP=10                                     !! 10 floodplain layer
DT=24*60*60                                 !! dt = 1day (automatically set by adaptive time step)

NXIN=360                                    !! 1 degree input
NYIN=180
INPN=1                                      !! maximum number of input grids corresponding to one CaMa-Flood grid
DTIN=24*60*60                               !! default: /day -> /sec
DROFUNIT=1.D-3                              !! default: mm   -> m3/m2

WEST=-180.D0                                !! west, east, north, south edges of the domain
EAST= 180.D0
NORTH= 90.D0
SOUTH=-90.D0

LFLD=.TRUE.     ! floodplain scheme active     (.FALSE. for no floodplain )
LKINE=.FALSE.   ! use local inertial equation  (.TRUE.  for kinematic wave)
LADPSTP=.TRUE.  ! use adaptive time step       (.FALSE. for no aptive time step)
LFLDOUT=.TRUE.  ! floodplan flow active        (.FALSE. for no floodplain flow)
LPTHOUT=.FALSE. ! floodplan path flow          (.FALSE. for no floodplain path flow)

LDAM=.FALSE.    ! dam operation                (.TRUE.  for dam operation: NOT SUPPORTED)
LMAPEND=.FALSE. ! map data     in different endian (.FALSE. for no endian conversion)
LINPEND=.FALSE. ! input runoff in different endian (.FALSE. for no endian conversion)

LLEAPYR=.TRUE.  ! false to always use 365days/year

! * change
REWIND(NSETFILE)
READ(NSETFILE,NCONF)

IF( CDIMINFO/="NONE" )THEN
  WRITE(NULNAM,*) " Reading Dimention Info: ", TRIM(CDIMINFO)
  OPEN(TMPNAM,FILE=CDIMINFO,FORM='FORMATTED')
  READ(TMPNAM,*) NX
  READ(TMPNAM,*) NY
  READ(TMPNAM,*) NLFP
  READ(TMPNAM,*) NXIN
  READ(TMPNAM,*) NYIN
  READ(TMPNAM,*) INPN
  READ(TMPNAM,*) 
  READ(TMPNAM,*) WEST
  READ(TMPNAM,*) EAST
  READ(TMPNAM,*) NORTH
  READ(TMPNAM,*) SOUTH
  CLOSE(TMPNAM)
ENDIF

IF(REGIONTHIS==1)THEN
  WRITE(NULNAM,*) " "
  WRITE(NULNAM,*) "****INIT_INPUTNAM: 6. SET_CONF****"
  WRITE(NULNAM,*) "NX, NY,NLFP",NX,NY,NLFP
  WRITE(NULNAM,*) "DT",DT
  WRITE(NULNAM,*) "WEST,EAST,NORTH,SOUTH",WEST,EAST,NORTH,SOUTH
  WRITE(NULNAM,*) "NXIN,NYIN,DTIN,DROFUNIT",NXIN,NYIN,DTIN,DROFUNIT
  WRITE(NULNAM,*) "LFLD",LFLD
  WRITE(NULNAM,*) "LKINE",LKINE
  WRITE(NULNAM,*) "LADPSTP",LADPSTP
  WRITE(NULNAM,*) "LFLDOUT",LFLDOUT
  WRITE(NULNAM,*) "LPTHOUT",LPTHOUT
  WRITE(NULNAM,*) "LLEAPYR",LLEAPYR
ENDIF

DT_DEF=DT                                   !! the default DT for adaptope time step mode

! *CHECK
IF ( .not.LFLD .AND. .not.LKINE ) THEN
  WRITE(NULNAM,*) "LFLD=.false. & LKINE=.false."
  WRITE(NULNAM,*) "no floodplain run only available with kinematic wave (LKINE=.true.)"
  WRITE(NULNAM,*) "stop"
  STOP 9
ENDIF

IF ( LKINE .AND. LADPSTP ) THEN
  WRITE(NULNAM,*) "LKINE=.true. & LADPSTP=.true."
  WRITE(NULNAM,*) "adaptive time step only available with local inertial equation (LKINE=.false.)"
  WRITE(NULNAM,*) "STOP"
  STOP 9
ENDIF

IF ( LKINE .AND. LPTHOUT ) THEN
  WRITE(NULNAM,*) "LKINE=.true. & LPATHOUT=.true."
  WRITE(NULNAM,*) "bifurcation channel flow only available with local inertial equation (LKINE=.false.)"
  WRITE(NULNAM,*) "STOP"
  STOP 9
ENDIF

IF( LRESTCDF .and. LPTHOUT )THEN
  LSTOONLY=.TRUE.
  WRITE(NULNAM,*) "LRESTCDF=.true. & LPATHOUT=.true."
  WRITE(NULNAM,*) "LSTOONLY set to .TRUE."
ENDIF

IF( .not. LPTHOUT )THEN
  CPTHOUTDIR='NONE'
  CPTHFLWDIR='NONE'
ENDIF

IF( LOUTCDF )THEN
  CPTHFLWDIR='NONE'
ENDIF

COUTVARDIR(1)=TRIM(CRIVOUTDIR)
COUTVARDIR(2)=TRIM(CRIVSTODIR)
COUTVARDIR(3)=TRIM(CRIVDPHDIR)
COUTVARDIR(4)=TRIM(CRIVVELDIR)

COUTVARDIR(5)=TRIM(CFLDOUTDIR)
COUTVARDIR(6)=TRIM(CFLDSTODIR)
COUTVARDIR(7)=TRIM(CFLDDPHDIR)
COUTVARDIR(8)=TRIM(CFLDFRCDIR)
COUTVARDIR(9)=TRIM(CFLDAREDIR)

COUTVARDIR(10)=TRIM(CSFCELVDIR)
COUTVARDIR(11)=TRIM(COUTFLWDIR)
COUTVARDIR(12)=TRIM(CSTORGEDIR)
COUTVARDIR(13)=TRIM(COUTINSDIR)

COUTVARDIR(14)=TRIM(CPTHOUTDIR)
COUTVARDIR(15)=TRIM(CPTHFLWDIR)


! ==============================
! *** 7. set PARAM: parameters
! * defaults
PMANRIV=0.03D0                              !! manning coefficient river
PMANFLD=0.10D0                              !! manning coefficient floodplain
PGRV   =9.8D0                               !! gravity accerelation
PDSTMTH=10000.D0                            !! downstream distance at river mouth [m]
PCADP  =0.7D0                                !! CFL coefficient
PMINSLP=1.D-5                               !! minimum slope (kinematic wave)

! * change
REWIND(NSETFILE)
READ(NSETFILE,NPARAM)

IF(REGIONTHIS==1)THEN
  WRITE(NULNAM,*) " "
  WRITE(NULNAM,*) "****INIT_INPUTNAM: 7. PARAMETERS****"
  WRITE(NULNAM,*) "PMANRIV, PMANFLD", PMANRIV, PMANFLD
  WRITE(NULNAM,*) "PGRV   ",          PGRV
  WRITE(NULNAM,*) "PDSTMTH",          PDSTMTH
  WRITE(NULNAM,*) "PCADP  ",          PCADP
  WRITE(NULNAM,*) "PMINSLP",          PMINSLP
ENDIF


! ==============================
! *** CLOSE FILE 
CLOSE(NSETFILE)




END SUBROUTINE INIT_INPUTNAM

END MODULE INIT_INPUTNAM_MOD
