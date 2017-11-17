MODULE MOD_INPUT
USE PARKIND1 ,ONLY: JPIM, JPRB, JPRM
IMPLICIT NONE

!!==================================================
! **** INPUT VARIABLES FROM NAMELISTS

SAVE 

INTEGER(KIND=JPIM)     ::  NSETFILE        !! input namelist
INTEGER(KIND=JPIM)     ::  NULNAM          !! DEFAULT OUTPUT UNIT (6)
INTEGER(KIND=JPIM)     ::  TMPNAM          !! DEFAULT LOG    UNIT (11)
INTEGER(KIND=JPIM)     ::  LOGNAM          !! DEFAULT LOG    UNIT (12)
CHARACTER(LEN=256)     ::  CLOGOUT         !! LOGFILE
PARAMETER(NSETFILE=11)
PARAMETER(NULNAM=6)
PARAMETER(TMPNAM=7)
DATA LOGNAM /6/                            !! modified to 12 when MPI

!!=== NAMELIST/NRUNVER/
INTEGER(KIND=JPIM)     ::  IRESTART        !! 1=> RESTART;  2=>SPINUP 
CHARACTER(LEN=256)     ::  CRESTDIR        !! RESTART DIRECTORY
CHARACTER(LEN=256)     ::  CRESTSTO        !! RESTART STORAGE FILE 
LOGICAL                ::  LRESTCDF        !! TRUE FOR NETCDF RESTART IN AND OUT
LOGICAL                ::  LSTOONLY        !! TRUE FOR LOOSE RESTART (ONLY FROM STORAGE)
INTEGER(KIND=JPIM)     ::  RESTFREQ        !! RESTART WRITE FRQUENCY: 
                                           !! "0": AT THE END ; 1-DAYLY :(ONLY FOR NETCDF ) BINARY 1 PER YEAR ... 

!!=== NAMELIST/NSIMTIME/
INTEGER(KIND=JPIM)     ::  ISYEAR          !! START YEAR
INTEGER(KIND=JPIM)     ::  ISMON           !! START MONTH
INTEGER(KIND=JPIM)     ::  ISDAY           !! START DAY
INTEGER(KIND=JPIM)     ::  IEYEAR          !! END YEAR
INTEGER(KIND=JPIM)     ::  IEMON           !! END MONTH
INTEGER(KIND=JPIM)     ::  IEDAY           !! END DAY

!!=== NAMELIST/NMAP/
CHARACTER(LEN=256)     ::  CDIMINFO        !! Dimention Information
CHARACTER(LEN=256)     ::  CNEXTXY         !! NEXTX FILE
CHARACTER(LEN=256)     ::  CGRAREA         !! GRAREA FILE
CHARACTER(LEN=256)     ::  CELEVTN         !! ELEVATION FILE
CHARACTER(LEN=256)     ::  CNXTDST         !! NEXT DISTANCE FILE
CHARACTER(LEN=256)     ::  CRIVWTH         !! RIVER WIDTH
CHARACTER(LEN=256)     ::  CRIVLEN         !! RIVER LENGTH
CHARACTER(LEN=256)     ::  CRIVHGT         !! RIVER HEIGHT
CHARACTER(LEN=256)     ::  CFLDHGT         !! FLOODPLAIN HEIGHT
CHARACTER(LEN=256)     ::  CPTHOUT         !! FLOODPLAIN PATH
CHARACTER(LEN=256)     ::  CRIVCLINC       !! RIV CLIM NETCDF FILE ('none' FOR BINARY)
CHARACTER(LEN=256)     ::  CRIVPARNC       !! RIV PARAMETERS (WIDTH,HEIGHT) NETCDF FILE ('none' FOR BINARY)
LOGICAL                ::  LMAPCDF         !! true for netCDF map input
LOGICAL                ::  LDIKHGT         !! true for dike calculation
CHARACTER(LEN=256)     ::  CDIKHGT         !! RIVER DIKE HEIGHT

!!=== NAMELIST/NINPUT/
LOGICAL                ::  LINTERP             !! true : runoff interpolation using input matrix
CHARACTER(LEN=256)     ::  CINPMAT             !! runoff input matrix

LOGICAL                ::  LINPCDF             !! true : netCDF runoff input

CHARACTER(LEN=256)     ::  CRUNOFFDIR          !! INPUT  SURFACE RUNOFF DIRECTORY
CHARACTER(LEN=256)     ::  CRUNOFFPRE          !! INPUT  SURFACE RUNOFF PREFIX
CHARACTER(LEN=256)     ::  CRUNOFFSUF          !! INPUT  SURFACE RUNOFF SUFFIX

LOGICAL                ::  LMEANSL             !! true : mean sea level
CHARACTER(LEN=256)     ::  CMEANSL             !! INPUT MEAN SEA LEVEL

LOGICAL                ::  LBOUNDSL            !! true : boundary condition for variable sea level
CHARACTER(LEN=256)     ::  CBOUNDDIR           !! INPUT BOUNDARY SEA LEVEL DIRECTORY
CHARACTER(LEN=256)     ::  CBOUNDPRE           !! INPUT BOUNDARY SEA LEVEL PREFIX
CHARACTER(LEN=256)     ::  CBOUNDSUF           !! INPUT BOUNDARY SEA LEVEL SUFFIX

LOGICAL                ::  LINTERPCDF          !! true : runoff input matrix in netcdf 
CHARACTER(LEN=256)     ::  CRUNOFFCDF          !! INPUT  netcdf file with runoff "NONE" FOR BINARY
CHARACTER(LEN=256)     ::  CROFCDFVAR          !! VARNAME of runoff in netCDF input file "default runoff"
INTEGER(KIND=JPIM)     ::  SYEARIN             !! START YEAR IN NETCDF INPUT RUNOFF
INTEGER(KIND=JPIM)     ::  SMONIN              !! START MONTH IN NETCDF INPUT RUNOFF
INTEGER(KIND=JPIM)     ::  SDAYIN              !! START DAY IN NETCDF INPUT RUNOFF

!!=== NAMELIST/NOUTPUT/
CHARACTER(LEN=256)     ::  COUTDIR             !! OUTPUT DIRECTORY
CHARACTER(LEN=256)     ::  CRIVOUTDIR          !! OUTPUT RIVER OUTFLOW  DIRECTORY
CHARACTER(LEN=256)     ::  CRIVSTODIR          !! OUTPUT RIVER STORAGE  DIRECTORY
CHARACTER(LEN=256)     ::  CRIVVELDIR          !! OUTPUT RIVER VELOCITY DIRECTORY
CHARACTER(LEN=256)     ::  CRIVDPHDIR          !! OUTPUT RIVER DEPTH    DIRECTORY

CHARACTER(LEN=256)     ::  CFLDOUTDIR          !! OUTPUT RIVER OUTFLOW  DIRECTORY
CHARACTER(LEN=256)     ::  CFLDSTODIR          !! OUTPUT FLOODPLAIN STORAGE  DIRECTORY
CHARACTER(LEN=256)     ::  CFLDDPHDIR          !! OUTPUT FLOODPLAIN DEPTH    DIRECTORY
CHARACTER(LEN=256)     ::  CFLDFRCDIR          !! OUTPUT FLOODPLAIN FRACTION DIRECTORY
CHARACTER(LEN=256)     ::  CFLDAREDIR          !! OUTPUT FLOODPLAIN AREA     DIRECTORY

CHARACTER(LEN=256)     ::  CSFCELVDIR          !! OUTPUT SURFACE ELEVATION   DIRECTORY
CHARACTER(LEN=256)     ::  COUTFLWDIR          !! OUTPUT RIVER+FLOODP. FLOW  DIRECTORY
CHARACTER(LEN=256)     ::  CSTORGEDIR          !! OUTPUT RIVER+FLOODP. STO.  DIRECTORY
CHARACTER(LEN=256)     ::  COUTINSDIR          !! OUTPUT RIVER OUTFLOW  INSTA. DIRECTORY

CHARACTER(LEN=256)     ::  CPTHOUTDIR          !! OUTPUT PATHWAY NET OUTFLOW DIRECTORY
CHARACTER(LEN=256)     ::  CPTHFLWDIR          !! OUTPUT FLOOD PATHWAY FLOW  DIRECTORY

LOGICAL                ::  LOUTVEC             !! TRUE FOR VECTORIAL OUTPUT, FALSE FOR NX,NY OUTPUT
LOGICAL                ::  LOUTCDF             !! true for netcdf outptu false for binary
LOGICAL                ::  LOUTYYYY            !! yearly output file (netCDF)

!!=== NAMELIST/NCONF/
INTEGER(KIND=JPIM)     ::  NX                  !! NUMBER OF GRIDS IN HORIZONTAL
INTEGER(KIND=JPIM)     ::  NY                  !! NUMBER OF GRIDS IN VERTICAL
INTEGER(KIND=JPIM)     ::  NLFP                !! NUMBER OF VERTICAL LEVELS DEFINING FLOODPLAIN 
REAL(KIND=JPRB)        ::  DT                  !! TIME STEP LENGTH [SEC]
REAL(KIND=JPRB)        ::  DT_DEF              !! Default TIME STEP LENGTH [SEC] for adaptive time step
REAL(KIND=JPRB)        ::  WEST                !! west, east, north, south edge of the domain [deg]
REAL(KIND=JPRB)        ::  EAST
REAL(KIND=JPRB)        ::  NORTH
REAL(KIND=JPRB)        ::  SOUTH

INTEGER(KIND=JPIM)     ::  NXIN                !! NUMBER OF GRIDS IN HORIZONTAL
INTEGER(KIND=JPIM)     ::  NYIN                !! NUMBER OF GRIDS IN VERTICAL
INTEGER(KIND=JPIM)     ::  INPN                !! MAX INPUT NUMBER
REAL(KIND=JPRB)        ::  DTIN                !! SECOND IN INPUT TIME STEP [SEC]
REAL(KIND=JPRB)        ::  DROFUNIT            !! UNIT CONVERSION COEFFICIENT (1.D-3 when input [mm] is converted to [m3/m2])
LOGICAL                ::  LFLD                !! true  floodplain active
LOGICAL                ::  LKINE               !! true: use kinematic wave
LOGICAL                ::  LADPSTP             !! true  use adaptive time step
LOGICAL                ::  LFLDOUT             !! true  floodplain flow active
LOGICAL                ::  LPTHOUT             !! true  flood path flow active
LOGICAL                ::  LDAM                !! true  for dam operation
LOGICAL                ::  LMAPEND             !! true  for map data endian conversion
LOGICAL                ::  LINPEND             !! true  for input    endian conversion
LOGICAL                ::  LLEAPYR             !! true: leap year calculation. false: 365 days / year

!!=== NAMELIST/PARAM/
REAL(KIND=JPRB)        ::  PMANRIV              !! manning roughness (river)
REAL(KIND=JPRB)        ::  PMANFLD              !! manning roughness (floodplain)
REAL(KIND=JPRB)        ::  PGRV                 !! gravity acceleration [m/s2]
REAL(KIND=JPRB)        ::  PDSTMTH              !! downstream distance at river mouth [m]
REAL(KIND=JPRB)        ::  PCADP                !! CFL coefficient
REAL(KIND=JPRB)        ::  PMINSLP              !! minimum topographic slope (kinematic wave) [m/m]

INTEGER(KIND=JPIM)     ::  IMIS                !! integer undefined
REAL(KIND=JPRM)        ::  RMIS                !! real    undefined
REAL(KIND=JPRB)        ::  DMIS                !! double  undefined
PARAMETER                 (IMIS=-9999_JPIM)
PARAMETER                 (RMIS=1.E20_JPRM)
PARAMETER                 (DMIS=1.E20_JPRB)



NAMELIST/NRUNVER/  IRESTART,CRESTDIR,CRESTSTO,LRESTCDF,LSTOONLY,RESTFREQ
NAMELIST/NSIMTIME/ ISYEAR,  ISMON,   ISDAY,   IEYEAR,  IEMON,   IEDAY
NAMELIST/NMAP/     CDIMINFO,CNEXTXY, CGRAREA, CELEVTN, CNXTDST, CRIVWTH, CRIVLEN, CRIVHGT, CFLDHGT,&
                            CPTHOUT, CRIVCLINC,CRIVPARNC, LMAPCDF,&
                            LDIKHGT, CDIKHGT
NAMELIST/NINPUT/   LINTERP, CINPMAT, LINPCDF, CRUNOFFDIR,CRUNOFFPRE,CRUNOFFSUF,CRUNOFFCDF,LINTERPCDF,&
                            LMEANSL, CMEANSL, LBOUNDSL, CBOUNDDIR, CBOUNDPRE, CBOUNDSUF,&
                            SYEARIN, SMONIN,   SDAYIN, CROFCDFVAR
NAMELIST/NOUTPUT/  COUTDIR, CRIVOUTDIR, CRIVSTODIR, CRIVDPHDIR, CRIVVELDIR,&
                            CFLDOUTDIR, CFLDSTODIR, CFLDDPHDIR, CFLDAREDIR, CFLDFRCDIR,&
                            CSFCELVDIR, COUTFLWDIR, CSTORGEDIR, CPTHOUTDIR, CPTHFLWDIR,&
                            COUTINSDIR, LOUTCDF,    LOUTVEC,    LOUTYYYY
NAMELIST/NCONF/    NX, NY,  NLFP, DT,  NXIN, NYIN, DTIN, DROFUNIT, INPN,WEST,EAST,NORTH,SOUTH,&
                   LFLD,    LKINE,  LADPSTP, LFLDOUT, LPTHOUT, LDAM,&
                   LMAPEND, LINPEND,LLEAPYR
NAMELIST/NPARAM/   PMANRIV, PMANFLD,  PGRV,   PDSTMTH, PCADP,  PMINSLP

END MODULE MOD_INPUT
