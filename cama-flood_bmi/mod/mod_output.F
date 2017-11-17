MODULE MOD_OUTPUT
USE PARKIND1 ,ONLY: JPIM, JPRB
IMPLICIT NONE

!**** OUTPUT RELATED 

SAVE 


INTEGER(KIND=JPIM)   ::  NVARS                    !! variable number for netCDF output 
PARAMETER               (NVARS=15)
INTEGER(KIND=JPIM)   ::  FUNITS(NVARS)            !! file number for binary output
INTEGER(KIND=JPIM)   ::  CDFUNITS(NVARS)          !! file number for netCDF output

CHARACTER(LEN=256)   ::  COUTVARDIR(NVARS)        !! output file directory
CHARACTER(LEN=128)   ::  COUTVN(NVARS)            !! output variable name
CHARACTER(LEN=128)   ::  COUTLN(NVARS)            !! output variable long name
CHARACTER(LEN=128)   ::  COUTUNIT(NVARS)          !! output variable unit

CHARACTER(LEN=128)   ::  CVNRESTART

LOGICAL              ::  LOPEN

CHARACTER(LEN=64)    ::  CSUFBIN             !! SUFFIX FOR GRID   BINARY
CHARACTER(LEN=64)    ::  CSUFVEC             !! SUFFIX FOR VECTOR BINARY
CHARACTER(LEN=64)    ::  CSUFCDF             !! SUFFIX FOR netCDF
CHARACTER(LEN=64)    ::  CSUFPTH             !! SUFFIX FOR PATHWAY

DATA FUNITS       /101,102,103,104,  105,106,107,108,109,  110,111,112,113,  114,115/

DATA COUTVN  /'rivout','rivsto','rivdph','rivvel',  'fldout','fldsto','flddph','fldfrc','fldare', &
                   'sfcelv','outflw','storge','outins',  'pthout','pthflw'/

DATA COUTLN /'river discharge','river storage','river depth','river velocity', &
                   'flood discharge','flood storage','flood depth','flood fraction','flood area', &
                   'water surface elevation','discharge (rivout+fldout)','storge (rivsto+fldsto)','inst discharge (not routed)', &
                   'net bifurcation discharge','bifurcation discharge (channel)'/

DATA COUTUNIT     /'m3/s','m3','m','m/s',  'm3/s','m3','m','m2/m2','m2', &
                   'm','m3/s','m3','m3/s', 'm3/s','m3/s'/

DATA CVNRESTART  /'restart'/

DATA LOPEN   /.FALSE./

DATA csufbin /'.bin'/
DATA csufvec /'.vec'/
DATA csufcdf /'.nc'/
DATA csufpth /'.pth'/

END MODULE MOD_OUTPUT