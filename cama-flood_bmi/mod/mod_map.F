MODULE MOD_MAP
USE PARKIND1 ,ONLY: JPIM, JPRB
IMPLICIT NONE

!**** MAP fields

SAVE 

INTEGER(KIND=JPIM),ALLOCATABLE           ::  I2NEXTX(:,:)       !! POINT DOWNSTREAM HORIZONTAL
INTEGER(KIND=JPIM),ALLOCATABLE           ::  I2NEXTY(:,:)       !! POINT DOWNSTREAM VERTICAL
INTEGER(KIND=JPIM),ALLOCATABLE           ::  I2RIVSEQ(:,:)      !! river sequence map:# of cellls from the hilltop cell

INTEGER(KIND=JPIM),ALLOCATABLE           ::  I1SEQX(:)          !! 1D SEQUENCE HORIZONTAL
INTEGER(KIND=JPIM),ALLOCATABLE           ::  I1SEQY(:)          !! 1D SEQUENCE VERTICAL
INTEGER(KIND=JPIM),ALLOCATABLE           ::  I1NEXT(:)          !! 1D DOWNSTREAM
INTEGER(KIND=JPIM)                       ::  NSEQRIV            !! LENGTH OF 1D SEQUNECE FOR RIVER
INTEGER(KIND=JPIM)                       ::  NSEQALL            !! LENGTH OF 1D SEQUNECE FOR RIVER AND MOUTH
INTEGER(KIND=JPIM)                       ::  NSEQMAX            !! MAX OF NSEQALL (PARALLEL)

INTEGER(KIND=JPIM),ALLOCATABLE           ::  I2VECTOR(:,:)      !! VECTOR INDEX
INTEGER(KIND=JPIM),ALLOCATABLE           ::  I2REGION(:,:)      !! REGION INDEX
INTEGER(KIND=JPIM)                       ::  REGIONALL          !! REGION TOTAL
INTEGER(KIND=JPIM)                       ::  REGIONTHIS         !! REGION THIS CPU

INTEGER(KIND=JPIM)                       ::  RIVSEQMAX          !! MAX. OF CATCH. IN A RIVER

REAL(KIND=JPRB),ALLOCATABLE              ::  D1LON(:)           !! longitude [degree_east]
REAL(KIND=JPRB),ALLOCATABLE              ::  D1LAT(:)           !! latitude  [degree_north]

REAL(KIND=JPRB),ALLOCATABLE              ::  D2GRAREA(:,:)      !! GRID AREA [M2]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2ELEVTN(:,:)      !! ELEVATION [M]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2NXTDST(:,:)      !! DISTANCE TO THE NEXT GRID [M]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2RIVLEN(:,:)      !! RIVER LENGTH [M]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2RIVWTH(:,:)      !! RIVER WIDTH [M]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2RIVHGT(:,:)      !! RIVER HEIGHT [M]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2FLDHGT(:,:,:)    !! FLOODPLAIN HEIGHT [M]

REAL(KIND=JPRB),ALLOCATABLE              ::  D2RIVSTOMAX(:,:)   !! maximum river storage [m3]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2RIVELV(:,:)      !! elevation of river bed [m3]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2FLDSTOMAX(:,:,:) !! MAXIMUM FLOODPLAIN STORAGE [M3]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2FLDGRD(:,:,:)    !! FLOODPLAIN GRADIENT
REAL(KIND=JPRB)                          ::  DFRCINC            !! FLOODPLAIN FRACTION INCREMENT [-] (1/NLFP)

REAL(KIND=JPRB),ALLOCATABLE              ::  D2MEANSL(:,:)      !! MEAN SEA LEVEL

REAL(KIND=JPRB),ALLOCATABLE              ::  D2FLDSTORES(:,:,:) !! each floodplain storage [m3]

REAL(KIND=JPRB),ALLOCATABLE              ::  D2DIKHGT(:,:)      !! DIKE HEIGHT [m]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2DIKSTOMAX(:,:)   !! maximum dike storage [m3]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2DIKFLDSTO(:,:)   !! dike+floodplain storage [m3]

CHARACTER(LEN=256)                       ::  CINPMAT            !! INPUT MATRIX
INTEGER(KIND=JPIM),ALLOCATABLE           ::  INPX(:,:,:)        !! INPUT GRID XIN
INTEGER(KIND=JPIM),ALLOCATABLE           ::  INPY(:,:,:)        !! INPUT GRID YIN
REAL(KIND=JPRB),ALLOCATABLE              ::  INPA(:,:,:)        !! INPUT AREA

INTEGER(KIND=JPIM)                       ::  NPTHOUT            !! NUMBER OF FLOODPLAIN PATH
INTEGER(KIND=JPIM)                       ::  NPTHLEV            !! NUMBER OF FLOODPLAIN PATH LAYER
INTEGER(KIND=JPIM),ALLOCATABLE           ::  PTH_UPST(:)        !! FLOOD PATHWAY UPSTREAM   ISEQ
INTEGER(KIND=JPIM),ALLOCATABLE           ::  PTH_DOWN(:)        !! FLOOD PATHWAY DOWNSTREAM JSEQ
REAL(KIND=JPRB),ALLOCATABLE              ::  PTH_DST(:)         !! FLOOD PATHWAY DISTANCE [m]
REAL(KIND=JPRB),ALLOCATABLE              ::  PTH_ELV(:,:)         !! FLOOD PATHWAY ELEVATION [m]
REAL(KIND=JPRB),ALLOCATABLE              ::  PTH_WTH(:,:)         !! FLOOD PATHWAY WIDTH [m]
REAL(KIND=JPRB),ALLOCATABLE              ::  PTH_MAN(:)         !! FLOOD PATHWAY Manning


DATA REGIONALL  /1/
DATA REGIONTHIS /1/

END MODULE MOD_MAP
