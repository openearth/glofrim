MODULE MOD_DIAG
USE PARKIND1 ,ONLY: JPIM, JPRB, JPRM
IMPLICIT NONE

SAVE

!! Exposed variables for BMI
REAL(KIND=JPRB),ALLOCATABLE, TARGET  ::  R2RIVINF(:,:)         !! river      inflow   [m3/s] (from upstream)
REAL(KIND=JPRB),ALLOCATABLE, TARGET  ::  R2RIVDPH(:,:)         !! river      depth    [m]
REAL(KIND=JPRB),ALLOCATABLE, TARGET  ::  R2RIVVEL(:,:)         !! flow velocity       [m/s]

REAL(KIND=JPRB),ALLOCATABLE, TARGET  ::  R2FLDINF(:,:)         !! floodplain inflow   [m3/s]
REAL(KIND=JPRB),ALLOCATABLE, TARGET  ::  R2FLDDPH(:,:)         !! floodplain depth    [m]
REAL(KIND=JPRB),ALLOCATABLE, TARGET  ::  R2FLDFRC(:,:)         !! flooded    fractipn [m2/m2]
REAL(KIND=JPRB),ALLOCATABLE, TARGET  ::  R2FLDARE(:,:)         !! flooded    area     [m2]

REAL(KIND=JPRB),ALLOCATABLE, TARGET  ::  R2PTHOUT(:,:)         !! flood path outflow   [m3/s]
REAL(KIND=JPRB),ALLOCATABLE, TARGET  ::  R2PTHINF(:,:)         !! flood path inflow   [m3/s]

REAL(KIND=JPRB),ALLOCATABLE, TARGET  ::  R2SFCELV(:,:)         !! water surface elev  [m]    (elevtn - rivhgt + rivdph)
REAL(KIND=JPRB),ALLOCATABLE, TARGET  ::  R2OUTFLW(:,:)         !! total outflow       [m3/s] (rivout + fldout)
REAL(KIND=JPRB),ALLOCATABLE, TARGET  ::  R2STORGE(:,:)         !! total storage       [m3]   (rivsto + fldsto)

REAL(KIND=JPRB),ALLOCATABLE, TARGET  ::  R2RIVOUT_AVG(:,:)      !! average river       discharge
REAL(KIND=JPRB),ALLOCATABLE, TARGET  ::  R2OUTFLW_AVG(:,:)      !! average total outflow       [m3/s] (rivout + fldout)  !! bugfix v362
REAL(KIND=JPRB),ALLOCATABLE, TARGET  ::  R2FLDOUT_AVG(:,:)      !! average floodplain  discharge
REAL(KIND=JPRB),ALLOCATABLE, TARGET  ::  R2RIVVEL_AVG(:,:)      !! average flow velocity
REAL(KIND=JPRB),ALLOCATABLE, TARGET  ::  R2PTHOUT_AVG(:,:)      !! flood pathway net outflow (2D)
REAL(KIND=JPRB),ALLOCATABLE, TARGET  ::  R1PTHFLW_AVG(:,:) !! flood pathway flow (1D)

!! diagnostics / state variables initial conditions

REAL(KIND=JPRB),ALLOCATABLE  ::  D2RIVINF(:,:)         !! river      inflow   [m3/s] (from upstream)
REAL(KIND=JPRB),ALLOCATABLE  ::  D2RIVDPH(:,:)         !! river      depth    [m]
REAL(KIND=JPRB),ALLOCATABLE  ::  D2RIVVEL(:,:)         !! flow velocity       [m/s]

REAL(KIND=JPRB),ALLOCATABLE  ::  D2FLDINF(:,:)         !! floodplain inflow   [m3/s]
REAL(KIND=JPRB),ALLOCATABLE  ::  D2FLDDPH(:,:)         !! floodplain depth    [m]
REAL(KIND=JPRB),ALLOCATABLE  ::  D2FLDFRC(:,:)         !! flooded    fractipn [m2/m2]
REAL(KIND=JPRB),ALLOCATABLE  ::  D2FLDARE(:,:)         !! flooded    area     [m2]

REAL(KIND=JPRB),ALLOCATABLE  ::  D2PTHOUT(:,:)         !! flood path outflow   [m3/s]
REAL(KIND=JPRB),ALLOCATABLE  ::  D2PTHINF(:,:)         !! flood path inflow   [m3/s]

REAL(KIND=JPRB),ALLOCATABLE  ::  D2SFCELV(:,:)         !! water surface elev  [m]    (elevtn - rivhgt + rivdph)
REAL(KIND=JPRB),ALLOCATABLE  ::  D2OUTFLW(:,:)         !! total outflow       [m3/s] (rivout + fldout)
REAL(KIND=JPRB),ALLOCATABLE  ::  D2STORGE(:,:)         !! total storage       [m3]   (rivsto + fldsto)

REAL(KIND=JPRB),ALLOCATABLE  ::  D2RIVOUT_AVG(:,:)      !! average river       discharge
REAL(KIND=JPRB),ALLOCATABLE  ::  D2OUTFLW_AVG(:,:)      !! average total outflow       [m3/s] (rivout + fldout)  !! bugfix v362
REAL(KIND=JPRB),ALLOCATABLE  ::  D2FLDOUT_AVG(:,:)      !! average floodplain  discharge
REAL(KIND=JPRB),ALLOCATABLE  ::  D2RIVVEL_AVG(:,:)      !! average flow velocity
REAL(KIND=JPRB),ALLOCATABLE  ::  D2PTHOUT_AVG(:,:)      !! flood pathway net outflow (2D)
REAL(KIND=JPRB),ALLOCATABLE  ::  D1PTHFLW_AVG(:,:) !! flood pathway flow (1D)

INTEGER(KIND=JPIM)           :: NADD                   !! count to average

REAL(KIND=JPRB),ALLOCATABLE  :: D2OUTINS(:,:)          !! instantaneous discharge [m3/s] (unrouted runoff)

REAL(KIND=JPRB)  ::  DGLBRIVSTO          !! global river storage      [m3]
REAL(KIND=JPRB)  ::  DGLBFLDSTO          !! global floodplain storage [m3]
REAL(KIND=JPRB)  ::  DGLBSTOPRE          !! global water storage      [m3] (befre flow calculation)
REAL(KIND=JPRB)  ::  DGLBSTONXT          !! global water storage      [m3] (after flow calculation)
REAL(KIND=JPRB)  ::  DGLBSTONEW          !! global water storage      [m3] (after runoff input)
REAL(KIND=JPRB)  ::  DGLBRIVINF          !! global inflow             [m3] (rivinf + fldinf)
REAL(KIND=JPRB)  ::  DGLBRIVOUT          !! global outflow            [m3] (rivout + fldout)
REAL(KIND=JPRB)  ::  DGLBFLDARE          !! global flooded area       [m2]

END MODULE MOD_DIAG
