MODULE CALC_DAMOUT_MOD
CONTAINS
SUBROUTINE CALC_DAMOUT
! ================================================
! to   Calculate discharge, Inertial Equation
! ================================================
USE PARKIND1   ,ONLY: JPIM, JPRB
USE MOD_MAP    ,ONLY: I1NEXT,   I2VECTOR
USE MOD_PROG   ,ONLY: D2RIVOUT_OUT, D2FLDOUT_OUT
USE MOD_DIAG   ,ONLY: D2RIVINF,     D2FLDINF

implicit none

!*** Local
      INTEGER(KIND=JPIM) ::  IX, IY, ISEQ, JSEQ
      REAL(KIND=JPRB)    ::  DOUTFLW, DEXCESS, DRIVRED, DFLDRED       !! Total outflw, Excess flow, river flow reduction, floodplain flow reduction

!*** DAM Parameter
      REAL(KIND=JPRB)    ::  DMAXOUT
! ================================================

!! for Stung Toreng (E105,94,N13.55)

IX=160   !! please check (ix,iy) of the target point
IY=215
ISEQ=I2VECTOR(IX,IY)  ! convert 2D-Map (ix,iy) to 1D-vector (iseq,1) 
JSEQ=I1NEXT(ISEQ)

!! dam operation (modify rivout & fldout)

DOUTFLW=D2RIVOUT_OUT(ISEQ,1)+D2FLDOUT_OUT(ISEQ,1)
DMAXOUT=40000.

IF( DOUTFLW > DMAXOUT )THEN
  DEXCESS=DOUTFLW-DMAXOUT
  IF( D2FLDOUT_OUT(ISEQ,1)>DEXCESS )THEN
    D2FLDOUT_OUT(ISEQ,1)=D2FLDOUT_OUT(1,ISEQ)-DEXCESS
    D2FLDINF(JSEQ,1)=D2FLDINF(JSEQ,1)-DEXCESS
  ELSE
    DRIVRED=D2RIVOUT_OUT(ISEQ,1)-DMAXOUT
    DFLDRED=D2FLDOUT_OUT(ISEQ,1)

    D2RIVOUT_OUT(ISEQ,1)=DMAXOUT
    D2FLDOUT_OUT(ISEQ,1)=0.0

    D2RIVINF(JSEQ,1)=D2RIVINF(JSEQ,1)-DRIVRED
    D2FLDINF(JSEQ,1)=D2FLDINF(JSEQ,1)-DFLDRED
  ENDIF
ENDIF


END SUBROUTINE CALC_DAMOUT
END MODULE CALC_DAMOUT_MOD
