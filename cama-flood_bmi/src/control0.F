MODULE CONTROL0_MOD
CONTAINS
SUBROUTINE CONTROL0
! ================================================
! Controls the job at level 0, the lowest level.
! ================================================
USE PARKIND1 ,ONLY: JPIM, JPRB
USE MOD_INPUT,ONLY: LOGNAM, CLOGOUT, ISYEAR
USE MOD_MAP  ,ONLY: REGIONALL, REGIONTHIS

USE INIT_INPUTNAM_MOD
USE INIT_MAP_MOD
USE INIT_TOPO_MOD
USE INIT_COND_MOD
USE INIT_TIME_MOD
USE CONTROL_TSTP_MOD

implicit none
! ================================================
! *** 1. read input namelist 
CALL INIT_INPUTNAM

IF (REGIONALL>=2 )then                                                !! for MPI
  WRITE(CLOGOUT,'(a5,i4.4,a1,i2.2,a4)') './out', ISYEAR, '-', REGIONTHIS, '.log'
  LOGNAM=12
  OPEN(LOGNAM,FILE=CLOGOUT,FORM='FORMATTED')
ENDIF

! *** 2. read input river map
CALL INIT_MAP

! *** 3. set topography
CALL INIT_TOPO

! *** 4. set initial conditions (zero-storage / restart)
CALL INIT_COND

! *** 5. set time
CALL INIT_TIME

! *** 6. time-step loop
CALL CONTROL_TSTP

IF (REGIONALL>=2 )then
  CLOSE(LOGNAM)
ENDIF

END SUBROUTINE CONTROL0

END MODULE CONTROL0_MOD
