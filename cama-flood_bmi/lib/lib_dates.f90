

MODULE LIB_DATES
USE PARKIND1   ,ONLY: JPIM, JPRB, JPRM, JPIB
USE MOD_INPUT  ,ONLY: LLEAPYR

IMPLICIT NONE

INTEGER(KIND=JPIM)    :: YYYY0, MM0, DD0
PARAMETER               (YYYY0=1850)
PARAMETER               (MM0=1)
PARAMETER               (DD0=1)

CONTAINS

!!==================================================
SUBROUTINE MIN2DATE(INMI,YYYYMMDD,HHMM)     !!  Return YYYYMMDD and HHMM for INMI
USE MOD_INPUT  ,ONLY: NULNAM

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: INMI      !!  input minutes
INTEGER(KIND=JPIM),INTENT(OUT) :: YYYYMMDD
INTEGER(KIND=JPIM),INTENT(OUT) :: HHMM

INTEGER(KIND=JPIM)             :: YYYY,MM,DD,HH,MI,NDAYS,NDM,ID
!!====================
YYYYMMDD = 0
HHMM     = 0

NDAYS = INMI/1440              !! days  in INMI : 1440 = (minutes in a day)
MI    = MOD(INMI,1440)
HH    = MI/60                  !! hours in INMI
MI    = MOD(MI,60)             !! mins  in INMI

YYYY     = YYYY0
MM       = MM0
DD       = DD0
NDM      = IMDAYS(YYYY,MM)     !! number of days in a month

! PRINT*,YYYY,MM,DD
DO ID=1,NDAYS
  DD=DD+1
  IF ( DD .GT. NDM ) THEN
    MM=MM+1
    DD=1
    IF ( MM .GT. 12 ) THEN
      MM=1
      YYYY=YYYY+1
    ENDIF
    NDM=IMDAYS(YYYY,MM)
  ENDIF
ENDDO

HHMM     = HH*100+MI
YYYYMMDD = YYYY*10000+MM*100+DD

! PRINT *, 'MIN2DATE====='
! PRINT *, INMI,YYYYMMDD,HHMM
! PRINT *, NDAYS,HH,MI
! PRINT *, YYYY,MM,DD

END SUBROUTINE MIN2DATE


!!==================================================
FUNCTION DATE2MIN(YYYYMMDD,HHMM)
USE MOD_INPUT  ,ONLY: NULNAM
IMPLICIT NONE

INTEGER(KIND=JPIM)            :: DATE2MIN

INTEGER(KIND=JPIM),INTENT(IN) :: YYYYMMDD
INTEGER(KIND=JPIM),INTENT(IN) :: HHMM

INTEGER(KIND=JPIM)            :: YYYY,MM,DD,HH,MI,D2MIN
INTEGER(KIND=JPIM)            :: IY,IM

!!====================
D2MIN    = 1440
DATE2MIN = 0
CALL SPLITDATE(YYYYMMDD,YYYY,MM,DD)
HH = HHMM/100
MI = HHMM-HH*100

! PRINT*,'DATE2MIN====='
! PRINT*,YYYYMMDD,HHMM,YYYY,MM,DD,HH,MM

!!====================
!! SOME CHECKS
IF (YYYY .LT. YYYY0) THEN
  WRITE(NULNAM,*) 'DATE2MIN: YYYY .LT. YYYY0: CHANGE YYYY0 IN MOD_DATES',YYYY,YYYY0
  STOP
ENDIF
IF (MM .GT. 12) THEN
  WRITE(NULNAM,*) 'DATE2MIN: MM .GT. 12: SOME PROBLEM',MM
  STOP
ENDIF
IF (DD .GT. IMDAYS(YYYY,MM)) THEN
  WRITE(NULNAM,*) 'DATE2MIN: DD .GT. : SOME PROBLEM',DD,IMDAYS(YYYY,MM)
  STOP
ENDIF
IF (HH .GT. 24) THEN
  WRITE(NULNAM,*) 'DATE2MIN: HH .GT. 24 : SOME PROBLEM',HH
  STOP
ENDIF
IF (MI .GT. 60) THEN
  WRITE(NULNAM,*) 'DATE2MIN: MI .GT. 60 : SOME PROBLEM',MI
  STOP
ENDIF

IY=YYYY0
DO WHILE (IY .LT. YYYY)
  DO IM=1,12
    DATE2MIN=DATE2MIN+IMDAYS(IY,IM)*D2MIN
  ENDDO
  IY=IY+1
ENDDO
IM=1
DO WHILE (IM .LT. MM )
  DATE2MIN=DATE2MIN+IMDAYS(IY,IM)*D2MIN
  IM=IM+1
ENDDO

DATE2MIN = DATE2MIN + (DD-1)*D2MIN
DATE2MIN = DATE2MIN + HH*60 + MI

END FUNCTION DATE2MIN


! =================================================
SUBROUTINE SPLITDATE(YYYYMMDD,YYYY,MM,DD)
IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: YYYYMMDD
INTEGER(KIND=JPIM),INTENT(OUT) :: YYYY,MM,DD

!!====================
YYYY =  YYYYMMDD/10000
MM   = (YYYYMMDD-YYYY*10000)/100
DD   =  YYYYMMDD-(YYYY*10000+MM*100)

END SUBROUTINE SPLITDATE


! ==================================================
FUNCTION IMDAYS(IYEAR,IMON)   !! days in month
IMPLICIT NONE
 
INTEGER(KIND=JPIM)            :: IMDAYS

INTEGER(KIND=JPIM),INTENT(IN) :: IYEAR
INTEGER(KIND=JPIM),INTENT(IN) :: IMON

INTEGER(KIND=JPIM)            :: ND(12)
DATA ND /31,28,31,30,31,30,31,31,30,31,30,31/

!!====================
IMDAYS=ND(IMON)
IF ( IMON == 2 .and. LLEAPYR ) THEN
  IF ( MOD(IYEAR,400) == 0 .OR. (MOD(IYEAR,100) .NE. 0 .AND. MOD(IYEAR,4) .EQ. 0 )) IMDAYS=29
ENDIF

END FUNCTION IMDAYS


END MODULE LIB_DATES