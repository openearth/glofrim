MODULE INIT_TIME_MOD
CONTAINS 
SUBROUTINE INIT_TIME
!**** *INIT_COND*  - INIT TIME

!     Purpose.
!     --------
!         .

!***  Interface.
!     ----------
!        *CALL* *INIT_TIME

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!        EMANUEL DUTRA

!     Modifications.
!     --------------
!        Original : 2011-05-04
!     ------------------------------------------------------------------

USE PARKIND1   ,ONLY: JPIM, JPRB, JPRM
USE MOD_INPUT  ,ONLY: LOGNAM,&
                      ISYEAR,ISMON,ISDAY,&
                      IEYEAR,IEMON,IEDAY,&
                      DT
USE MOD_TIME   ,ONLY: KSTEP,NSTEPS,KYYYYMMDD,KHHMM,KMIN,KMINSTART

USE LIB_DATES  ,ONLY: DATE2MIN
IMPLICIT NONE

INTEGER(KIND=JPIM)     :: SYYYYMMDD,SHHMM,EYYYYMMDD,EHHMM
INTEGER(KIND=JPIM)     :: EMIN

!!===========================
WRITE(LOGNAM,*) " "
WRITE(LOGNAM,*) "****INIT_TIME: ****"

KSTEP=0 

SYYYYMMDD=ISYEAR*10000+ISMON*100+ISDAY
SHHMM=0_JPIM
EYYYYMMDD=IEYEAR*10000+IEMON*100+IEDAY
EHHMM=0_JPIM

KMIN=DATE2MIN(SYYYYMMDD,SHHMM)
EMIN=DATE2MIN(EYYYYMMDD,EHHMM)
NSTEPS=int ( ( (EMIN-KMIN)*60_JPIM )/DT )

WRITE(LOGNAM,*) 'START DATE:',SYYYYMMDD,SHHMM,KMIN
WRITE(LOGNAM,*) 'END DATE:',EYYYYMMDD,EHHMM,EMIN
WRITE(LOGNAM,*) 'NSTEPS:',NSTEPS

KYYYYMMDD=SYYYYMMDD
KHHMM=SHHMM
KMINSTART=KMIN


END SUBROUTINE INIT_TIME
END MODULE INIT_TIME_MOD
