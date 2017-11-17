PROGRAM MAIN_day                                           !! Main program for CaMa-Flood

USE PARKIND1      ,ONLY: JPIM
USE MOD_INPUT     ,ONLY: ISYEAR
USE MOD_MAP       ,ONLY: REGIONALL, REGIONTHIS
USE MOD_OUTPUT    ,ONLY: CSUFBIN, CSUFVEC, CSUFCDF
USE CONTROL0_MOD

IMPLICIT NONE

#ifdef UseMPI
INCLUDE 'mpif.h'
INTEGER(KIND=JPIM)    :: ierr, Nproc, Nid
#endif
!!==================================================


#ifdef UseMPI
CALL MPI_Init(ierr)
CALL MPI_Comm_size(MPI_COMM_WORLD, Nproc, ierr)
CALL MPI_Comm_rank(MPI_COMM_WORLD, Nid, ierr)
REGIONALL =Nproc
REGIONTHIS=Nid+1
#endif

IF (REGIONALL>=2 )then                                      !! Regional Output for MPI run
  WRITE(CSUFBIN,'(a5,i2.2)') '.bin-', REGIONTHIS            !! Change suffix of output file for each calculation node
  WRITE(CSUFVEC,'(a5,i2.2)') '.vec-', REGIONTHIS
  WRITE(CSUFCDF,'(a4,i2.2)') '.nc-',  REGIONTHIS
ENDIF

CALL CONTROL0

#ifdef UseMPI
CALL MPI_Finalize(ierr)
#endif

!!==================================================

END PROGRAM MAIN_day
