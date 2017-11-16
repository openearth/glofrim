MODULE PARKIND1
  !
  !     *** Define usual kinds for strong typing ***
  !
  use iso_c_binding
  IMPLICIT NONE
  SAVE
  !
  !     Integer Kinds
  !     -------------
  !
  INTEGER, PARAMETER :: JPIT = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: JPIS = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: JPIM = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: JPIB = SELECTED_INT_KIND(12)

  !
  !     Real Kinds
  !     ----------
  !
  INTEGER, PARAMETER :: JPRT = SELECTED_REAL_KIND(2,1)
  INTEGER, PARAMETER :: JPRS = SELECTED_REAL_KIND(4,2)
  ! equivalent to SELECTED_REAL_KIND(6,37)
  INTEGER, PARAMETER :: JPRM = c_float
  ! equivalent to SELECTED_REAL_KIND(13,300)
  INTEGER, PARAMETER :: JPRB = c_double
  !
END MODULE PARKIND1
