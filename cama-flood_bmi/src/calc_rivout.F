MODULE CALC_RIVOUT_MOD
CONTAINS
SUBROUTINE CALC_RIVOUT
! ================================================
! to   Calculate discharge, Inertial Equation
! ================================================
USE PARKIND1   ,ONLY: JPIM, JPRB
USE MOD_INPUT  ,ONLY: DT
USE MOD_INPUT  ,ONLY: PMANRIV,  PDSTMTH,  PGRV
USE MOD_INPUT  ,ONLY: LMEANSL, LBOUNDSL
USE MOD_MAP    ,ONLY: I1NEXT,   NSEQALL,  NSEQRIV,  NSEQMAX
USE MOD_MAP    ,ONLY: D2RIVELV, D2ELEVTN, D2NXTDST, D2RIVWTH
USE MOD_MAP    ,ONLY: D2MEANSL
USE MOD_PROG   ,ONLY: D2RIVSTO_OUT, D2RIVOUT_OUT
USE MOD_PROG   ,ONLY: D2RIVOUT_PRE, D2RIVDPH_PRE
USE MOD_PROG   ,ONLY: D2RIVSTO_PRE
USE MOD_PROG   ,ONLY: D2BOUND 
USE MOD_DIAG   ,ONLY: D2RIVDPH, D2RIVVEL, D2RIVINF, D2SFCELV

implicit none

!*** Local
      REAL(KIND=JPRB)    ::  D2SFCELV_PRE(NSEQMAX,1)                  !! water surface elevation (t-1) [m]
      REAL(KIND=JPRB)    ::  D2STOOUT(NSEQMAX,1)                      !! total outflow from a grid     [m3]
      REAL(KIND=JPRB)    ::  D2RATE(NSEQMAX,1)                        !! outflow correction
!$ SAVE
      INTEGER(KIND=JPIM) ::  ISEQ, JSEQ
      REAL(KIND=JPRB)    ::  DSLOPE, DAREA, DFLW, DOUT_PRE, DFLW_PRE, DFLW_IMP
!$OMP THREADPRIVATE   (JSEQ, DSLOPE, DAREA, DFLW, DOUT_PRE, DFLW_PRE, DFLW_IMP)
! ================================================
!$OMP PARALLEL DO
DO ISEQ=1, NSEQALL
  D2SFCELV(ISEQ,1)     = D2RIVELV(ISEQ,1) + D2RIVDPH(ISEQ,1)
  D2SFCELV_PRE(ISEQ,1) = D2RIVELV(ISEQ,1) + D2RIVDPH_PRE(ISEQ,1)
  D2RIVINF(ISEQ,1) = 0.D0
  D2STOOUT(ISEQ,1) = 0.D0
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
DO ISEQ=1, NSEQRIV                                                    !! for normal cells
  JSEQ=I1NEXT(ISEQ)

  DSLOPE =  (D2SFCELV(ISEQ,1)-D2SFCELV(JSEQ,1)) * D2NXTDST(ISEQ,1)**(-1.D0)
  DFLW = MAX(D2SFCELV(ISEQ,1),D2SFCELV(JSEQ,1)) - D2RIVELV(ISEQ,1)
  DAREA = D2RIVWTH(ISEQ,1) * DFLW

  DFLW_PRE=MAX(D2SFCELV_PRE(ISEQ,1),D2SFCELV_PRE(JSEQ,1)) - D2RIVELV(ISEQ,1)
  DFLW_IMP=(DFLW*DFLW_PRE)**0.5D0                                       !! semi implicit flow depth
  IF( DFLW_IMP<=0.D0 ) DFLW_IMP=DFLW

  IF( DFLW_IMP>1.D-5 .and. DAREA>1.D-5 )THEN                                                  !! local inertial equation [Bates et al., 2010, J.Hydro.]
    DOUT_PRE = D2RIVOUT_PRE(ISEQ,1) * D2RIVWTH(ISEQ,1)**(-1.D0)                         !! outflow (t-1) [m2/s] (unit width)
    D2RIVOUT_OUT(ISEQ,1) = D2RIVWTH(ISEQ,1) * ( DOUT_PRE + PGRV*DT*DFLW_IMP*DSLOPE ) &
                           * ( 1.D0 + PGRV*DT*PMANRIV**2.D0*abs(DOUT_PRE)*DFLW_IMP**(-7.D0/3.D0) )**(-1.D0)
    D2RIVVEL(ISEQ,1) = D2RIVOUT_OUT(ISEQ,1) * DAREA**(-1.D0)
  ELSE
    D2RIVOUT_OUT(ISEQ,1) = 0.D0
    D2RIVVEL(ISEQ,1) = 0.D0
  ENDIF

  IF( D2RIVOUT_OUT(ISEQ,1) >= 0.D0 )THEN                                  !! total outflow from a grid
    D2STOOUT(ISEQ,1) = D2STOOUT(ISEQ,1) + D2RIVOUT_OUT(ISEQ,1)*DT
  ELSE
!$OMP ATOMIC
    D2STOOUT(JSEQ,1) = D2STOOUT(JSEQ,1) - D2RIVOUT_OUT(ISEQ,1)*DT
  ENDIF
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO                                                     !! for river mouth grids
DO ISEQ=NSEQRIV+1, NSEQALL
  IF( LMEANSL ) THEN
    IF( LBOUNDSL )THEN
      DSLOPE = ( D2SFCELV(ISEQ,1)-(D2ELEVTN(ISEQ,1)+D2BOUND(ISEQ,1)+D2MEANSL(ISEQ,1)) ) * PDSTMTH**(-1.D0)   !! variable rivermouth boundary + mean sea level
    ELSE
      DSLOPE = ( D2SFCELV(ISEQ,1)-(D2ELEVTN(ISEQ,1)+D2MEANSL(ISEQ,1)) ) * PDSTMTH**(-1.D0)   !! mean sea level
    ENDIF
  ELSE
    IF( LBOUNDSL )THEN
      DSLOPE = ( D2SFCELV(ISEQ,1)-(D2ELEVTN(ISEQ,1)+D2BOUND(ISEQ,1)) ) * PDSTMTH**(-1.D0)   !! variable rivermouth boundary
    ELSE
    DSLOPE = ( D2SFCELV(ISEQ,1)-D2ELEVTN(ISEQ,1) ) * PDSTMTH**(-1.D0)
    ENDIF
  ENDIF


!!  DSLOPE = MAX( DSLOPE, 0.D0 )                                        !! no backfrom from ocean
  DFLW   = D2RIVDPH(ISEQ,1)
  DAREA  = D2RIVWTH(ISEQ,1) * DFLW

  DFLW_PRE=D2RIVDPH_PRE(ISEQ,1)
  DFLW_IMP=(DFLW*DFLW_PRE)**0.5D0                                     !! semi implicit flow depth
  IF( DFLW_IMP<=0.D0 ) DFLW_IMP=DFLW

  IF( DFLW_IMP>1.D-5 .and. DAREA>1.D-5 )THEN 
    DOUT_PRE = D2RIVOUT_PRE(ISEQ,1) * D2RIVWTH(ISEQ,1)**(-1.D0)
    D2RIVOUT_OUT(ISEQ,1) = D2RIVWTH(ISEQ,1) * ( DOUT_PRE + PGRV*DT*DFLW_IMP*DSLOPE ) &
                           * ( 1.D0 + PGRV*DT*PMANRIV**2.D0*abs(DOUT_PRE)*DFLW_IMP**(-7.D0/3.D0) )**(-1.D0)
    D2RIVVEL(ISEQ,1) = D2RIVOUT_OUT(ISEQ,1) * DAREA**(-1.D0)
  ELSE
    D2RIVOUT_OUT(ISEQ,1) = 0.D0
    D2RIVVEL(ISEQ,1) = 0.D0
  ENDIF
  IF( D2RIVOUT_OUT(ISEQ,1)>=0.D0 ) D2STOOUT(ISEQ,1) = D2STOOUT(ISEQ,1) + D2RIVOUT_OUT(ISEQ,1)*DT
END DO
!$OMP END PARALLEL DO


!$OMP PARALLEL DO                                                     !! outflow correcttion if total outflow > storage
DO ISEQ=1, NSEQALL
  IF( D2STOOUT(ISEQ,1) > 0.D0)THEN
    D2RATE(ISEQ,1) = MIN( D2RIVSTO_OUT(ISEQ,1) * D2STOOUT(ISEQ,1)**(-1.D0), 1.D0 )
  ELSE
    D2RATE(ISEQ,1) = 1.D0
  ENDIF
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
DO ISEQ=1, NSEQRIV
  JSEQ=I1NEXT(ISEQ)
  IF( D2RIVOUT_OUT(ISEQ,1) >= 0.D0 )THEN
    D2RIVOUT_OUT(ISEQ,1) = D2RIVOUT_OUT(ISEQ,1)*D2RATE(ISEQ,1)
  ELSE
    D2RIVOUT_OUT(ISEQ,1) = D2RIVOUT_OUT(ISEQ,1)*D2RATE(JSEQ,1)
  ENDIF
!$OMP ATOMIC
  D2RIVINF(JSEQ,1) = D2RIVINF(JSEQ,1) + D2RIVOUT_OUT(ISEQ,1)             !! total inflow to a grid (from upstream)

  D2RIVOUT_PRE(ISEQ,1)=D2RIVOUT_OUT(ISEQ,1)                              !! save outflow (t)
  D2RIVSTO_PRE(ISEQ,1)=D2RIVSTO_OUT(ISEQ,1)
  D2RIVDPH_PRE(ISEQ,1)=D2RIVDPH(ISEQ,1)                              !! save depth   (t)
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
DO ISEQ=NSEQRIV+1, NSEQALL
  D2RIVOUT_OUT(ISEQ,1) = D2RIVOUT_OUT(ISEQ,1)*D2RATE(ISEQ,1)
  D2RIVOUT_PRE(ISEQ,1)=D2RIVOUT_OUT(ISEQ,1)                              !! save outflow (t)
  D2RIVSTO_PRE(ISEQ,1)=D2RIVSTO_OUT(ISEQ,1)
  D2RIVDPH_PRE(ISEQ,1)=D2RIVDPH(ISEQ,1)                              !! save depth   (t)
END DO
!$OMP END PARALLEL DO

END SUBROUTINE CALC_RIVOUT
END MODULE CALC_RIVOUT_MOD
