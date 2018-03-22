module demo_model
  use iso_c_binding

  use iso_c_utils
  use logging

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
  USE MOD_OUTPUT    ,ONLY: CSUFBIN, CSUFVEC, CSUFCDF
  USE CONTROL0_MOD

  ! ================================================
  ! CONTROL TIME-STEP INTEGRATION
  ! ================================================
  USE PARKIND1   ,ONLY: JPIM, JPRB, JPRM
  USE MOD_INPUT  ,ONLY: LOGNAM, NX, NY, DT, LOUTCDF
  USE MOD_INPUT  ,ONLY: LBMIROF, LBOUNDSL
  USE MOD_INPUT  ,ONLY: COUTDIR
  USE MOD_TIME   ,ONLY: KSTEP,  NSTEPS, KYYYYMMDD, KHHMM, KMIN, KMINSTART, KYYYYMMDDP, KHHMMP, KMINP
  USE MOD_MAP    ,ONLY: NSEQMAX, NPTHOUT, NPTHLEV
  USE MOD_DIAG   ,ONLY: D2RIVINF, D2RIVDPH, D2RIVVEL, D2FLDINF, D2FLDDPH, D2FLDFRC, D2FLDARE,&
       D2PTHOUT, D2PTHINF, D2SFCELV, D2OUTFLW, D2STORGE
  USE MOD_DIAG   ,ONLY: D2OUTFLW_AVG, D2RIVOUT_AVG, D2FLDOUT_AVG, D2RIVVEL_AVG, D2PTHOUT_AVG, D1PTHFLW_AVG, D2OUTINS, NADD
  USE MOD_DIAG   ,ONLY: R2RIVINF, R2RIVDPH, R2RIVVEL, R2FLDINF, R2FLDDPH, R2FLDFRC, R2FLDARE,&
       R2PTHOUT, R2PTHINF, R2SFCELV, R2OUTFLW, R2STORGE
  USE MOD_DIAG   ,ONLY: R2OUTFLW_AVG, R2RIVOUT_AVG, R2FLDOUT_AVG, R2RIVVEL_AVG, R2PTHOUT_AVG, R1PTHFLW_AVG
  USE MOD_PROG   ,ONLY: D2RIVSTO_IN, D2RIVSTO_OUT
  USE MOD_PROG   ,ONLY: R2ROFFIN, R2RUNOFF, R2RUNOFF_M, R2RIVSTO_OUT, R2FLDSTO_OUT
  USE LIB_DATES  ,ONLY: MIN2DATE, SPLITDATE

  ! USE MOD_PROG   ,ONLY: D2RIVSTO_OUT, D2FLDSTO_OUT, D2RIVOUT_OUT, D2FLDOUT_OUT, D1PTHFLW_OUT
  ! USE MOD_PROG   ,ONLY: D2RIVSTO_IN, D2FLDSTO_IN, D2RIVOUT_IN , D2FLDOUT_IN, D1PTHFLW_IN

  USE CONTROL_INP_MOD
  USE CONTROL_INP_BMI_MOD
  USE CONTROL_BOUND_MOD
  USE CONTROL_PHY_MOD
  USE CONTROL_OUT_MOD
  USE CONTROL_REST_MOD

  implicit none
  INTEGER(KIND=JPIM)    ::    ISTEP                              !! local time step
  INTEGER(KIND=JPIM)    ::    IYYYY, IMM, IDD                    !! date (end of time step)
  REAL(KIND=JPRM)       ::    ELAPSED(2),TOTAL,TOTAL1,TSEC !! time elapsed related

  ! This is assumed.....
  integer(c_int), parameter :: MAXDIMS = 6

  double precision, target :: t
  double precision, target :: t_end
  double precision, target :: t_start

  double precision, target :: arr1(3)
  integer, target :: arr2(2,3)
  logical(c_bool), target :: arr3(2,2,3)
contains



  integer(c_int) function finalize() result(ierr) bind(C, name="finalize")
    !DEC$ ATTRIBUTES DLLEXPORT::finalize
    ierr = 0
#ifdef UseMPI
    CALL MPI_Finalize(ierr)
#endif

    call log(LEVEL_INFO, 'Finalize')
  end function finalize


  integer(c_int) function initialize(c_configfile) result(ierr) bind(C, name="initialize")
    !DEC$ ATTRIBUTES DLLEXPORT::initialize

    implicit none

    ! Variables
    character(kind=c_char), intent(in) :: c_configfile(*)
    character(len=strlen(c_configfile)) :: configfile
    CHARACTER(LEN=4)      :: IYEAR
    CHARACTER(LEN=256)    :: CLOGTXT


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

    ! Convert c string to fortran string
    ierr = 0
    t = 0.0d0
    t_end = 10.0d0
    configfile = char_array_to_string(c_configfile)

    write(msgbuf,*) 'Initializing with ', configfile
    call log(LEVEL_INFO, trim(msgbuf))

    ! ================================================
    ! *** 1. read input namelist
    CALL INIT_INPUTNAM

    IF (REGIONALL>=2 )then                                                !! for MPI
       WRITE(CLOGOUT,'(a5,i4.4,a1,i2.2,a4)') './out', ISYEAR, '-', REGIONTHIS, '.log'
       LOGNAM=12
       OPEN(LOGNAM,FILE=CLOGOUT,FORM='FORMATTED')
    ENDIF
    
    LOGNAM=12
    WRITE(IYEAR, '(i4)') ISYEAR
    CLOGTXT=TRIM(COUTDIR)//'run_'//TRIM(IYEAR)//'.log'
    OPEN(LOGNAM, FILE=CLOGTXT, FORM='FORMATTED')

    ! *** 2. read input river map
    CALL INIT_MAP

    ! *** 3. set topography
    CALL INIT_TOPO

    ! *** 4. set initial conditions (zero-storage / restart)
    CALL INIT_COND

    ! *** 5. set time
    CALL INIT_TIME

    ! test array
    arr1 = (/3,  2, 1/)

    call log(LEVEL_INFO, 'Initialized level 1')

    ! ================================================
    WRITE(LOGNAM,*) '*** CONTROLTSTP: ***'

    CALL INIT_DIAG
    CALL SPLITDATE(KYYYYMMDD,IYYYY,IMM,IDD)

#ifdef UseCDF
    KMINP=KMIN
!    CALL CONTROL_INP(IYYYY,IMM,IDD)
#endif

    write(LOGNAM,*) '[main] calculation start'
    write (LOGNAM,*) 'YYYY/MM/DD_HH:MM  ADP_IT  Sto(t)  Sto(flw)  Sto(t+1)  Sto(riv)  Sto(fld)  Inf  Out   Err   Fare'

    TOTAL = ETIME(ELAPSED)

    call log(LEVEL_INFO, 'Initialized level 2')


    istep = 1
  end function initialize


  !> Performs a single timestep with the current model.
  integer(c_int) function update(dt) result(ierr) bind(C,name="update")
    !DEC$ ATTRIBUTES DLLEXPORT::update

    !< Custom timestep size, use -1 to use model default.
    real(c_double), value, intent(in) :: dt
    INTEGER(KIND=JPIM)    ::    IYYYYP,IMMP,IDDP                   !! date (start of time step)

    ierr = 0

    !!========
    !! 1. UPDATE TIME
    IYYYYP=IYYYY                                !! IYYYP, IMMP, IDDP , KYYYMMDDP, KHHMMP, KMINP: (start of the current time step)
    IMMP=IMM
    IDDP=IDD
    KYYYYMMDDP=KYYYYMMDD
    KHHMMP=KHHMM
    KMINP=KMIN
    KSTEP=ISTEP
    KMIN=INT ( KMINSTART+KSTEP*(DT/60_JPIM) )   !! update KMIN (end of the current time step)
    CALL MIN2DATE(KMIN,KYYYYMMDD,KHHMM)         !! IYYY, IMM, IDD , KYYYMMDD, KHHMM, KMIN: (end of the current time step)
    CALL SPLITDATE(KYYYYMMDD,IYYYY,IMM,IDD)

    ! ==========================================
    !  2. UPDATE BOUNDARY SEA LEVEL AT 10 min
    IF (( LBOUNDSL ) .and. MOD(KHHMMP,10)==0 ) THEN
       CALL CONTROL_BOUND(IYYYYP,IMMP,IDDP,KHHMMP)
    ENDIF

    ! ==========================================
    !  3. UPDATE RUNOFF AT 00H
    IF ( KHHMMP==0 ) THEN
       IF ( LBMIROF ) THEN
         CALL CONTROL_INP_BMI
       ELSE
         CALL CONTROL_INP(IYYYYP,IMMP,IDDP)
       ENDIF
    ENDIF

    ! ==========================================
    !  4. MAIN CALCULATION
    CALL CONTROL_PHY

    ! ============================================
    !  5. UPDATE DIAGNOSTICS, AND WRITE TO OUTPUT IF THATS THE CASE
    CALL CONTROL_OUT(IYYYYP,IMMP,IDDP)

    !  ===============================
    !  6. WRTE_RESTART
    CALL CONTROL_REST(ISTEP,IYYYY,IMM,IDD)
    !  ! ============================================
    !  !  7. UPDATE PROGNOSTICS T+1 > T
    !  D2RIVSTO_IN(:,:)=D2RIVSTO_OUT(:,:)
    !  D2FLDSTO_IN(:,:)=D2FLDSTO_OUT(:,:)
    !  D2RIVOUT_IN(:,:)=D2RIVOUT_OUT(:,:)
    !  D2FLDOUT_IN(:,:)=D2FLDOUT_OUT(:,:)
    !  D1PTHFLW_IN(:,:)  =D1PTHFLW_OUT(:,:)

    IF (IDDP .NE. IDD) THEN
          TOTAL1=ETIME(ELAPSED)
          CALL CPU_TIME(TOTAL1)
          TSEC=TOTAL1-TOTAL
          TOTAL=TOTAL1
       WRITE(LOGNAM,'(A6,I5,I3,I3,I5,I10,F7.2,F7.2,A4)') 'TIME: ',IYYYY,IMM,IDD,KHHMM,KSTEP,100.*REAL(KSTEP)/REAL(NSTEPS),TSEC,'sec'
    ENDIF

    ISTEP = ISTEP + 1


  end function update


  ! Void function is a subroutine
  subroutine get_var_type(c_var_name, c_type_name)  bind(C, name="get_var_type")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_var_type

    character(kind=c_char), intent(in) :: c_var_name(*)
    character(kind=c_char), intent(out) :: c_type_name(MAXSTRINGLEN)

    character(len=strlen(c_var_name)) :: var_name
    character(len=MAXSTRINGLEN) :: type_name

    var_name = char_array_to_string(c_var_name)

    select case(var_name)
    case('roffin')
       type_name = 'float'
    case('runoff')
       type_name = 'float'
    case('runoff_m')
       type_name = 'float'
    case('rivinf')
       type_name = 'float'
    case('rivdph')
       type_name = 'float'
    case('rivvel')
       type_name = 'float'
    case('fldinf')
       type_name = 'float'
    case('flddph')
       type_name = 'float'
    case('fldfrc')
       type_name = 'float'
    case('fldare')
       type_name = 'float'
    case('pthout')
       type_name = 'float'
    case('pthinf')
       type_name = 'float'
    case('sfcelv')
       type_name = 'float'
    case('outflw')
       type_name = 'float'
    case('storge')
       type_name = 'float'
    case('rivout_avg')
       type_name = 'float'
    case('outflw_avg')
       type_name = 'float'
    case('fldout_avg')
       type_name = 'float'
    case('rivvel_avg')
      type_name = 'float'
    case('pthout_avg')
       type_name = 'float'
    case('pthflw_avg')
       type_name = 'float'
    case('rivsto_out')
       type_name = 'float'
    case('fldsto_out')
       type_name = 'float'
    case default
       type_name = 'unknown'
    end select

    c_type_name = string_to_char_array(trim(type_name))

  end subroutine get_var_type

  subroutine get_var_rank(c_var_name, rank) bind(C, name="get_var_rank")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_var_rank

    character(kind=c_char), intent(in) :: c_var_name(*)
    integer(c_int), intent(out) :: rank

    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    ! Store the name
    var_name = char_array_to_string(c_var_name)

    select case(var_name)
    case("roffin")
       rank = size(shape(R2ROFFIN))
    case("runoff")
       rank = size(shape(R2RUNOFF))
    case("runoff_m")
       rank = size(shape(R2RUNOFF_M))
    case("rivinf")
       rank = size(shape(R2RIVINF))
    case("rivdph")
       rank = size(shape(R2RIVDPH))
    case("rivvel")
       rank = size(shape(R2RIVVEL))
    case("fldinf")
       rank = size(shape(R2FLDINF))
    case("flddph")
       rank = size(shape(R2FLDDPH))
    case("fldfrc")
       rank = size(shape(R2FLDFRC))
    case("fldare")
       rank = size(shape(R2FLDARE))
    case("pthout")
       rank = size(shape(R2PTHOUT))
    case("pthinf")
       rank = size(shape(R2PTHINF))
    case("sfcelv")
       rank = size(shape(R2SFCELV))
    case("outflw")
       rank = size(shape(R2OUTFLW))
    case("storge")
       rank = size(shape(R2STORGE))
    case("rivout_avg")
       rank = size(shape(R2RIVOUT_AVG))
    case("outflw_avg")
       rank = size(shape(R2OUTFLW_AVG))
    case("fldout_avg")
       rank = size(shape(R2FLDOUT_AVG))
    case("rivvel_avg")
       rank = size(shape(R2RIVVEL_AVG))
    case("pthout_avg")
       rank = size(shape(R2PTHOUT_AVG))
    case("pthflw_avg")
       rank = size(shape(R1PTHFLW_AVG))
    case("rivsto_out")
       rank = size(shape(R2RIVSTO_OUT))
    case("fldsto_out")
       rank = size(shape(R2FLDSTO_OUT))
    case default
       rank = 0
    end select
  end subroutine get_var_rank

  subroutine get_var_shape(c_var_name, the_shape) bind(C, name="get_var_shape")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_var_shape

    character(kind=c_char), intent(in) :: c_var_name(*)
    integer(c_int), intent(inout) :: the_shape(MAXDIMS)

    character(len=strlen(c_var_name)) :: var_name

    var_name = char_array_to_string(c_var_name)
    the_shape = (/0, 0, 0, 0, 0, 0/)

    select case(var_name)
    case("roffin")
       the_shape(1:size(shape(R2ROFFIN)))     = shape(R2ROFFIN)
    case("runoff")
       the_shape(1:size(shape(R2RUNOFF)))     = shape(R2RUNOFF)
    case("runoff_m")
       the_shape(1:size(shape(R2RUNOFF_M)))   = shape(R2RUNOFF_M)
    case("rivinf")
       the_shape(1:size(shape(R2RIVINF)))     = shape(R2RIVINF)
    case("rivdph")
       the_shape(1:size(shape(R2RIVDPH)))     = shape(R2RIVDPH)
    case("rivvel")
       the_shape(1:size(shape(R2RIVVEL)))     = shape(R2RIVVEL)
    case("fldinf")
       the_shape(1:size(shape(R2FLDINF)))     = shape(R2FLDINF)
    case("flddph")
       the_shape(1:size(shape(R2FLDDPH)))     = shape(R2FLDDPH)
    case("fldfrc")
       the_shape(1:size(shape(R2FLDFRC)))     = shape(R2FLDFRC)
    case("fldare")
       the_shape(1:size(shape(R2FLDARE)))     = shape(R2FLDARE)
    case("pthout")
       the_shape(1:size(shape(R2PTHOUT)))     = shape(R2PTHOUT)
    case("pthinf")
       the_shape(1:size(shape(R2PTHINF)))     = shape(R2PTHINF)
    case("sfcelv")
       the_shape(1:size(shape(R2SFCELV)))     = shape(R2SFCELV)
    case("outflw")
       the_shape(1:size(shape(R2OUTFLW)))     = shape(R2OUTFLW)
    case("storge")
       the_shape(1:size(shape(R2STORGE)))     = shape(R2STORGE)
    case("rivout_avg")
       the_shape(1:size(shape(R2RIVOUT_AVG))) = shape(R2RIVOUT_AVG)
    case("outflw_avg")
       the_shape(1:size(shape(R2OUTFLW_AVG))) = shape(R2OUTFLW_AVG)
    case("fldout_avg")
       the_shape(1:size(shape(R2FLDOUT_AVG))) = shape(R2FLDOUT_AVG)
    case("rivvel_avg")
       the_shape(1:size(shape(R2RIVVEL_AVG))) = shape(R2RIVVEL_AVG)
    case("pthout_avg")
       the_shape(1:size(shape(R2PTHOUT_AVG))) = shape(R2PTHOUT_AVG)
    case("pthflw_avg")
       the_shape(1:size(shape(R1PTHFLW_AVG))) = shape(R1PTHFLW_AVG)
    case("rivsto_out")
       the_shape(1:size(shape(R2RIVSTO_OUT))) = shape(R2RIVSTO_OUT)
    case("fldsto_out")
       the_shape(1:size(shape(R2FLDSTO_OUT))) = shape(R2FLDSTO_OUT)
    end select
  end subroutine get_var_shape


  subroutine get_var(c_var_name, x) bind(C, name="get_var")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_var

    ! Return a pointer to the variable

    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), intent(inout) :: x

    character(len=strlen(c_var_name)) :: var_name
    ! Store the name

    var_name = char_array_to_string(c_var_name)

    select case(var_name)
    case("roffin")
       x = c_loc(R2ROFFIN)
    case("runoff")
       x = c_loc(R2RUNOFF)
    case("runoff_m")
       x = c_loc(R2RUNOFF_M)
    case("rivinf")
       x = c_loc(R2RIVINF)
    case("rivdph")
       x = c_loc(R2RIVDPH)
    case("rivvel")
       x = c_loc(R2RIVVEL)
    case("fldinf")
       x = c_loc(R2FLDINF)
    case("flddph")
       x = c_loc(R2FLDDPH)
    case("fldfrc")
       x = c_loc(R2FLDFRC)
    case("fldare")
       x = c_loc(R2FLDARE)
    case("pthout")
       x = c_loc(R2PTHOUT)
    case("pthinf")
       x = c_loc(R2PTHINF)
    case("sfcelv")
       x = c_loc(R2SFCELV)
    case("outflw")
       x = c_loc(R2OUTFLW)
    case("storge")
       x = c_loc(R2STORGE)
    case("rivout_avg")
       x = c_loc(R2RIVOUT_AVG)
    case("outflw_avg")
       x = c_loc(R2OUTFLW_AVG)
    case("fldout_avg")
       x = c_loc(R2FLDOUT_AVG)
    case("rivvel_avg")
       x = c_loc(R2RIVVEL_AVG)
    case("pthout_avg")
       x = c_loc(R2PTHOUT_AVG)
    case("pthflw_avg")
       x = c_loc(R1PTHFLW_AVG)
    case("rivsto_out")
       x = c_loc(R2RIVSTO_OUT)
    case("fldsto_out")
       x = c_loc(R2FLDSTO_OUT)
    case default
       x = c_null_ptr
    end select

  end subroutine get_var

  subroutine set_var(c_var_name, xptr) bind(C, name="set_var")
    !DEC$ ATTRIBUTES DLLEXPORT :: set_var
    ! Return a pointer to the variable
    use iso_c_binding, only: c_double, c_char, c_loc, c_f_pointer

    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), value, intent(in) :: xptr


    real(c_double), pointer :: x_1d_double_ptr(:)
    real(c_double), pointer :: x_2d_double_ptr(:,:)
    real(c_double), pointer :: x_3d_double_ptr(:,:,:)
    integer(c_int), pointer :: x_1d_int_ptr(:)
    integer(c_int), pointer :: x_2d_int_ptr(:,:)
    integer(c_int), pointer :: x_3d_int_ptr(:,:,:)
    real(c_float), pointer  :: x_1d_float_ptr(:)
    real(c_float), pointer  :: x_2d_float_ptr(:,:)
    real(c_float), pointer  :: x_3d_float_ptr(:,:,:)
    logical(c_bool), pointer  :: x_1d_bool_ptr(:)
    logical(c_bool), pointer  :: x_2d_bool_ptr(:,:)
    logical(c_bool), pointer  :: x_3d_bool_ptr(:,:,:)

    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    ! Store the name
    var_name = char_array_to_string(c_var_name)

    select case(var_name)
    case("roffin")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2ROFFIN))
       R2ROFFIN(:, :) = x_2d_float_ptr
    case("runoff")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2RUNOFF))
       R2RUNOFF(:, :) = x_2d_float_ptr
    case("runoff_m")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2RUNOFF_M))
       R2RUNOFF_M(:, :) = x_2d_float_ptr
    case("rivinf")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2RIVINF))
       R2RIVINF(:, :) = x_2d_float_ptr
    case("rivdph")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2RIVDPH))
       R2RIVDPH(:, :) = x_2d_float_ptr
    case("rivvel")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2RIVVEL))
       R2RIVVEL(:, :) = x_2d_float_ptr
    case("fldinf")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2FLDINF))
       R2FLDINF(:, :) = x_2d_float_ptr
    case("flddph")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2FLDDPH))
       R2FLDDPH(:, :) = x_2d_float_ptr
    case("fldfrc")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2FLDFRC))
       R2FLDFRC(:, :) = x_2d_float_ptr
    case("fldare")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2FLDARE))
       R2FLDARE(:, :) = x_2d_float_ptr
    case("pthout")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2PTHOUT))
       R2PTHOUT(:, :) = x_2d_float_ptr
    case("pthinf")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2PTHINF))
       R2PTHINF(:, :) = x_2d_float_ptr
    case("sfcelv")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2SFCELV))
       R2SFCELV(:, :) = x_2d_float_ptr
    case("outflw")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2OUTFLW))
       R2OUTFLW(:, :) = x_2d_float_ptr
    case("storge")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2STORGE))
       R2STORGE(:, :) = x_2d_float_ptr
    case("rivout_avg")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2RIVOUT_AVG))
       R2RIVOUT_AVG(:, :) = x_2d_float_ptr
    case("outflw_avg")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2OUTFLW_AVG))
       R2OUTFLW_AVG(:, :) = x_2d_float_ptr
    case("fldout_avg")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2FLDOUT_AVG))
       R2FLDOUT_AVG(:, :) = x_2d_float_ptr
    case("rivvel_avg")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2RIVVEL_AVG))
       R2RIVVEL_AVG(:, :) = x_2d_float_ptr
    case("pthout_avg")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2PTHOUT_AVG))
       R2PTHOUT_AVG(:, :) = x_2d_float_ptr
    case("pthflw_avg")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R1PTHFLW_AVG))
       R1PTHFLW_AVG(:, :) = x_2d_float_ptr
    case("rivsto_out")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2RIVSTO_OUT))
       R2RIVSTO_OUT(:, :) = x_2d_float_ptr
    case("fldsto_out")
       call c_f_pointer(xptr, x_2d_float_ptr, shape(R2FLDSTO_OUT))
       R2FLDSTO_OUT(:, :) = x_2d_float_ptr
    case default
       ! do nothing
    end select

  end subroutine set_var

  subroutine set_var_slice(c_var_name, start, count, xptr) bind(C, name="set_var_slice")

    character(kind=c_char), intent(in) :: c_var_name(*)
    integer(c_int), intent(inout) :: start(MAXDIMS)
    integer(c_int), intent(inout) :: count(MAXDIMS)
    type(c_ptr), value, intent(in) :: xptr

    real(c_double), pointer :: x_1d_double_ptr(:)
    real(c_double), pointer :: x_2d_double_ptr(:,:)
    real(c_double), pointer :: x_3d_double_ptr(:,:,:)
    integer(c_int), pointer :: x_1d_int_ptr(:)
    integer(c_int), pointer :: x_2d_int_ptr(:,:)
    integer(c_int), pointer :: x_3d_int_ptr(:,:,:)
    real(c_float), pointer  :: x_1d_float_ptr(:)
    real(c_float), pointer  :: x_2d_float_ptr(:,:)
    real(c_float), pointer  :: x_3d_float_ptr(:,:,:)
    logical(c_bool), pointer  :: x_1d_bool_ptr(:)
    logical(c_bool), pointer  :: x_2d_bool_ptr(:,:)
    logical(c_bool), pointer  :: x_3d_bool_ptr(:,:,:)

    integer(c_int) :: start1(MAXDIMS)
    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    var_name = char_array_to_string(c_var_name)

    start1 = start + 1
    select case(var_name)
    case("sfcelv")
       call c_f_pointer(xptr, x_1d_double_ptr, (/count(1)/))
       arr1(start1(1):(start1(1)+count(1)-1)) = x_1d_double_ptr
    end select


  end subroutine set_var_slice

  subroutine get_current_time(time) bind(C, name="get_current_time")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_current_time
    use mod_time

    real(c_double) :: time
    time = kminstart+ kstep * (dt/60)

  end subroutine get_current_time

  subroutine get_start_time(time) bind(C, name="get_start_time")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_start_time

    use mod_time

    real(c_double) :: time
    time = kminstart

  end subroutine get_start_time

  subroutine get_end_time(time) bind(C, name="get_end_time")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_end_time

    use mod_time

    real(c_double) :: time
    time = kminstart + nsteps * (dt/60)

  end subroutine get_end_time

  subroutine get_time_step(dt_out) bind(C, name="get_time_step")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_time_step

    use mod_time

    real(c_double), intent(out) :: dt_out

    dt_out = dt

  end subroutine get_time_step

  SUBROUTINE INIT_DIAG

    ALLOCATE(D2RIVINF(NSEQMAX,1))
    ALLOCATE(D2RIVDPH(NSEQMAX,1))
    ALLOCATE(D2RIVVEL(NSEQMAX,1))

    ALLOCATE(D2FLDINF(NSEQMAX,1))
    ALLOCATE(D2FLDDPH(NSEQMAX,1))
    ALLOCATE(D2FLDFRC(NSEQMAX,1))
    ALLOCATE(D2FLDARE(NSEQMAX,1))

    ALLOCATE(D2PTHOUT(NSEQMAX,1))
    ALLOCATE(D2PTHINF(NSEQMAX,1))

    ALLOCATE(D2SFCELV(NSEQMAX,1))
    ALLOCATE(D2OUTFLW(NSEQMAX,1))
    ALLOCATE(D2STORGE(NSEQMAX,1))

    ALLOCATE(D2RIVOUT_AVG(NSEQMAX,1))
    ALLOCATE(D2FLDOUT_AVG(NSEQMAX,1))
    ALLOCATE(D2RIVVEL_AVG(NSEQMAX,1))
    ALLOCATE(D2PTHOUT_AVG(NSEQMAX,1))
    ALLOCATE(D2OUTFLW_AVG(NSEQMAX,1))    !! bugfix v362
    ALLOCATE(D1PTHFLW_AVG(NPTHOUT,NPTHLEV))
    ALLOCATE(D2OUTINS(NSEQMAX,1))

    D2RIVINF(:,:)=0.D0
    D2RIVDPH(:,:)=0.D0
    D2RIVVEL(:,:)=0.D0

    D2FLDINF(:,:)=0.D0
    D2FLDDPH(:,:)=0.D0
    D2FLDFRC(:,:)=0.D0
    D2FLDARE(:,:)=0.D0

    D2PTHOUT(:,:)=0.D0
    D2PTHINF(:,:)=0.D0

    D2SFCELV(:,:)=0.D0
    D2OUTFLW(:,:)=0.D0
    D2STORGE(:,:)=0.D0

    D2RIVOUT_AVG(:,:)=0.D0
    D2FLDOUT_AVG(:,:)=0.D0
    D2RIVVEL_AVG(:,:)=0.D0
    D2PTHOUT_AVG(:,:)=0.D0
    D1PTHFLW_AVG(:,:)=0.D0
    D2OUTFLW_AVG(:,:)=0.D0     !! bugfix v362
    NADD=0

    D2OUTINS(:,:)=0.D0

    ALLOCATE(R2RIVINF(NX,NY))
    ALLOCATE(R2RIVDPH(NX,NY))
    ALLOCATE(R2RIVVEL(NX,NY))
    ALLOCATE(R2FLDINF(NX,NY))
    ALLOCATE(R2FLDDPH(NX,NY))
    ALLOCATE(R2FLDFRC(NX,NY))
    ALLOCATE(R2FLDARE(NX,NY))
    ALLOCATE(R2PTHOUT(NX,NY))
    ALLOCATE(R2PTHINF(NX,NY))
    ALLOCATE(R2SFCELV(NX,NY))
    ALLOCATE(R2OUTFLW(NX,NY))
    ALLOCATE(R2STORGE(NX,NY))
    ALLOCATE(R2RIVOUT_AVG(NX,NY))
    ALLOCATE(R2OUTFLW_AVG(NX,NY))
    ALLOCATE(R2FLDOUT_AVG(NX,NY))
    ALLOCATE(R2RIVVEL_AVG(NX,NY))
    ALLOCATE(R2PTHOUT_AVG(NX,NY))
    ALLOCATE(R1PTHFLW_AVG(NX,NY))
    ALLOCATE(R2RIVSTO_OUT(NX,NY))
    ALLOCATE(R2FLDSTO_OUT(NX,NY))
  
    R2RIVINF(:,:)    =0.D0
    R2RIVDPH(:,:)    =0.D0
    R2RIVVEL(:,:)    =0.D0
    R2FLDINF(:,:)    =0.D0
    R2FLDDPH(:,:)    =0.D0
    R2FLDFRC(:,:)    =0.D0
    R2FLDARE(:,:)    =0.D0
    R2PTHOUT(:,:)    =0.D0
    R2PTHINF(:,:)    =0.D0
    R2SFCELV(:,:)    =0.D0
    R2OUTFLW(:,:)    =0.D0
    R2STORGE(:,:)    =0.D0
    R2RIVOUT_AVG(:,:)=0.D0
    R2OUTFLW_AVG(:,:)=0.D0
    R2FLDOUT_AVG(:,:)=0.D0
    R2RIVVEL_AVG(:,:)=0.D0
    R2PTHOUT_AVG(:,:)=0.D0
    R1PTHFLW_AVG(:,:)=0.D0
    R2RIVSTO_OUT(:,:)=0.D0
    R2FLDSTO_OUT(:,:)=0.D0

  END SUBROUTINE INIT_DIAG
end module demo_model
