      program calc_rivwth
! ================================================
      implicit none
! ===================
! river network map parameters
      integer             ::  ix, iy
      integer             ::  nx, ny          !! river map grid number
! river netwrok map
      integer,allocatable ::  nextx(:,:)      !! downstream x
      integer,allocatable ::  nexty(:,:)      !! downstream y
      real,allocatable    ::  lon(:), lat(:)  !! longitude, latitude [deg]

! variable
      real,allocatable    ::  rivwth(:,:)     !! channel width (empirical) [m]
      real,allocatable    ::  gwdlr(:,:)      !! channel width (GWD-LR)    [m]
! file
      character*128       ::  diminfo
!
      character*128       ::  cnextxy, crivwth, cwidth, cgwdlr
      parameter              (cnextxy='./nextxy.bin')
      parameter              (crivwth='./rivwth.bin')
      parameter              (cwidth='./width.bin')
      parameter              (cgwdlr='./rivwth_gwdlr.bin')
! Undefined Values
      integer             ::  imis                !! integer undefined value
      real                ::  rmis                !! real    undefined value
      parameter              (imis = -9999)
      parameter              (rmis = 1.e+20)
! ================================================
      call getarg(1,diminfo)

      open(11,file=diminfo,form='formatted')
      read(11,*    ) nx
      read(11,*    ) ny
print *, nx, ny

      allocate(nextx(nx,ny),nexty(nx,ny))
      allocate(rivwth(nx,ny),gwdlr(nx,ny))
      allocate(lon(nx),lat(ny))

print *, cnextxy
      open(11,file=cnextxy,form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) nextx
      read(11,rec=2) nexty
      close(11)

print *, crivwth
      open(13,file=crivwth,form='unformatted',access='direct',recl=4*nx*ny)
      read(13,rec=1) rivwth
      close(13)

print *, cwidth
      open(13,file=cwidth,form='unformatted',access='direct',recl=4*nx*ny)
      read(13,rec=1) gwdlr
      close(13)


! ********************************************************************
! please modify this part to calibrate river width and height parameters
! ********************************************************************
      do iy=1, ny
        do ix=1, nx
          if( nextx(ix,iy)/=imis )then
            if( gwdlr(ix,iy)<300 )then
              gwdlr(ix,iy)=max(gwdlr(ix,iy),rivwth(ix,iy))
            endif
            if( gwdlr(ix,iy)>10000. )then
              gwdlr(ix,iy)=10000.
            endif
          else
            gwdlr(ix,iy)=-9999
          endif
        end do
      end do
! ********************************************************************
print *, cgwdlr
      open(21,file=cgwdlr,form='unformatted',access='direct',recl=4*nx*ny)
      write(21,rec=1) gwdlr
      close(21)


      end program calc_rivwth


