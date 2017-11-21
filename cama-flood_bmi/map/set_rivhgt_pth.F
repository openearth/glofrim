      program calc_rivhgt_pth
! ================================================
      implicit none
! ===================
! calculation type
      character*128       ::  type            !! 'bin' for binary, 'cdf' for netCDF
! river network map parameters
      integer             ::  ix, iy, jx, jy
      integer             ::  nx, ny          !! river map grid number
! river netwrok map
      integer,allocatable ::  nextx(:,:)      !! downstream x
      integer,allocatable ::  nexty(:,:)      !! downstream y
      real,allocatable    ::  elevtn(:,:)     !! elevation                 [m]
      real,allocatable    ::  lon(:), lat(:)  !! longitude, latitude       [deg]
! variable
      real,allocatable    ::  rivhgt(:,:)     !! channel depth (empirical) [m]

! file
      character*128       ::  diminfo
      parameter              (diminfo='./diminfo_1deg.txt')
!
      character*128       ::  cnextxy, celevtn, crivhgt, crivhgt_pth
      parameter              (cnextxy='./nextxy.bin')
      parameter              (celevtn='./elevtn.bin')
      parameter              (crivhgt='./rivhgt.bin')
      parameter              (crivhgt_pth='./rivhgt_pth.bin')
!
      character*128       ::  cfldpth
      parameter              (cfldpth='./fldpth.txt')

      integer             ::  ipath, npath
      integer             ::  pth_ix, pth_iy, pth_jx, pth_jy
      real                ::  pth_elv

! Undefined Values
      integer             ::  imis                !! integer undefined value
      real                ::  rmis                !! real    undefined value
      parameter              (imis = -9999)
      parameter              (rmis = 1.e+20)

      character*128       ::  buf
! ================================================
      call getarg(1,type)

      if( type=='cdf')then
        print *, 'calculation for netCDF map'
      else
        type='bin'
        print *, 'calculation for binary map'
      endif

      open(11,file=diminfo,form='formatted')
      read(11,*    ) nx
      read(11,*    ) ny
print *, nx, ny

      allocate(nextx(nx,ny),nexty(nx,ny),elevtn(nx,ny),rivhgt(nx,ny))
      allocate(lon(nx),lat(ny))

print *, cnextxy
      open(11,file=cnextxy,form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) nextx
      read(11,rec=2) nexty
      close(11)

print *, celevtn
      open(11,file=celevtn,form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) elevtn
      close(11)

print *, crivhgt
      open(11,file=crivhgt,form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) rivhgt
      close(11)

! ===============================

      open(21,file=cfldpth,form='formatted')
      read(21,*) npath

      do ipath=1, npath
        read(21,*) pth_ix, pth_iy, pth_jx, pth_jy, buf, pth_elv
        if( pth_elv/=9999. )then
          ix=pth_ix
          iy=pth_iy
          do while( elevtn(ix,iy)-rivhgt(ix,iy)>pth_elv )
            rivhgt(ix,iy)=elevtn(ix,iy)-pth_elv
            if( nextx(ix,iy)<0 ) exit
            jx=nextx(ix,iy)
            jy=nexty(ix,iy)
            ix=jx
            iy=jy
          end do

          ix=pth_jx
          iy=pth_jy
          do while( elevtn(ix,iy)-rivhgt(ix,iy)>pth_elv )
            rivhgt(ix,iy)=elevtn(ix,iy)-pth_elv
            if( nextx(ix,iy)<0 ) exit
            jx=nextx(ix,iy)
            jy=nexty(ix,iy)
            ix=jx
            iy=jy
          end do
        endif
      end do

! ==============================

print *, crivhgt_pth
      open(11,file=crivhgt_pth,form='unformatted',access='direct',recl=4*nx*ny)
      write(11,rec=1) rivhgt
      close(11)

      end program calc_rivhgt_pth


