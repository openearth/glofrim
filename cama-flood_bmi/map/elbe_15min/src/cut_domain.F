      program cut_domain
! ==================================================
      implicit none
!
      character*128          ::  region_info
      parameter                 (region_info='./region_info.txt')
      character*128          ::  global_dir
      real                   ::  west, east, north, south
      real                   ::  west2, east2, north2, south2

      character*128          ::  global_param
      real                   ::  lon_ori, lat_ori
      real                   ::  glon, glat, d1, d2

      character*128          ::  region_dir, region_param
      parameter                 (region_dir='../')
      parameter                 (region_param='../params.txt')
!
      integer                ::  ix, iy, jx, jy, kx, ky
      integer                ::  nx, ny, mx, my, dx, dy
      real                   ::  gsize, csize
      integer                ::  iflp, nflp                  !! number of floodplain layer
!
      character*128          ::  fnextxy0, felevtn0, ffldhgt0, fgrarea0, fuparea0, flonlat0, fnxtdst0, frivlen0, fwidth0
      character*128          ::  fnextxy,  felevtn,  ffldhgt,  fgrarea,  fuparea,  flonlat,  fnxtdst,  frivlen,  fwidth
      character*128          ::  flsmask

      integer,allocatable    ::  nextx0(:,:),  nexty0(:,:)      !!  global maps
      real,allocatable       ::  elevtn0(:,:), fldhgt0(:,:,:)
      real,allocatable       ::  grarea0(:,:), uparea0(:,:)
      real,allocatable       ::  lon0(:,:)   , lat0(:,:)
      real,allocatable       ::  nxtdst0(:,:), rivlen0(:,:)
      real,allocatable       ::  width0(:,:)

      integer,allocatable    ::  nextx(:,:),  nexty(:,:)        !!  regional maps
      real,allocatable       ::  elevtn(:,:), fldhgt(:,:,:)
      real,allocatable       ::  grarea(:,:), uparea(:,:)
      real,allocatable       ::  lon(:,:),    lat(:,:)
      real,allocatable       ::  nxtdst(:,:), rivlen(:,:)
      real,allocatable       ::  width(:,:)
      integer,allocatable    ::  lsmask(:,:)

      real,allocatable       ::  tmp0(:,:), tmp(:,:)            !! for fldhgt I/O

! ==================================================
      open(11,file=region_info,form='formatted')
      read(11,*) global_dir
      read(11,*) west
      read(11,*) east
      read(11,*) north
      read(11,*) south
      close(11)

      global_param=trim(global_dir)//'params.txt'
      open(11,file=global_param,form='formatted')
      read(11,*) lon_ori
      read(11,*) lat_ori
      read(11,*) nx
      read(11,*) ny
      read(11,*) gsize
      read(11,*) nflp
      read(11,*) 
      read(11,*) csize
      close(11)

      allocate(nextx0(nx,ny),nexty0(nx,ny),elevtn0(nx,ny),fldhgt0(nx,ny,nflp),grarea0(nx,ny),uparea0(nx,ny))
      allocate(lon0(nx,ny),  lat0(nx,ny),  nxtdst0(nx,ny),rivlen0(nx,ny), width0(nx,ny))

      print *, 'input:  ', west, east, north, south
      west2=west
      east2=east
      north2=north
      south2=south

      d1=1.e20
      d2=1.e20
      do ix=1, nx
        glon=lon_ori+real(ix-1)*gsize
        if( west2>=glon .and. west2< glon+gsize .and. abs(west2-glon)<d1 )then
          west=glon
          d1=abs(west2-glon)
        endif
        if( east2> glon .and. east2<=glon+gsize .and. abs(east2-glon-gsize)<d2 )then
          east=glon+gsize
          d2=abs(east2-glon-gsize)
        endif
      end do

      d1=1.e20
      d2=1.e20
      do iy=1, ny
        glat=lat_ori-real(iy-1)*gsize
        if( north2> glat-gsize .and. north2<=glat .and. abs(north2-glat)<d1 )then
          north=glat
          d1=abs(north2-glat)
        endif
        if( south2>=glat-gsize .and. south2< glat .and. abs(south2-glat+gsize)<d2 )then
          south=glat-gsize
          d2=abs(south2-glat+gsize)
        endif
      end do
      if( abs(west- west2 )<0.001 ) west =west2
      if( abs(east- east2 )<0.001 ) east =east2
      if( abs(north-north2)<0.001 ) north=north2
      if( abs(south-south2)<0.001 ) south=south2

      print *, 'output: ', west, east, north, south

      mx=int( (east-west)    /gsize )
      my=int( (north-south)  /gsize )
      dx=int( (west-lon_ori) /gsize )
      dy=int( (lat_ori-north)/gsize )

      open(11,file=region_param,form='formatted')
      write(11,'(f12.3,a35)') west,  'west  edge [deg]'
      write(11,'(f12.3,a35)') north, 'north edge [deg]'
      write(11,'(i12  ,a35)') mx,    'grid number (east-west)'
      write(11,'(i12  ,a35)') my,    'grid number (north-south)'
      write(11,'(f12.8,a35)') gsize, 'grid size  [deg]'
      write(11,'(i12  ,a35)') nflp,  'floodplain layers'
      write(11,'(i12  ,a35)') 1,      'area number of hires data'
      write(11,'(f12.8,a35)') csize, 'pixel size  of hires data [deg]'
      close(11)

      write(6,'(f12.3,a35)') west,  'west  edge [deg]'
      write(6,'(f12.3,a35)') north, 'north edge [deg]'
      write(6,'(i12  ,a35)') mx,    'grid number (east-west)'
      write(6,'(i12  ,a35)') my,    'grid number (north-south)'
      write(6,'(f12.8,a35)') gsize, 'grid size  [deg]'
      write(6,'(i12  ,a35)') nflp,  'floodplain layers'
      write(6,'(i12  ,a35)') 1,      'area number of hires data'
      write(6,'(f12.8,a35)') csize, 'pixel size  of hires data [deg]'

      allocate(nextx(mx,my),nexty(mx,my),elevtn(mx,my),fldhgt(mx,my,nflp),grarea(mx,my),uparea(mx,my))
      allocate(lon(mx,my),  lat(mx,my),  nxtdst(mx,my),rivlen(mx,my), width(mx,my))
      allocate(lsmask(mx,my))

      allocate(tmp0(nx,ny), tmp(mx,my))

      fnextxy0=trim(global_dir)//'nextxy.bin'
      felevtn0=trim(global_dir)//'elevtn.bin'
      ffldhgt0=trim(global_dir)//'fldhgt.bin'
      fgrarea0=trim(global_dir)//'grarea.bin'
      fuparea0=trim(global_dir)//'uparea.bin'
      flonlat0=trim(global_dir)//'lonlat.bin'
      fnxtdst0=trim(global_dir)//'nxtdst.bin'
      frivlen0=trim(global_dir)//'rivlen.bin'
      fwidth0=trim(global_dir)//'width.bin'

      fnextxy=trim(region_dir)//'nextxy.bin'
      felevtn=trim(region_dir)//'elevtn.bin'
      ffldhgt=trim(region_dir)//'fldhgt.bin'
      fgrarea=trim(region_dir)//'grarea.bin'
      fuparea=trim(region_dir)//'uparea.bin'
      flonlat=trim(region_dir)//'lonlat.bin'
      fnxtdst=trim(region_dir)//'nxtdst.bin'
      frivlen=trim(region_dir)//'rivlen.bin'
      fwidth=trim(region_dir)//'width.bin'
      flsmask=trim(region_dir)//'lsmask.bin'

print *, 'read global maps'
      print *, fnextxy0
      open(11,file=fnextxy0,form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) nextx0
      read(11,rec=2) nexty0
      close(11)

      print *, felevtn0
      open(11,file=felevtn0,form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) elevtn0
      close(11)

      print *, ffldhgt0
      open(11,file=ffldhgt0,form='unformatted',access='direct',recl=4*nx*ny)
      do iflp=1, nflp
        read(11,rec=iflp) tmp0
        fldhgt0(:,:,iflp)=tmp0(:,:)
      end do
      close(11)

      print *, fgrarea0
      open(11,file=fgrarea0,form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) grarea0
      close(11)

      print *, fuparea0
      open(11,file=fuparea0,form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) uparea0
      close(11)

      print *, flonlat0
      open(11,file=flonlat0,form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) lon0
      read(11,rec=2) lat0
      close(11)

      print *, fnxtdst0
      open(11,file=fnxtdst0,form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) nxtdst0
      close(11)

      print *, frivlen0
      open(11,file=frivlen0,form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) rivlen0
      close(11)

      print *, fwidth0
      open(11,file=fwidth0,form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) width0
      close(11)

print *, 'cut domain'

      nextx(:,:)=-9999
      nexty(:,:)=-9999
      elevtn(:,:)=-9999.
      fldhgt(:,:,:)=-9999.
      grarea(:,:)=-9999
      uparea(:,:)=-9999
      lon(:,:)=-9999
      lat(:,:)=-9999
      nxtdst(:,:)=-9999
      rivlen(:,:)=-9999
      width(:,:)=-9999
      lsmask(:,:)=0

      do iy=1, my
        do ix=1, mx
          jx=ix+dx
          jy=iy+dy

          if( nextx0(jx,jy)>0 )then
            nextx(ix,iy)=nextx0(jx,jy)-dx
            nexty(ix,iy)=nexty0(jx,jy)-dy
            kx=nextx(ix,iy)
            ky=nexty(ix,iy)
            if( kx<1 .or. kx>mx .or. ky<1 .or. ky>my )then  !! if downstream is outside the domain
              nextx(ix,iy)=-9
              nexty(ix,iy)=-9
            endif
          elseif( nextx0(jx,jy)/=-9999 )then     !! if river mouth
            nextx(ix,iy)=nextx0(jx,jy)
            nexty(ix,iy)=nexty0(jx,jy)
          endif

          if( nextx(ix,iy)/=-9999 )then
            elevtn(ix,iy)  =elevtn0(jx,jy)
            fldhgt(ix,iy,:)=fldhgt0(jx,jy,:)
            grarea(ix,iy)  =grarea0(jx,jy)
            uparea(ix,iy)  =uparea0(jx,jy)
            lon(ix,iy)     =lon0(jx,jy)
            lat(ix,iy)     =lat0(jx,jy)
            nxtdst(ix,iy)  =nxtdst0(jx,jy)
            rivlen(ix,iy)  =rivlen0(jx,jy)
            width(ix,iy)  =width0(jx,jy)
            lsmask(ix,iy) =1
          endif
        end do
      end do

print *, 'write reagional maps'
      print *, fnextxy
      open(21,file=fnextxy,form='unformatted',access='direct',recl=4*mx*my)
      write(21,rec=1) nextx
      write(21,rec=2) nexty
      close(21)

      print *, felevtn
      open(21,file=felevtn,form='unformatted',access='direct',recl=4*mx*my)
      write(21,rec=1) elevtn
      close(21)

      print *, ffldhgt
      open(21,file=ffldhgt,form='unformatted',access='direct',recl=4*mx*my)
      do iflp=1, nflp
        tmp(:,:)=fldhgt(:,:,iflp)
        write(21,rec=iflp) tmp
      end do
      close(21)

      print *, fgrarea
      open(21,file=fgrarea,form='unformatted',access='direct',recl=4*mx*my)
      write(21,rec=1) grarea
      close(21)

      print *, fuparea
      open(21,file=fuparea,form='unformatted',access='direct',recl=4*mx*my)
      write(21,rec=1) uparea
      close(21)

      print *, flonlat
      open(21,file=flonlat,form='unformatted',access='direct',recl=4*mx*my)
      write(21,rec=1) lon
      write(21,rec=2) lat
      close(21)

      print *, fnxtdst
      open(21,file=fnxtdst,form='unformatted',access='direct',recl=4*mx*my)
      write(21,rec=1) nxtdst
      close(21)

      print *, frivlen
      open(21,file=frivlen,form='unformatted',access='direct',recl=4*mx*my)
      write(21,rec=1) rivlen
      close(21)

      print *, fwidth
      open(21,file=fwidth,form='unformatted',access='direct',recl=4*mx*my)
      write(21,rec=1) width
      close(21)

      print *, flsmask
      open(21,file=flsmask,form='unformatted',access='direct',recl=4*mx*my)
      write(21,rec=1) lsmask
      close(21)

! ==================================================
      end program cut_domain