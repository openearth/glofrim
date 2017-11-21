      program COMBINE_SHEDS
! ===============================================
      implicit none
      character*128          ::  region_info                   !! regionalization info
      parameter                 (region_info='./region_info.txt')
      character*128          ::  global_dir                    !! global map directory

      character*128          ::  global_param                  !! global map parameter file
      real                   ::  lon_ori, lat_ori              !! west and north edge of global map
      integer                ::  narea                         !! global map number of sheds map area
      real                   ::  csize                         !! sheds map pixel size [deg]

      character*128          ::  region_param, region_list_loc
      parameter                 (region_param='../params.txt')
      parameter                 (region_list_loc='../sheds/location.txt')
      integer                ::  nXX, nYY                      !! river map grid numbers
      integer                ::  dXX, dYY                      !! river map domain shift from global map
      real                   ::  gsize                         !! river map grid size [deg]
      real                   ::  west, east, north, south      !! domain

      integer                ::  ix, iy, dx, dy                        
      integer                ::  nx, ny                        !! pixel number of regional sheds map
!
      character*64           ::  list_loc
      character*3            ::  area                          !! area code
      integer                ::  i, iarea                      !! area ID
      character*3,allocatable::  list_area(:)                 !! area code
      integer,allocatable    ::  list_nx(:),  list_ny(:)
      real,allocatable       ::  list_lon(:), list_lat(:)

      integer                ::  ix0, iy0, nx0, ny0            !! input sheds map dimention
      real                   ::  west0, north0, east0, south0

! regional sheds map
      integer*2,allocatable  ::  nextx(:,:), nexty(:,:)
      real,allocatable       ::  elevtn(:,:)
      real,allocatable       ::  uparea(:,:),rivwth(:,:)

! input sheds map
      integer*2,allocatable  ::  nextx0(:,:), nexty0(:,:)
      real,allocatable       ::  elevtn0(:,:)
      real,allocatable       ::  uparea0(:,:),rivwth0(:,:)
      real,allocatable       ::  lon0(:), lat0(:)

! files
      character*64           ::  rfile1, rfile2, rfile3, rfile4
      character*64           ::  wfile1, wfile2, wfile3, wfile4
      character*64           ::  sheds
      integer                ::  ios
      character*64           ::  buf
! ===============================================
      open(11,file=region_info,form='formatted')
      read(11,*) global_dir
      close(11)

      sheds=trim(global_dir)//'sheds/'
      list_loc=trim(sheds)//'location.txt'

      global_param=trim(global_dir)//'params.txt'
      open(11,file=global_param,form='formatted')
      read(11,*) lon_ori
      read(11,*) lat_ori
      read(11,*) 
      read(11,*) 
      read(11,*) gsize
      read(11,*) 
      read(11,*) narea
      read(11,*) csize
      close(11)

      open(11,file=region_param,form='formatted')
      read(11,*) west
      read(11,*) north
      read(11,*) nXX
      read(11,*) nYY
      read(11,*) 
      read(11,*) 
      read(11,*) 
      read(11,*) csize
      close(11)

      east =west +real(nXX)*gsize
      south=north-real(nYY)*gsize

      nx=int( (east-west)  /csize )
      ny=int( (north-south)/csize )

      dXX=int( (west-lon_ori)  /gsize )
      dYY=int( (lat_ori-north) /gsize )

      print '(a8,4f8.2)', 'W E N S ', west, east, north, south
      print '(a16,4i8)' , 'nx ny dXX dYY ', nx, ny, dXX, dYY

      allocate(nextx(nx,ny),nexty(nx,ny),elevtn(nx,ny),uparea(nx,ny),rivwth(nx,ny))
      nextx(:,:)=-9999
      nexty(:,:)=-9999
      elevtn(:,:)=-9999
      uparea(:,:)=-9999
      rivwth(:,:)=-9999

      allocate(list_area(narea),list_lon(narea),list_lat(narea),list_nx(narea),list_ny(narea))

      open(11,file=list_loc,form='formatted')
      read(11,*)
      read(11,*) buf, (list_area(i) ,i=1,narea)
      read(11,*) buf, (list_lon(i)  ,i=1,narea)
      read(11,*) buf, (list_lat(i)  ,i=1,narea)
      read(11,*) buf, (list_nx(i)   ,i=1,narea)
      read(11,*) buf, (list_ny(i)   ,i=1,narea)
      close(11)

      do i=1, narea
        area=list_area(i)
        iarea=i
        nx0=list_nx(i)
        ny0=list_ny(i)
        west0 =list_lon(i)
        north0=list_lat(i)

        east0 =west0 +real(nx0)*csize
        south0=north0-real(ny0)*csize

        if( west0>east .or. east0<west .or. north0<south .or. south0>north )then
          print *, '      out of domain: ', trim(area)
          goto 2000
        endif

        allocate(nextx0(nx0,ny0),nexty0(nx0,ny0),elevtn0(nx0,ny0),uparea0(nx0,ny0),rivwth0(nx0,ny0))
        allocate(lon0(nx0),lat0(ny0))

        rfile1=trim(sheds)//trim(area)//'.nextxy.bsq'
        rfile2=trim(sheds)//trim(area)//'.elevtn.flt'
        rfile3=trim(sheds)//trim(area)//'.uparea.flt'
        rfile4=trim(sheds)//trim(area)//'.rivwth.flt'

        print *, rfile1
        open(21,file=rfile1,form='unformatted',access='direct',recl=2*nx0*ny0,status='old',iostat=ios)
        if( ios==0 )then
          read(21,rec=1) nextx0
          read(21,rec=2) nexty0
          close(21)
        else
          print *, '*******************'
          print *, 'no data: ', rfile1
          goto 1000
        endif

        print *, rfile2
        open(21,file=rfile2,form='unformatted',access='direct',recl=4*nx0*ny0,status='old',iostat=ios)
        if( ios==0 )then
          read(21,rec=1) elevtn0
          close(21)
        else
          print *, '*******************'
          print *, 'no data: ', rfile2
          goto 1000
        endif

        print *, rfile3
        open(21,file=rfile3,form='unformatted',access='direct',recl=4*nx0*ny0,status='old',iostat=ios)
        if( ios==0 )then
          read(21,rec=1) uparea0
          close(21)
        else
          print *, '*******************'
          print *, 'no data: ', rfile3
          goto 1000
        endif

        print *, rfile4
        open(21,file=rfile4,form='unformatted',access='direct',recl=4*nx0*ny0,status='old',iostat=ios)
        if( ios==0 )then
          read(21,rec=1) rivwth0
          close(21)
        else
          print *, '*******************'
          print *, 'no data: ', rfile4
          goto 1000
        endif

        do ix0=1, nx0
          lon0(ix0)=west0+(real(ix0)-0.5)*csize
        end do

        do iy0=1, ny0
          lat0(iy0)=north0-(real(iy0)-0.5)*csize
        end do

        do iy0=1, ny0
          do ix0=1, nx0

            if( lon0(ix0)>west .and. lon0(ix0)<east .and. lat0(iy0)>south .and. lat0(iy0)<north )then

              ix=int( (lon0(ix0) -west)/csize )+1
              iy=int( (north-lat0(iy0))/csize )+1

              if( nextx(ix,iy)==-9999 )then
                if( nextx0(ix0,iy0)>0 )then
                  dx=nextx0(ix0,iy0)-ix0
                  dy=nexty0(ix0,iy0)-iy0
                  nextx(ix,iy)=ix+dx
                  nexty(ix,iy)=iy+dy
                  if( nextx(ix,iy)<1 .or. nextx(ix,iy)>nx .or. nexty(ix,iy)<1 .or. nexty(ix,iy)>ny )then
                    nextx(ix,iy)=-9
                    nexty(ix,iy)=-9
                  endif
                else
                  nextx(ix,iy)=nextx0(ix0,iy0)
                  nexty(ix,iy)=nexty0(ix0,iy0)
                endif
              endif

              if( nextx0(ix0,iy0)/=-9999 )then
                elevtn(ix,iy)=elevtn0(ix0,iy0)
                uparea(ix,iy)=uparea0(ix0,iy0)
                rivwth(ix,iy)=rivwth0(ix0,iy0)
              endif

            endif
          end do
        end do
 1000   continue
        deallocate(nextx0,nexty0,elevtn0,uparea0,rivwth0,lon0,lat0)
 2000   continue
      end do

      wfile1='../sheds/reg.nextxy.bsq'
      open(21,file=wfile1,form='unformatted',access='direct',recl=2*nx*ny)
      write(21,rec=1) nextx
      write(21,rec=2) nexty
      close(21)

      wfile2='../sheds/reg.elevtn.flt'
      open(21,file=wfile2,form='unformatted',access='direct',recl=4*nx*ny)
      write(21,rec=1) elevtn
      close(21)

      wfile3='../sheds/reg.uparea.flt'
      open(21,file=wfile3,form='unformatted',access='direct',recl=4*nx*ny)
      write(21,rec=1) uparea
      close(21)

      wfile4='../sheds/reg.rivwth.flt'
      open(21,file=wfile4,form='unformatted',access='direct',recl=4*nx*ny)
      write(21,rec=1) rivwth
      close(21)

      open(11,file=region_list_loc,form='formatted')
      write(11,'(a9)' )      'code    1'
      write(11,'(a11)')      'area    reg'
      write(11,'(a8,f12.3)') 'west    ', west
      write(11,'(a8,f12.3)') 'north   ', north
      write(11,'(a8,i12)')   'nx      ', nx
      write(11,'(a8,i12)')   'ny      ', ny
      write(11,'(a8,f12.8)') 'csize   ', csize
      close(11)

! ====================

      end program COMBINE_SHEDS

