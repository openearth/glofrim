      program COMBINE_HIRES
! ===============================================
      implicit none
      character*128          ::  region_info                   !! regionalization info
      parameter                 (region_info='./region_info.txt')
      character*128          ::  global_dir                    !! global map directory

      character*128          ::  global_param                  !! global map parameter file
      real                   ::  lon_ori, lat_ori              !! west and north edge of global map
      integer                ::  narea                         !! global map number of hires map area
      real                   ::  csize                         !! hires map pixel size [deg]

      character*128          ::  region_param, region_list_loc
      parameter                 (region_param='../params.txt')
      parameter                 (region_list_loc='../hires/location.txt')
      integer                ::  nXX, nYY                      !! river map grid numbers
      integer                ::  dXX, dYY                      !! river map domain shift from global map
      real                   ::  gsize                         !! river map grid size [deg]
      real                   ::  west, east, north, south      !! domain

      integer                ::  ix, iy                        
      integer                ::  nx, ny                        !! pixel number of regional hires map
!
      character*64           ::  list_loc
      character*3            ::  area                          !! area code
      integer                ::  i, iarea                      !! area ID
      character*3,allocatable::  list_area(:)                 !! area code
      integer,allocatable    ::  list_nx(:),  list_ny(:)
      real,allocatable       ::  list_lon(:), list_lat(:)

      integer                ::  ix0, iy0, nx0, ny0            !! input hires map dimention
      real                   ::  west0, north0, east0, south0

! regional hires map
      integer*2,allocatable  ::  catmXX(:,:), catmYY(:,:)
      real,allocatable       ::  flddif(:,:)

! input hires map
      integer*2,allocatable  ::  catmXX0(:,:), catmYY0(:,:)
      real,allocatable       ::  flddif0(:,:)
      real,allocatable       ::  lon0(:), lat0(:)

! files
      character*64           ::  rfile1, rfile2
      character*64           ::  wfile1, wfile2
      character*64           ::  hires
      integer                ::  ios
      character*64           ::  buf
! ===============================================
      open(11,file=region_info,form='formatted')
      read(11,*) global_dir
      close(11)

      hires=trim(global_dir)//'hires/'
      list_loc=trim(hires)//'location.txt'

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

      allocate(catmXX(nx,ny),catmYY(nx,ny),flddif(nx,ny))
      catmXX(:,:)=-9999
      catmYY(:,:)=-9999
      flddif(:,:)=-9999

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

        allocate(catmXX0(nx0,ny0),catmYY0(nx0,ny0),flddif0(nx0,ny0))
        allocate(lon0(nx0),lat0(ny0))

        rfile1=trim(hires)//trim(area)//'.catmxy'
        rfile2=trim(hires)//trim(area)//'.flddif'

        print *, rfile1
        open(21,file=rfile1,form='unformatted',access='direct',recl=2*nx0*ny0,status='old',iostat=ios)
        if( ios==0 )then
          read(21,rec=1) catmXX0
          read(21,rec=2) catmYY0
          close(21)
        else
          print *, '*******************'
          print *, 'no data: ', rfile1
          goto 1000
        endif

        print *, rfile2
        open(21,file=rfile2,form='unformatted',access='direct',recl=4*nx0*ny0,status='old',iostat=ios)
        if( ios==0 )then
          read(21,rec=1) flddif0
          close(21)
        else
          print *, '*******************'
          print *, 'no data: ', rfile2
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

              if( catmXX(ix,iy)==-9999 )then
                if( catmXX0(ix0,iy0)>0 )then
                  catmXX(ix,iy)=catmXX0(ix0,iy0)-(dXX)
                  catmYY(ix,iy)=catmYY0(ix0,iy0)-(dYY)
                  if( catmXX(ix,iy)<1 .or. catmXX(ix,iy)>nXX .or. catmYY(ix,iy)<1 .or. catmYY(ix,iy)>nYY )then
                    catmXX(ix,iy)=0
                    catmYY(ix,iy)=0
                  endif
                else
                  catmXX(ix,iy)=catmXX0(ix0,iy0)
                  catmYY(ix,iy)=catmYY0(ix0,iy0)
                endif
              endif

              if( flddif(ix,iy)==-9999 )then
                flddif(ix,iy)=flddif0(ix0,iy0)
              endif

            endif
          end do
        end do
 1000   continue
        deallocate(catmXX0,catmYY0,flddif0,lon0,lat0)
 2000   continue
      end do

      wfile1='../hires/reg.catmxy'
      open(21,file=wfile1,form='unformatted',access='direct',recl=2*nx*ny)
      write(21,rec=1) catmXX
      write(21,rec=2) catmYY
      close(21)

      wfile2='../hires/reg.flddif'
      open(21,file=wfile2,form='unformatted',access='direct',recl=4*nx*ny)
      write(21,rec=1) flddif
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

      end program COMBINE_HIRES

