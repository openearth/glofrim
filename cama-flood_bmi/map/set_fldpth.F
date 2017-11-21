      program SET_FLDPTH
! ===============================================
      implicit none

! River Map (nXX,nYY)
      character*64           ::  param                            !! river map parameter
      parameter                 (param='./params.txt')

      integer                ::  nXX, nYY                         !! river map grid number
      integer                ::  iXX, iYY, jXX, jYY, kXX, kYY, lXX, lYY
      real                   ::  gsize                            !! river map grid size

      integer,allocatable    ::  lsmask(:,:)                      !! hires map ID         of grid (iXX,iYY)
      integer,allocatable    ::  nextXX(:,:)                      !! downstream (jXX,jYY) of grid(iXX,iYY)
      integer,allocatable    ::  nextYY(:,:)
      real,allocatable       ::  GRDupa(:,:)                      !! upstream area (outlet) [m2]
      real,allocatable       ::  GRDelv(:,:)                      !! elevation     (outlet) [m]
      real,allocatable       ::  nxtdst(:,:)                      !! downstream distance    [m]
      real,allocatable       ::  rivhgt(:,:)                      !! channel depth          [m]

! hires map parameters
      integer                ::  i, j, iarea, narea
      character*64           ::  list_loc                         !! hires map parameter file
      parameter                 (list_loc='./hires/location.txt')
      character*3,allocatable::  list_area(:)                     !! area tag
      integer,allocatable    ::  list_nx(:),  list_ny(:)          !! grid numbers
      real,allocatable       ::  list_lon(:), list_lat(:)         !! west, north edge

! Hi-res map map "extended (mx,my)" <- "original(nx,ny)"
      character*3            ::  area                             !! area tag (sa1, ca1, na1, af1, eu1, eu2, as1, as2, as3, oc1, na2, eu3, si2, si2)
      integer                ::  nx, ny, mx, my                   !! hires map (domain+margin)
      integer                ::  ix, iy, jx, jy, kx, ky, dx, dy

      real                   ::  csize                            !! hires map pixel size
      real                   ::  lon_ori, lat_ori                 !! west, north edge [deg]
      real                   ::  lon_end, lat_end                 !! east, south edge [deg]
      real                   ::  margin                           !! extended margin [deg]
      parameter                 (margin=5.0)

      integer*2,allocatable  ::  nextx(:,:),  nexty(:,:)          !! downstream (ix,iy)
      integer*2,allocatable  ::  catmXX(:,:), catmYY(:,:)         !! grid(iXX,iYY) of pixel(ix,iy)
      real,allocatable       ::  elevtn(:,:)                      !! elevation [m]
      real,allocatable       ::  uparea(:,:)                      !! drainage area [km2]
      integer*1,allocatable  ::  water(:,:)                       !! water mask
      real,allocatable       ::  lon(:), lat(:)
! calculation
      real,allocatable       ::  chndst(:,:)                      !! distance to river channel [m]
      integer,allocatable    ::  adjacXX(:,:), adjacYY(:,:)       !! bifurcation adjacent catchment
      real,allocatable       ::  pthdst(:,:)                      !! bifurcation pathway distance
      integer*2,allocatable  ::  blayer(:,:)                      !! bifurcation pathway layer (1:river, >2: overland)

      integer,allocatable    ::  check(:,:)                       !! to calculate river bifurcation depth

      integer                ::  elvdif                           !! floodplain height above river channel
      integer                ::  adjac_jx, adjac_jy               !! bifurcation adjacent pixel (jx,jy)
      real                   ::  dst, dst_this                    !! distance calculation
      real                   ::  elv_min                          !! minimum elevation among adjacent pixels in different unit-catchment
      real                   ::  cwidth                           !! pixel width (for pathway width calculation)


! parameters
      real                   ::  uparea_min
      real                   ::  pth_dif
      parameter                 (uparea_min=100.)                 !! threshold for channel definition [km2]
      parameter                 (pth_dif=3.)                      !! maximum bifurcation height  [m]

      real                   ::  bifu_hgt0, bifu_hgt1, bifu_hgt2
      real                   ::  brivdst
      real                   ::  brivhgt
      real                   ::  bhgtmod
      parameter                 (brivdst=100000.)                 !! search downstream untill 100km     to decide bifurcation river depth
      parameter                 (brivhgt=3.0)                      !! search downstream untill 3m rivhgt to decide bifurcation river depth
      parameter                 (bhgtmod=0.5)                      !! bifurcation channel depth, modification ratio

! bifurcation channel topography
      integer                ::  ipath, npath, ipath_max
      integer                ::  ilev, nlev                       !! bifurcation channel layer
      parameter                 (nlev=int(pth_dif)+2)             !! bifurcation channel max layer

      integer,allocatable    ::  upst_XX(:), upst_YY(:)           !! pathway upstream downstream grid
      integer,allocatable    ::  down_XX(:), down_YY(:)

      real,allocatable       ::  bifu_dst(:)                      !! bifurcation channel distance  [m]
      real,allocatable       ::  bifu_elv(:,:)                    !! bifurcation channel elevation [m]
      real,allocatable       ::  bifu_wth(:,:)                    !! bifurcation channel width     [m]

      real,allocatable       ::  bifu_lon(:)                      !! bifurcation channel lon lat
      real,allocatable       ::  bifu_lat(:)

      real,allocatable       ::  bifu_num(:)                      !! bifurcation point number (only for calculation)
      integer,allocatable    ::  bifu_minlev(:)                   !! lowest bifurcation layer (only for calculation)
! file
      character*128          ::  clsmask, crivmap, cgrdupa, cgrdelv, cnxtdst, crivhgt
      character*128          ::  wfile1, wfile2
      character*128          ::  clen, cfmt
      character*128          ::  buf
! function
      real                   ::  rgetlen
! ===============================================
print *, 'SET_FLDPTH:'

      open(11,file=param,form='formatted')
      read(11,*)          !! west
      read(11,*)          !! north
      read(11,*) nXX
      read(11,*) nYY
      read(11,*) gsize
      read(11,*)          !! nflp
      read(11,*) narea
      read(11,*) csize
      close(11)

      allocate(lsmask(nXX,nYY),nextXX(nXX,nYY),nextYY(nXX,nYY))
      allocate(GRDupa(nXX,nYY),GRDelv(nXX,nYY),nxtdst(nXX,nYY),rivhgt(nXX,nYY))
      allocate(list_area(narea),list_nx(narea),list_ny(narea),list_lon(narea),list_lat(narea))
      allocate(check(nXX,nYY))

      npath=nXX*nYY
      allocate(upst_XX(npath), upst_YY(npath),down_XX(npath),down_YY(npath))
      allocate(bifu_dst(npath),bifu_elv(npath,nlev),bifu_wth(npath,nlev))
      allocate(bifu_lon(npath),bifu_lat(npath),bifu_num(npath),bifu_minlev(npath))



      clsmask='./lsmask.bin'
      crivmap='./nextxy.bin'
      cgrdupa='./uparea.bin'
      cgrdelv='./elevtn.bin'
      cnxtdst='./nxtdst.bin'
      crivhgt='./rivhgt.bin'

      wfile1='./fldpth.txt'
      wfile2='./pathll.txt'

      print *, '  reading lsmask: ', clsmask
      open(11, file=clsmask, form='unformatted', access='direct', recl=4*nXX*nYY)
      read(11,rec=1) lsmask
      close(11)

      print *, '  reading river netwrok map: ', crivmap
      open(11, file=crivmap, form='unformatted', access='direct', recl=4*nXX*nYY)
      read(11,rec=1) nextXX
      read(11,rec=2) nextYY
      close(11)

      print *, '  reading drainage area: ', cgrdupa
      open(11, file=cgrdupa, form='unformatted', access='direct', recl=4*nXX*nYY)
      read(11,rec=1) GRDupa
      close(11)

      print *, '  reading drainage area: ', cgrdelv
      open(11, file=cgrdelv, form='unformatted', access='direct', recl=4*nXX*nYY)
      read(11,rec=1) GRDelv
      close(11)

      print *, '  reading drainage area: ', cnxtdst
      open(11, file=cnxtdst, form='unformatted', access='direct', recl=4*nXX*nYY)
      read(11,rec=1) nxtdst
      close(11)

      print *, '  reading drainage area: ', crivhgt
      open(11, file=crivhgt, form='unformatted', access='direct', recl=4*nXX*nYY)
      read(11,rec=1) rivhgt
      close(11)

      open(21,file=wfile2,form='formatted')

      ipath_max=0
      upst_XX(:)=-999
      upst_YY(:)=-999
      down_XX(:)=-999
      down_YY(:)=-999

      bifu_dst(:)=0
      bifu_elv(:,:)=9999
      bifu_wth(:,:)=0

      bifu_num(:)=0
      bifu_lon(:)=-999
      bifu_lat(:)=-999
      bifu_minlev(:)=999

! ===============================================
      open(11,file=list_loc,form='formatted')
      read(11,*)
      read(11,*) buf, (list_area(i) ,i=1,narea)
      read(11,*) buf, (list_lon(i)  ,i=1,narea)
      read(11,*) buf, (list_lat(i)  ,i=1,narea)
      read(11,*) buf, (list_nx(i)   ,i=1,narea)
      read(11,*) buf, (list_ny(i)   ,i=1,narea)
      close(11)

! ===============================================
      do i=1, narea                              !! set up bifurcation channel for each area
        area=list_area(i)
        iarea=i
        nx=list_nx(i)
        ny=list_ny(i)
        lon_ori=list_lon(i)
        lat_ori=list_lat(i)

        lon_end=lon_ori+real(nx)*csize
        lat_end=lat_ori-real(ny)*csize

        lon_ori=lon_ori-margin                  !! set margin for each area
        lon_end=lon_end+margin
        lat_ori=lat_ori+margin
        lat_end=lat_end-margin
        mx=nx+int(margin*2./csize)
        my=ny+int(margin*2./csize)

        print *, '  **********'
        print *, '    calculation for area: ', area

        allocate(nextx(mx,my),nexty(mx,my))
        allocate(catmXX(mx,my),catmYY(mx,my))
        allocate(elevtn(mx,my),uparea(mx,my),water(mx,my))
        allocate(lon(mx),lat(my))

        allocate(chndst(mx,my))
        allocate(adjacXX(mx,my),adjacYY(mx,my),pthdst(mx,my),blayer(mx,my))

        CALL SET_HIRES_MAP  !! set hires map for area

! ===============================================
print *, '  *** main calculation: ', area
print *, '    define river channel pixels'     !! set channel distance chndst=0 for main river channels 

        chndst(:,:)=-9999

        do iy=1, my                   !! outlet pixel and river mouth
          do ix=1, mx
            if( nextx(ix,iy)>0 .and. catmXX(ix,iy)>0 )then
              jx=nextx(ix,iy)
              jy=nexty(ix,iy)
              if( catmXX(ix,iy)/=catmXX(jx,jy) .or. catmYY(ix,iy)/=catmYY(jx,jy) )then
                chndst(ix,iy)=0
              endif
            elseif( nextx(ix,iy)<0 .and. nextx(ix,iy)/=-9999 )then
              chndst(ix,iy)=0
            endif
          end do
        end do

        do iy=1, my                   !! downstrem of outlet pixel
          do ix=1, mx
            if( chndst(ix,iy)==0 .and. nextx(ix,iy)>0 .and. catmXX(ix,iy)>0 )then
              jx=nextx(ix,iy)
              jy=nexty(ix,iy)
              do while( chndst(jx,jy)==-9999 .and. nextx(jx,jy)>0 )
                chndst(jx,jy)=0
                kx=nextx(jx,jy)
                ky=nexty(jx,jy)
                jx=kx
                jy=ky
              end do
            endif
          end do
        end do

        do iy=1, my                   !! large channels
          do ix=1, mx
            if( uparea(ix,iy)>uparea_min ) chndst(ix,iy)=0
          end do
        end do
! =====
print *, '    calculate distance to river channel'

        do iy=1, my
          do ix=1, mx
            if( chndst(ix,iy)==-9999 .and. nextx(ix,iy)>0 .and. catmXX(ix,iy)>0 )then
              dst=0
              jx=ix
              jy=iy
              do while( chndst(jx,jy)==-9999 .and. nextx(jx,jy)>0 )
                kx=nextx(jx,jy)
                ky=nexty(jx,jy)
                dst=dst+rgetlen(lon(jx),lat(jy),lon(kx),lat(ky))
                jx=kx
                jy=ky
              end do
              dst_this=dst+chndst(jx,jy)

              jx=ix
              jy=iy
              do while( chndst(jx,jy)==-9999 )
                chndst(jx,jy)=dst_this
                kx=nextx(jx,jy)
                ky=nexty(jx,jy)
                dst_this=dst_this-rgetlen(lon(jx),lat(jy),lon(kx),lat(ky))
                jx=kx
                jy=ky
              end do
            endif
          end do
        end do

! ===============================================
print *, '    calculate boundary'
        adjacXX(:,:)=-9999
        adjacYY(:,:)=-9999
        pthdst(:,:)=-9999
        blayer(:,:)=0

        do iy=1, my
          do ix=1, mx
            if( catmXX(ix,iy)>0 )then
              iXX=catmXX(ix,iy)
              iYY=catmYY(ix,iy)
              elvdif=int(elevtn(ix,iy)-GRDelv(iXX,iYY))
              adjac_jx=-9
              adjac_jy=-9

              if( elvdif>=0 .and. elvdif<=pth_dif )then   !!  potential bifurcation point

                elv_min=1.e20
                do dy=-1, 1       !! search neighboring pixel in different unit-catchment with lowest elevation
                  do dx=-1, 1
                    if( abs(dx)+abs(dy)>=1 ) then
                      jx=ix+dx
                      jy=iy+dy
                      if( jx<=0 .or. jx>mx .or. jy<=0 .or. jy>my ) goto 1100
                      jXX=catmXX(jx,jy)
                      jYY=catmYY(jx,jy)
                                                                                      !! reject in the cases below
                      if( iXX<=0 .or. iXX>nXX .or. iYY<=0 .or. iYY>nYY     ) goto 1100
                      if( jXX<=0 .or. jXX>nXX .or. jYY<=0 .or. jYY>nYY     ) goto 1100
                      if( lsmask(iXX,iYY)>iarea .or. lsmask(jXX,jYY)>iarea ) goto 1100
                      if( catmXX(jx,jy)<=0           ) goto 1100

                      if( (iXX/=jXX .or. iYY/=jYY) .and. jXX>0 )then   !! if neighbouring pixels are in different grids

                        if( elevtn(ix,iy)<elevtn(jx,jy) )                     goto 1000        !! decide flow pathway downstream (jXX,jYY)

                        if( nextXX(iXX,iYY)==jXX .and. nextYY(iXX,iYY)==jYY ) goto 1000        !! if downstream, reject
                        if( nextXX(jXX,jYY)==iXX .and. nextYY(jXX,jYY)==iYY ) goto 1000

                        if( GRDelv(jXX,jYY)>GRDelv(iXX,iYY) )then
                          elvdif=int(elevtn(ix,iy)-GRDelv(jXX,jYY))
                        endif
                        if( elvdif>pth_dif ) goto 1000

                        if( elevtn(ix,iy)==elevtn(jx,jy) )then                !! if elevation is same, decide upstream/downstream by GRDupa or (iXX,iYY)
                          if(  GRDupa(iXX,iYY) < GRDupa(jXX,jYY) ) goto 1000
                          if( (GRDupa(iXX,iYY)==GRDupa(jXX,jYY) .and. iXX+nXX*iYY>jXX+nXX*jYY ) ) goto 1000
                        endif

                  !! === decide bifurcation flow direction toward lowest neighbouting pixel in different unit-catchment

                        if( elevtn(jx,jy)<elv_min )then
                          elv_min=elevtn(jx,jy)
                          adjacXX(ix,iy)=jXX
                          adjacYY(ix,iy)=jYY
                          pthdst(ix,iy)=chndst(ix,iy)+chndst(jx,jy)+rgetlen(lon(ix),lat(iy),lon(jx),lat(jy))
                          adjac_jx=jx
                          adjac_jy=jy
                          if( water(ix,iy)==1 .and. water(jx,jy)==1 ) blayer(ix,iy)=1
                        endif

                        if( blayer(ix,iy)==0 .and. water(ix,iy)==1 .and. water(jx,jy)==1 )then
                          elv_min=elevtn(jx,jy)
                          adjacXX(ix,iy)=jXX
                          adjacYY(ix,iy)=jYY
                          pthdst(ix,iy)=chndst(ix,iy)+chndst(jx,jy)+rgetlen(lon(ix),lat(iy),lon(jx),lat(jy))
                          adjac_jx=jx
                          adjac_jy=jy
                          if( water(ix,iy)==1 .and. water(jx,jy)==1 ) blayer(ix,iy)=1
                        endif

 1000                   continue
                      endif
 1100                 continue
                    endif
                  end do
                end do

                if( pthdst(ix,iy)>0 )then               !! if pixel(ix,iy) is bifurcation point
                  if( blayer(ix,iy)==0 )then              !! blayer=1 (layer=1): river bifurcation pathway
                    blayer(ix,iy)=elvdif+2                  !! for overland bifurcation, (layer=elvdif+ 2)
                  endif
                  jx=adjac_jx
                  jy=adjac_jy
                  write(21,*) ix, iy, jx, jy, lon(ix), lat(iy), lon(jx), lat(jy), elevtn(ix,iy), blayer(ix,iy)
                endif

              endif

            end if
          end do
        end do

print *, '    set floodplain path'

        check(:,:)=0
        do iy=1, my
          do ix=1, mx
            if( adjacXX(ix,iy)>0 )then
              iXX=catmXX(ix,iy)
              iYY=catmYY(ix,iy)
              jXX=adjacXX(ix,iy)
              jYY=adjacYY(ix,iy)
              if( GRDupa(iXX,iYY) <GRDupa(jXX,jYY) .or. &                                  !! reverse upstream & downstream relation
                 (GRDupa(iXX,iYY)==GRDupa(jXX,jYY) .and. iXX+nXX*iYY>jXX+nXX*jYY ) )then
                iXX=adjacXX(ix,iy)
                iYY=adjacYY(ix,iy)
                jXX=catmXX(ix,iy)
                jYY=catmYY(ix,iy)
              endif

              do ipath=1, npath
                if( iXX==upst_XX(ipath) .and. iYY==upst_YY(ipath) .and. jXX==down_XX(ipath) .and. jYY==down_YY(ipath) )then

                  bifu_dst(ipath)=bifu_dst(ipath)+pthdst(ix,iy)
                  bifu_num(ipath)=bifu_num(ipath)+1.

                  ilev=blayer(ix,iy)

                  cwidth=( rgetlen(0.,lat(iy),csize,lat(iy))+rgetlen(0.,lat(iy),0.,lat(iy)+csize) ) * 0.5
                  bifu_wth(ipath,ilev)=bifu_wth(ipath,ilev)+cwidth

                  if( bifu_elv(ipath,ilev)==9999 )then
                    if( ilev==1 ) then  !! for river pathway

                      kXX=iXX
                      kYY=iYY
                      dst=0
                      do while( dst<brivdst .and. rivhgt(kXX,kYY)<brivhgt .and. nextXX(kXX,kYY)>0 )
                        check(kXX,kYY)=ipath
                        lXX=nextXX(kXX,kYY)
                        lYY=nextYY(kXX,kYY)
                        dst=dst+nxtdst(kXX,kYY)
                        kXX=lXX
                        kYY=lYY
                      end do
                      check(kXX,kYY)=ipath
                      bifu_hgt1=rivhgt(kXX,kYY)

                      kXX=jXX
                      kYY=jYY
                      dst=0
                      do while( dst<brivdst .and. rivhgt(kXX,kYY)<brivhgt .and. nextXX(kXX,kYY)>0 .and. check(kXX,kYY)/=ipath )
                        lXX=nextXX(kXX,kYY)
                        lYY=nextYY(kXX,kYY)
                        dst=dst+nxtdst(kXX,kYY)
                        kXX=lXX
                        kYY=lYY
                      end do
                      bifu_hgt2=rivhgt(kXX,kYY)

                      if( check(kXX,kYY)==ipath )then
                        bifu_hgt0=rivhgt(kXX,kYY)
                      else
                        bifu_hgt0=max(bifu_hgt1,bifu_hgt2)
                      endif
                      bifu_elv(ipath,ilev)=max(GRDelv(iXX,iYY),GRDelv(jXX,jYY))-bifu_hgt0*bhgtmod

                    else   !! for overland pathway

                      bifu_elv(ipath,ilev)=elevtn(ix,iy)
                    endif
                  endif

                  if( ilev<bifu_minlev(ipath) )then
                    bifu_lon(ipath)=lon(ix)
                    bifu_lat(ipath)=lat(iy)
                    bifu_minlev(ipath)=ilev
                  endif
                  goto 2000
                endif

                if( upst_XX(ipath)==-999 )then
                  ipath_max=ipath_max+1

                  upst_XX(ipath)=iXX
                  upst_YY(ipath)=iYY
                  down_XX(ipath)=jXX
                  down_YY(ipath)=jYY

                  bifu_dst(ipath)=pthdst(ix,iy)
                  bifu_num(ipath)=1.

                  ilev=blayer(ix,iy)
                  cwidth=( rgetlen(0.,lat(iy),csize,lat(iy))+rgetlen(0.,lat(iy),0.,lat(iy)+csize) ) * 0.5
                  bifu_wth(ipath,ilev)=cwidth
                  print *, ilev
                  if( ilev==1 ) then  !! for river pathway

                    kXX=iXX
                    kYY=iYY
                    dst=0
                    do while( dst<brivdst .and. rivhgt(kXX,kYY)<brivhgt .and. nextXX(kXX,kYY)>0 )
                      check(kXX,kYY)=ipath
                      lXX=nextXX(kXX,kYY)
                      lYY=nextYY(kXX,kYY)
                      dst=dst+nxtdst(kXX,kYY)
                      kXX=lXX
                      kYY=lYY
                    end do
                    check(kXX,kYY)=ipath
                    bifu_hgt1=rivhgt(kXX,kYY)

                    kXX=jXX
                    kYY=jYY
                    dst=0
                    do while( dst<brivdst .and. rivhgt(kXX,kYY)<brivhgt .and. nextXX(kXX,kYY)>0 .and. check(kXX,kYY)/=ipath )
                      lXX=nextXX(kXX,kYY)
                      lYY=nextYY(kXX,kYY)
                      dst=dst+nxtdst(kXX,kYY)
                      kXX=lXX
                      kYY=lYY
                    end do
                    bifu_hgt2=rivhgt(kXX,kYY)

                    if( check(kXX,kYY)==ipath )then
                      bifu_hgt0=rivhgt(kXX,kYY)
                    else
                      bifu_hgt0=max(bifu_hgt1,bifu_hgt2)
                    endif
                    bifu_elv(ipath,ilev)=max(GRDelv(iXX,iYY),GRDelv(jXX,jYY))-bifu_hgt0*bhgtmod

                  else   !! for overland pathway

                    bifu_elv(ipath,ilev)=elevtn(ix,iy)
                  endif

                  bifu_lon(ipath)=lon(ix)
                  bifu_lat(ipath)=lat(iy)
                  bifu_minlev(ipath)=ilev

                  goto 2000
                endif
                if( ipath==npath )then
                  print *, 'ipath > ', npath
                  stop
                endif
              end do
 2000         continue
            endif
          end do
        end do

        deallocate(nextx, nexty,catmXX,catmYY,elevtn,uparea,water,lon,lat)
        deallocate(chndst,adjacXX,adjacYY,pthdst,blayer)

      end do

      close(21)

! ===================================================
rivhgt(:,:)=-9999

      print *, 'write output: ', wfile1
      print *, '  number of flood path =', ipath_max

      open(31,file=wfile1,form='formatted')
      write(31,'(2i8,a)') ipath_max, nlev, &
               '   npath, npathlev, (ix,iy), (jx,jy), length, (elevtn1, width1), (elevtn2, width2), ... , (lon,lat)'

      do ipath=1, ipath_max
        bifu_dst(ipath)=bifu_dst(ipath)/bifu_num(ipath)

        iXX=upst_XX(ipath)
        iYY=upst_YY(ipath)
        jXX=down_XX(ipath)
        jYY=down_YY(ipath)

if( bifu_elv(ipath,1)/=9999. )then
rivhgt(iXX,iYY)=max( max(GRDelv(iXX,iYY),GRDelv(jXX,jYY)) -bifu_elv(ipath,1), rivhgt(iXX,iYY))
endif
        write(clen,'(i2)') 1+2*nlev
        cfmt='(4i8,'//trim(clen)//'f12.2,2f10.3)'

        write(31,cfmt), upst_XX(ipath), upst_YY(ipath), down_XX(ipath), down_YY(ipath), &
                        bifu_dst(ipath),(bifu_elv(ipath,ilev),bifu_wth(ipath,ilev),ilev=1,nlev), &
                        bifu_lon(ipath),bifu_lat(ipath)
      end do
      close(31)

open(11, file='tmp.bin', form='unformatted', access='direct', recl=4*nXX*nYY)
write(11,rec=1) rivhgt
close(11)

! ==================================================


! ====================================================================
      CONTAINS

! ====================================================================
      subroutine SET_HIRES_MAP

      implicit none

! Hi-res map input (nx2,ny2)
      character*3            ::  area2
      integer                ::  nx2, ny2
      real                   ::  lon_ori2, lat_ori2               !! west, north edge [deg]
      real                   ::  lon_end2, lat_end2               !! east, south edge [deg]

      integer*2,allocatable  ::  nextx2(:,:),  nexty2(:,:)        !! downstream (ix,iy)
      integer*2,allocatable  ::  catmXX2(:,:), catmYY2(:,:)       !! grid(iXX,iYY) of pixel(ix,iy)
      real,allocatable       ::  elevtn2(:,:)                     !! elevation [m]
      real,allocatable       ::  uparea2(:,:)                        !! drainage area [km2]
      real,allocatable       ::  rivwth2(:,:)                     !! river width (water mask) [m]
      real,allocatable       ::  lon2(:), lat2(:)

      character*128          ::  SHEDSdir, HIRESdir
      parameter                 (SHEDSdir='./sheds/')
      parameter                 (HIRESdir='./hires/')
      character*128          ::  rfile1, rfile2, rfile3, rfile4, rfile5
! ===============================================
        print *, '    read hires maps for area: ', area

        catmXX(:,:)=-9999
        catmYY(:,:)=-9999
        nextx(:,:)=-9999
        nexty(:,:)=-9999
        elevtn(:,:)=-9999
        water(:,:)=0
        uparea(:,:)=-9999

        do ix=1, mx
          lon(ix)=lon_ori+csize*(real(ix)-0.5)
        end do
        do iy=1, my
          lat(iy)=lat_ori-csize*(real(iy)-0.5)
        end do

        do j=1, iarea                            !! combine hires maps within each area (+margin)

          area2=list_area(j)
          nx2=list_nx(j)
          ny2=list_ny(j)
          lon_ori2=list_lon(j)
          lat_ori2=list_lat(j)
          lon_end2=lon_ori2+real(nx2)*csize
          lat_end2=lat_ori2-real(ny2)*csize

          if( lon_ori2>lon_end .or. lon_end2<lon_ori .or. lat_ori2<lat_end .or. lat_end2>lat_ori )then
            print *, '      out of domain: ', trim(area2)
          else
            print *, '      read hires maps: ', trim(area2)

            rfile1=trim(HIRESdir)//trim(area2)//'.catmxy'
            rfile2=trim(SHEDSdir)//trim(area2)//'.elevtn.flt'
            rfile3=trim(SHEDSdir)//trim(area2)//'.nextxy.bsq'
            rfile4=trim(SHEDSdir)//trim(area2)//'.uparea.flt'
            rfile5=trim(SHEDSdir)//trim(area2)//'.rivwth.flt'

            allocate(nextx2(nx2,ny2), nexty2(nx2,ny2))
            allocate(catmXX2(nx2,ny2),catmYY2(nx2,ny2))
            allocate(elevtn2(nx2,ny2),uparea2(nx2,ny2),rivwth2(nx2,ny2))
            allocate(lon2(nx2),       lat2(ny2))

            print *, '      - reading catchment xy: ', rfile1
            open(12, file=rfile1, form='unformatted', access='direct', recl=2*nx2*ny2)
            read(12,rec=1) catmXX2
            read(12,rec=2) catmYY2
            close(12)

            print *, '      - reading elevation: ', rfile2
            open(11, file=rfile2, form='unformatted', access='direct', recl=4*nx2*ny2)
            read(11,rec=1) elevtn2
            close(11)

            print *, '      - reading nextxy: ', rfile3
            open(11, file=rfile3, form='unformatted', access='direct', recl=2*nx2*ny2)
            read(11,rec=1) nextx2
            read(11,rec=2) nexty2
            close(11)

            print *, '      - reading elevation: ', rfile4
            open(11, file=rfile4, form='unformatted', access='direct', recl=4*nx2*ny2)
            read(11,rec=1) uparea2
            close(11)

            print *, '      - reading river width: ', rfile5
            open(11, file=rfile5, form='unformatted', access='direct', recl=4*nx2*ny2)
            read(11,rec=1) rivwth2
            close(11)

            do ix=1, nx2
              lon2(ix)=lon_ori2+csize*(real(ix)-0.5)
            end do
            do iy=1, ny2
              lat2(iy)=lat_ori2-csize*(real(iy)-0.5)
            end do

! ===================

            do iy=1, ny2
              do ix=1, nx2
                if( nextx2(ix,iy)/=-9999 )then
                  jx=int( (lon2(ix)-lon_ori)/csize )+1
                  jy=int( (lat_ori-lat2(iy))/csize )+1
                  if( jx>=1 .and. jx<=mx .and. jy>=1 .and. jy<=my )then
                    if( nextx2(ix,iy)>0 )then
                      dx=nextx2(ix,iy)-ix
                      dy=nexty2(ix,iy)-iy
                      nextx(jx,jy)=jx+dx
                      nexty(jx,jy)=jy+dy
                      if( nextx(jx,jy)<1 .or. nextx(jx,jy)>mx .or. nexty(jx,jy)<1 .or. nexty(jx,jy)>my )then
                        nextx(jx,jy)=-9
                        nexty(jx,jy)=-9
                      endif
                    else
                      nextx(jx,jy)=nextx2(ix,iy)
                      nexty(jx,jy)=nexty2(ix,iy)
                    endif
                    catmXX(jx,jy)=catmXX2(ix,iy)
                    catmYY(jx,jy)=catmYY2(ix,iy)
                    elevtn(jx,jy)=elevtn2(ix,iy)
                    uparea(jx,jy)=uparea2(ix,iy)
                    if( rivwth2(ix,iy)>0 .or. rivwth2(ix,iy)==-1 )then
                      water(jx,jy)=1
                    endif
                  endif
                endif
              end do
            end do

            deallocate(nextx2, nexty2)
            deallocate(catmXX2,catmYY2)
            deallocate(elevtn2,uparea2,rivwth2)
            deallocate(lon2,lat2)
          endif

        end do

print *, '    done: read hires maps'

      END SUBROUTINE SET_HIRES_MAP
! ====================================================================
      END PROGRAM SET_FLDPTH






      real function rgetlen(rlon1, rlat1, rlon2, rlat2)
! ================================================
! to   get the length (m) between (rlon1, rlat1) to (rlon2, rlat2)
! by   nhanasaki
! on   1st Nov 2003
! at   IIS,UT
!
!     see page 643 of Rika-Nenpyo (2000)
!     at the final calculation, earth is assumed to be a sphere
! ================================================
      implicit none
      real                ::  rpi                !! Pi
      double precision    ::  de2                !! eccentricity powered by 2
      double precision    ::  da                 !! the radius of the earth
!
      real                ::  rlon1              !! longitude of the origin
      real                ::  rlon2              !! longitude of the destination
      real                ::  rlat1              !! latitude of the origin
      real                ::  rlat2              !! latitude of the destination
      double precision    ::  dsinlat1           !! sin(lat1)
      double precision    ::  dsinlon1           !! sin(lon1)
      double precision    ::  dcoslat1           !! cos(lat1)
      double precision    ::  dcoslon1           !! cos(lon1)
      double precision    ::  dsinlat2           !! sin(lat2) 
      double precision    ::  dsinlon2           !! sin(lon2)
      double precision    ::  dcoslat2           !! cos(lat2)
      double precision    ::  dcoslon2           !! cos(lon2)
      double precision    ::  dh1                !! hegiht of the origin
      double precision    ::  dn1                !! intermediate val of calculation
      double precision    ::  dx1                !! X coordinate of the origin
      double precision    ::  dy1                !! Y coordinate of the origin
      double precision    ::  dz1                !! Z coordinate of the origin
      double precision    ::  dh2                !! height of the destination
      double precision    ::  dn2                !! intermediate val of calculation
      double precision    ::  dx2                !! X coordinate of the destination
      double precision    ::  dy2                !! Y coordinate of the destination
      double precision    ::  dz2                !! Z coordinate of the destination
!
      double precision    ::  dlen               !! length between origin and destination
      double precision    ::  drad               !! half of the angle
! parameters
      data             da/6378137.0/
      data             de2/0.006694470/
      data             rpi/3.141592/      
! ================================================
! (lon1,lat1) --> (x1,y1,z1)
! ================================================
      dh1=0
      dh2=0

      dsinlat1 = dble(sin(rlat1 * rpi/180))
      dsinlon1 = dble(sin(rlon1 * rpi/180))
      dcoslat1 = dble(cos(rlat1 * rpi/180))
      dcoslon1 = dble(cos(rlon1 * rpi/180))
!
      dn1 = da/(sqrt(1.0-de2*dsinlat1*dsinlat1))
      dx1 = (dn1+dh1)*dcoslat1*dcoslon1
      dy1 = (dn1+dh1)*dcoslat1*dsinlon1
      dz1 = (dn1*(1-de2)+dh1)*dsinlat1
! ================================================
! (lon2,lat2) --> (x2,y2,z2)
! ================================================
      dsinlat2 = dble(sin(rlat2 * rpi/180))
      dsinlon2 = dble(sin(rlon2 * rpi/180))
      dcoslat2 = dble(cos(rlat2 * rpi/180))
      dcoslon2 = dble(cos(rlon2 * rpi/180))
!
      dn2 = da/(sqrt(1.0-de2*dsinlat2*dsinlat2))
      dx2 = (dn2+dh2)*dcoslat2*dcoslon2
      dy2 = (dn2+dh2)*dcoslat2*dsinlon2
      dz2 = (dn2*(1-de2)+dh2)*dsinlat2      
! ================================================
! Calculate length
! ================================================
      dlen=sqrt((dx1-dx2)**2+(dy1-dy2)**2+(dz1-dz2)**2)
      drad=dble(asin(real(dlen/2/da)))
      rgetlen=real(drad*2*da)
!
      return
      end function rgetlen
