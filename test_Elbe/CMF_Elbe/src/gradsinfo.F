      program wrte_ctl
! ===============================================
      implicit none
      character*128      ::  params
      parameter             (params='../params.txt')

      integer            ::  nXX, nYY                      !! river netwrok domain (nXX,nYY)
      real               ::  gsize                         !! grid size [degree]
      real               ::  west, east, north, south

      real               ::  lon_ori, lat_ori
! ===============================================
      open(11,file=params,form='formatted')
      read(11,*) west
      read(11,*) north
      read(11,*) nXX
      read(11,*) nYY
      read(11,*) gsize
      close(11)

      east =west +real(nXX)*gsize
      south=north-real(nYY)*gsize

      lon_ori=west +gsize*0.5
      lat_ori=south+gsize*0.5

      write(*,'(2i6,f10.6,2f10.3)') nXX, nYY, gsize, lon_ori, lat_ori

      end program wrte_ctl
