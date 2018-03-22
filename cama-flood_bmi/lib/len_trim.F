      integer function len_trim(char)
! ================================================
! to get length of characters
! ================================================
      character*(*)       ::  char
      integer             ::  i
!
      do i=1,len(char)
        if(char(i:i).eq.' ')then
          len_trim=i-1
          return
        end if
      end do
      len_trim=len(char)
!
      return
      end function len_trim
