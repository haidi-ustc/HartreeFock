subroutine timestamp(str)

!*********************************************************************72
!
! TIMESTAMP prints out the current YMDHMS date as a timestamp.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none
  character(len=*),optional::  str
  character * ( 8 ) ampm
  integer d
  character * ( 8 ) date
  integer h
  integer m
  integer mm
  character * ( 9 ) month(12)
  integer n
  integer s
  character * ( 10 ) time
  integer y

  save month

  data month / &
   'January  ', 'February ', 'March    ', 'April    ',&
   'May      ', 'June     ', 'July     ', 'August   ',&
   'September', 'October  ', 'November ', 'December ' /

  call date_and_time ( date, time )

  read ( date, '(i4,i2,i2)' ) y, m, d
  read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

  if ( h .lt. 12 ) then
    ampm = 'AM'
  else if ( h .eq. 12 ) then
    if ( n .eq. 0 .and. s .eq. 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h .lt. 12 ) then
      ampm = 'PM'
    else if ( h .eq. 12 ) then
      if ( n .eq. 0 .and. s .eq. 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if
  if(present(str))then
   write ( *, &
     '(a20,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )&
     trim(str),d, month(m), y, h, ':', n, ':', s, '.', mm, ampm
  else 
    write ( *, &
     '(a20,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )&
     d, month(m), y, h, ':', n, ':', s, '.', mm, ampm
  end if
  return
end  subroutine
