! Module to defer logging to a function that can be set using the set_logger function

module logging

  use iso_c_utils
  implicit none

  abstract interface
     subroutine ilogger(level, msg)
       use iso_c_binding
       use iso_c_utils
       integer(c_int), value, intent(in) :: level !< severity
       character(c_char), intent(in) :: msg(MAXSTRINGLEN) !< c message null terminated
     end subroutine ilogger
  end interface


  procedure(ilogger), pointer :: logging_callback => null()

  ! Levels correspond to log4net and log4j

  integer, parameter, public :: LEVEL_ALL = 0
  integer, parameter, public :: LEVEL_DEBUG = 1
  integer, parameter, public :: LEVEL_INFO  = 2
  integer, parameter, public :: LEVEL_WARN  = 3
  integer, parameter, public :: LEVEL_ERROR = 4
  integer, parameter, public :: LEVEL_FATAL = 5
  integer, parameter, public :: LEVEL_OFF = 6


  ! message buffer that can be used to fill log message
  ! write(msgbuf, *) ''x=', x
  ! call log(LEVEL_DEBUG, trim(msgbuf))
  character(len=MAXSTRINGLEN) :: msgbuf

contains
  subroutine set_logger(c_callback) bind(C, name="set_logger")
    !DEC$ ATTRIBUTES DLLEXPORT::set_logger

    type(c_funptr), value :: c_callback

    ! Set a callback that will be cauled with new messages

    call c_f_procpointer(c_callback, logging_callback)
  end subroutine set_logger

  subroutine log(level, msg)
    integer(c_int), intent(in) :: level
    character(len=*), intent(in) :: msg

    character(c_char)             :: c_string(MAXSTRINGLEN)


    if (associated(logging_callback)) then
       c_string = string_to_char_array(msg)
       call logging_callback(level, c_string)
    end if
  end subroutine log


end module logging
