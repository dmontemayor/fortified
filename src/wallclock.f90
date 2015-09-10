!===========================================================
!> \brief
!! Timming Module
!! \details
!! Supplies the software with a standard method of computing
!! wall clock time in units of mils, sec, hours or even days.
!! \authors
!! Author: Daniel Montemayor
!<===========================================================
module wallclock
  use type_kinds
  implicit none
  
  private
  public:: date,time,zone,timearray
  
     character(len=8)::d
     character(len=10)::t
     character(len=5)::z
     integer(long) :: a(8)
  
contains
  !----------------------------
  function date()!<returns date in 8 char string format YYYYMMDD
    character(len=8)::date
    call date_and_time(d,t,z,a)
    date=d
  end function date
  !----------------------------
  function time()!<returns time in 10 char string format hhmmss.sss
    character(len=10)::time
    call date_and_time(d,t,z,a)
    time=t
  end function time
  !----------------------------
  function zone()!<returns zone relative to Coordinated Univeral Time (UTC) in 5 char string format (+-)hhmm
    character(len=5)::zone
    call date_and_time(d,t,z,a)
    zone=z
  end function zone
  !----------------------------
  function timearray()!<returns 8 element integer array representing time in format Yr,Mo,Day,UTCmin_dif,hr,min,sec,msec
    integer(long) :: timearray(8)
    call date_and_time(d,t,z,a)
    timearray=a
  end function timearray

end module wallclock
!===========================================================
