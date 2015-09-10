!===========================================================
module rand_class
!> \brief
!! Random number generators
!<

  use type_kinds
  implicit none

  private
  public :: seed,uran,gran

  integer(long) :: idum=0_long

contains
  !----------------------------
  function seed()
    IMPLICIT NONE
    integer(long) :: seed
    seed=idum
    RETURN
  END function seed
  !----------------------------
  FUNCTION URAN(seedin)
    IMPLICIT NONE
    REAL(double):: URAN
    integer(long), optional:: seedin
    INTEGER(long),PARAMETER:: IA=16807_long
    INTEGER(long),PARAMETER:: IM=2147483647_long
    INTEGER(long),PARAMETER:: IQ=127773_long
    INTEGER(long),PARAMETER:: IR=2836_long
    INTEGER(long),PARAMETER:: MASK=20410876_long
    REAL(double),PARAMETER:: AM=1.0_double/IM
    INTEGER(long):: k

    integer(long)::t,dt(8),pid

    if(present(seedin))idum=seedin

    if(idum.EQ.0)then
       !call system_clock(t)
       !if (t == 0) then
       call date_and_time(values=dt)
       t = (dt(1) - 1970) * 365_long * 24 * 60 * 60 * 1000 &
            + dt(2) * 31_long * 24 * 60 * 60 * 1000 &
            + dt(3) * 24_long * 60 * 60 * 1000 &
            + dt(5) * 60 * 60 * 1000 &
            + dt(6) * 60 * 1000 + dt(7) * 1000 &
            + dt(8)
       pid = getpid()
       idum = ieor(t, int(pid, kind(t)))
    end if

!!$    if(idum.EQ.0)then
!!$       do idum=0, myid*100000!E5
!!$          array=timearray()
!!$       end do
!!$       idum=myid*1E4
!!$       idum=idum+array(8)
!!$       idum=idum+array(7)*1E3
!!$       idum=idum+array(6)*6E4
!!$       idum=idum+array(5)*36E5
!!$       !write(*,*)"core ",myid,"generating new random seed ",idum)
!!$    end if

    idum=IEOR(idum,MASK)
    k=idum/IQ
    idum=IA*(idum-k*IQ)-IR*k
    IF (idum.LT.0) idum=idum+IM
    URAN=AM*idum
    idum=IEOR(idum,MASK)
    RETURN
  END FUNCTION URAN
  !----------------------------
  FUNCTION GRAN(seed)
    !         FUNCTION TO GENERATE GAUSSIAN RANDOM NUMBERS
    !         USING BOX-MUELLER/WEINER METHOD.
    IMPLICIT NONE
    integer(long), optional:: seed
    REAL(double) :: GRAN,w1,w2
    real(double) , parameter :: PI=3.1415926535898_double

    if(present(seed))idum=seed

    w1=URAN()
    w2=URAN()
    gran=SQRT(-2.*LOG(w1))*COS(2.*PI*w2)

  END FUNCTION GRAN

end module rand_class
!===========================================================---------
