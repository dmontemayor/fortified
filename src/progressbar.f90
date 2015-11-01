module progressbar
  implicit none
  private

  public::progress

  integer::icurr=0,imax=30
  logical::initialized=.false.

CONTAINS
  subroutine progress(current,end)
    integer,intent(in)::current,end
    icurr=current
    imax=end

    call delete_bar()
    call print_bar()
  end subroutine progress

  SUBROUTINE print_bar
    character(len=1) :: bar= '='
    integer::k
    if(.not.initialized)initialized=.true.
    ! print the percentage and the bar
    write(6,'(2x,1i3,1a1,2x,1a1,256a1)', advance='no') &
         100*icurr/imax,'%','|', (bar, k =1,50*icurr/imax)
    !close(6)
    !open(6)
    if(icurr.GE.imax)then
       initialized=.false.
       write(6,'(a)') '| done.'
    end if
  END SUBROUTINE print_bar

  SUBROUTINE delete_bar
    character(len=1) :: back = char(8)
    integer::k
    ! delete the bar and the percentage
    if(initialized)&
         write(6,'(256a1)', advance='no') (back, k =1,(50*icurr/imax)+9)
  END SUBROUTINE delete_bar
  
END module PROGRESSBAR
