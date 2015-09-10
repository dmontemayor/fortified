module testing_class
  use type_kinds
  implicit none
  private

  public::assert

  interface assert
     module procedure assert_logical
     module procedure assert_shortinteq
     module procedure assert_longinteq
     module procedure assert_singlerealeq
     module procedure assert_doublerealeq
     !module procedure assert_quadrealeq
     module procedure assert_streq
  end interface

contains
  subroutine assert_logical(cond,msg,iostat)
    logical,intent(in)::cond
    character(len=*),optional,intent(in)::msg
    integer(short),optional,intent(out)::iostat
    if(present(iostat))iostat=-1
    if(cond)then
       if(present(iostat))iostat=0
       return
    else
       if(present(iostat))iostat=1
       write(*,*)msg
       if(.not.present(iostat))stop
       return   
    end if
    stop

  end subroutine assert_logical
!-------------------------
  subroutine assert_shortinteq(a,b,msg,iostat)
    integer(short),intent(in)::a,b
    character(len=*),optional,intent(in)::msg
    integer(short),optional,intent(out)::iostat
    logical::inputOK=.true.
    if(present(iostat))iostat=-1   
    if(abs(a).GE.huge(a))inputOK=.false.
    if(abs(b).GE.huge(b))inputOK=.false.
    if(a.NE.a)inputOK=.false.
    if(b.NE.b)inputOK=.false.
    if(a.EQ.b.and.inputOK)then
       if(present(iostat))iostat=0
       return
    else
       if(present(iostat).and.inputOK)iostat=1
       write(*,*)msg
       if(.not.present(iostat))stop
       return   
    end if
    stop
 end subroutine assert_shortinteq
!-------------------------
  subroutine assert_longinteq(a,b,msg,iostat)
    integer(long),intent(in)::a,b
    character(len=*),optional,intent(in)::msg
    integer(short),optional,intent(out)::iostat
    logical::inputOK=.true.
    if(present(iostat))iostat=-1
    if(abs(a).GE.huge(a))inputOK=.false.
    if(abs(b).GE.huge(b))inputOK=.false.
    if(a.ne.a)inputOK=.false.
    if(b.ne.b)inputOK=.false.
    if(a.EQ.b.and.inputOK)then
       if(present(iostat))iostat=0
       return
    else
       write(*,*)msg
       if(present(iostat).and.inputOK)iostat=1
       if(present(iostat))return
    end if
    stop 
 end subroutine assert_longinteq
!-------------------------
  subroutine assert_singlerealeq(a,b,tol,msg,iostat)
    real(single),intent(in)::a,b
    real(single),optional,intent(in)::tol
    character(len=*),optional,intent(in)::msg
    integer(short),optional,intent(out)::iostat
    logical::inputOK=.false.
    stop
  end subroutine assert_singlerealeq
!-------------------------
  subroutine assert_doublerealeq(a,b,c,msg,iostat)
    real(double),intent(in)::a,b,c
    character(len=*),optional,intent(in)::msg
    integer(short),optional,intent(out)::iostat
    stop
  end subroutine assert_doublerealeq
!-------------------------
!  subroutine assert_quadrealeq(a,b,c,msg,iostat)
!    real(quad),intent(in)::a,b,c
!    character(len=*),intent(in)::msg
!    integer(short),optional,intent(out)::iostat
!  end subroutine assert_quadrealeq
!-------------------------
  subroutine assert_streq(a,b,msg,iostat)
    character(len=*),intent(in)::a,b
    character(len=*),optional,intent(in)::msg
    integer(short),optional,intent(out)::iostat
    logical::inputOK=.false.
    stop
  end subroutine assert_streq
!-------------------------


end module testing_class
