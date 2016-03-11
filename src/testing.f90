module testing_class
  use type_kinds
  implicit none
  private

  public::assert,check

  interface assert
     module procedure assert_logical
     module procedure assert_shortinteq
     module procedure assert_longinteq
     module procedure assert_doublerealeq
     !module procedure assert_singlerealeq
     !module procedure assert_singlerealeqtol
     !module procedure assert_quadrealeq
     !module procedure assert_streq
     module procedure assert_filepresent
  end interface assert

  interface check
     module procedure check_shortint
     module procedure check_longint
     module procedure check_doublereal
     module procedure check_doublerealpointer
  end interface check


  
contains
  integer(short) function check_doublerealpointer(ptr)
    real(double),pointer,intent(in)::ptr(:)
    !initiate with no problems found
    check_doublerealpointer=0
    !check if pointer points to something
    if(.not.associated(ptr))then
       check_doublerealpointer=-1
       return
    end if
    !!check that elements are well behaved
    !iptr(:)=check(ptr(:))
    !forall (i=1:size(ptr)) iptr(i)=check(ptr(i))
    !if(any(iptr))
!!$    !check if real is huge
!!$    if(abs(a).GE.huge(a))then
!!$       check_doublerealpointer=-1
!!$       return
!!$    end if
!!$    !check if real is NaN
!!$    if(a.NE.a)then
!!$       check_doublerealpointer=-1
!!$       return
!!$    end if
    return
  end function check_doublerealpointer
!-------------------------
  integer(short) function check_doublereal(a)
    real(double),intent(in)::a
    !initiate with no problems found
    check_doublereal=0
    !check if real is huge
    if(abs(a).GE.huge(a))then
       check_doublereal=-1
       return
    end if
    !check if real is NaN
    if(a.NE.a)then
       check_doublereal=-1
       return
    end if
    return
  end function check_doublereal
!-------------------------
  integer(short) function check_shortint(a)
    integer(short),intent(in)::a
    !initiate with no problems found
    check_shortint=0
    !check if integer is huge
    if(abs(a).GE.huge(a))then
       check_shortint=-1
       return
    end if
    !check if integer is NaN
    if(a.NE.a)then
       check_shortint=-1
       return
    end if
    return
  end function check_shortint
!-------------------------
  integer(short) function check_longint(a)
    integer(long),intent(in)::a
    !initiate with no problems found
    check_longint=0
    !check if integer is huge
    if(abs(a).GE.huge(a))then
       check_longint=-1
       return
    end if
    !check if integer is NaN
    if(a.NE.a)then
       check_longint=-1
       return
    end if
    return
  end function check_longint
!-------------------------  
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
       if(present(msg))write(*,*)msg
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
    if(check(a).NE.0)inputOK=.false.
    if(check(b).NE.0)inputOK=.false.
    if(a.EQ.b.and.inputOK)then
       if(present(iostat))iostat=0
       return
    else
       if(present(iostat).and.inputOK)iostat=1
       if(present(msg))write(*,*)msg
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
    if(check(a).NE.0)inputOK=.false.
    if(check(b).NE.0)inputOK=.false.
    if(a.EQ.b.and.inputOK)then
       if(present(iostat))iostat=0
       return
    else
       if(present(msg))write(*,*)msg
       if(present(iostat).and.inputOK)iostat=1
       if(present(iostat))return
    end if
    stop 
 end subroutine assert_longinteq
!-------------------------
  subroutine assert_doublerealeq(a,b,msg,iostat)
    real(double),intent(in)::a,b
    character(len=*),optional,intent(in)::msg
    integer(short),optional,intent(out)::iostat
    logical::inputOK=.true.
    if(present(iostat))iostat=-1
    if(check(a).NE.0)inputOK=.false.
    if(check(b).NE.0)inputOK=.false.
    if(abs(a-b).LT.epsilon(a).and.inputOK)then
       if(present(iostat))iostat=0
       return
    else
       if(present(msg))write(*,*)msg
       if(present(iostat).and.inputOK)iostat=1
       if(present(iostat))return
    end if
    stop
  end subroutine assert_doublerealeq
!-------------------------
  subroutine assert_filepresent(file,msg,iostat)
    character(len=*),intent(in)::file
    character(len=*),optional,intent(in)::msg
    integer(short),optional,intent(out)::iostat
    logical::there
    
    if(present(iostat))iostat=-1
    inquire(file=file,exist=there)
    if(there)then
       if(present(iostat))iostat=0
       return
    else
       if(present(iostat))iostat=1
       if(present(msg))write(*,*)msg
       if(.not.present(iostat))stop
       return   
    end if
    stop

  end subroutine assert_filepresent


end module testing_class
