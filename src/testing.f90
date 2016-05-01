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
     module procedure check_doublecomplex
     module procedure check_doublecomplexpointer
     module procedure check_doublecomplexpointermatrix
     module procedure check_doublerealpointermatrix
  end interface check

contains
  integer(short) function check_doublerealpointermatrix(ptr)
    real(double),pointer,intent(in)::ptr(:,:)
    real(double)::targ(size(ptr,1),size(ptr,2))
    integer(short)::ierr(size(ptr,1),size(ptr,2))
    integer(long)::i,j
    !initiate with no problems found
    check_doublerealpointermatrix=0
    !check if pointer points to something
    if(.not.associated(ptr))then
       check_doublerealpointermatrix=-1
       return
    end if
    !check that elements are well behaved
    targ=ptr
    ierr=0
    do i=1,size(ptr,1)
       do j=1,size(ptr,2)
          ierr(i,j)=check(targ(i,j))
       end do
    end do
    if(any(ierr.NE.0))then
       check_doublerealpointermatrix=-1
    end if
    return
  end function check_doublerealpointermatrix
!-------------------------
  integer(short) function check_doublecomplexpointermatrix(ptr)
    complex(double),pointer,intent(in)::ptr(:,:)
    complex(double)::targ(size(ptr,1),size(ptr,2))
    integer(short)::ierr(size(ptr,1),size(ptr,2))
    integer(long)::i,j
    !initiate with no problems found
    check_doublecomplexpointermatrix=0
    !check if pointer points to something
    if(.not.associated(ptr))then
       check_doublecomplexpointermatrix=-1
       return
    end if
    !check that elements are well behaved
    targ=ptr
    ierr=0
    do i=1,size(ptr,1)
       do j=1,size(ptr,2)
          ierr(i,j)=check(targ(i,j))
       end do
    end do
    if(any(ierr.NE.0))then
       check_doublecomplexpointermatrix=-1
    end if
    return
  end function check_doublecomplexpointermatrix
!-------------------------
  integer(short) function check_doublecomplexpointer(ptr)
    complex(double),pointer,intent(in)::ptr(:)
    complex(double)::targ(size(ptr))
    integer(short)::ierr(size(ptr))
    integer(long)::i
    !initiate with no problems found
    check_doublecomplexpointer=0
    !check if pointer points to something
    if(.not.associated(ptr))then
       check_doublecomplexpointer=-1
       return
    end if
    !check that elements are well behaved
    targ=ptr
    ierr=0
    do i=1,size(ptr)
       ierr(i)=check(targ(i))
    end do
    if(any(ierr.NE.0))then
       check_doublecomplexpointer=-1
    end if
    return
  end function check_doublecomplexpointer
!-------------------------
  integer(short) function check_doublerealpointer(ptr)
    real(double),pointer,intent(in)::ptr(:)
    real(double)::targ(size(ptr))
    integer(short)::ierr(size(ptr))
    integer(long)::i
    !initiate with no problems found
    check_doublerealpointer=0
    !check if pointer points to something
    if(.not.associated(ptr))then
       check_doublerealpointer=-1
       return
    end if
    !check that elements are well behaved
    targ=ptr
    ierr=0
    do i=1,size(ptr)
       ierr(i)=check(targ(i))
    end do
    if(any(ierr.NE.0))then
       check_doublerealpointer=-1
    end if
    return
  end function check_doublerealpointer
!-------------------------
  integer(short) function check_doublecomplex(a)
    complex(double),intent(in)::a
    !initiate with no problems found
    check_doublecomplex=0
    !check real component
    check_doublecomplex=check(real(a))
    if(check_doublecomplex.NE.0)return
    !check imaginary component
    check_doublecomplex=check(aimag(a))
    if(check_doublecomplex.NE.0)return
    return
  end function check_doublecomplex
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
    logical::inputOK
    inputOK=.true.
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
    logical::inputOK
    inputOK=.true.
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
  subroutine assert_doublerealeq(a,b,tol,msg,iostat)
    real(double),intent(in)::a,b
    real(double),optional,intent(in)::tol
    character(len=*),optional,intent(in)::msg
    integer(short),optional,intent(out)::iostat
    logical::inputOK
    real(double)::tolerance
    inputOK=.true.
    if(present(iostat))iostat=-1
    if(check(a).NE.0)inputOK=.false.
    if(check(b).NE.0)inputOK=.false.
    tolerance=epsilon(a)
    if(present(tol))tolerance=tol
    if(abs(a-b).LT.tolerance.and.inputOK)then
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
