!===========================================================
!> \brief
!! Definition of mathematical constants and operators.
!! \authors
!! Author: Daniel Montemayor 2010
!! \todo
!! * Add cross and dot product operators
!<===========================================================
module math
  use type_kinds
  implicit none

  real(double),parameter::pi=3.1415926535898_double
  real(double),parameter::twopi=pi*2.0_double
  complex(double),parameter::eye=(0.0_double,1.0_double)
  public::iden,trace,cnorm,fnorm,krondelta,diag

  interface krondelta
     module procedure krondelta_real
     module procedure krondelta_int
  end interface krondelta

contains
!---------------------------
  !> \brief diagonal elements of a square matirx
  !! \param[in] A square real matrix
  !!<
  function diag(A)
    real(double),intent(in)::A(:,:)
    real(double)::diag(size(A,1))
    integer(long)::i
    diag=huge(diag)
    if(size(A,1).NE.size(A,2))then
       write(*,*)'Error Diag function: input matrix is not square!'       
       return
    end if
    do i=1,size(A,1)
       diag(i)=A(i,i)
    end do
  end function diag
!---------------------------
  !> \brief kronecker delta function
  !!<
  function krondelta_real(X,Y)
    real(double),intent(in)::X,Y
    real(double)::krondelta_real
    
    krondelta_real=0._double
    if(x.EQ.y)krondelta_real=1.0_double
  end function krondelta_real
!---------------------------
  !> \brief kronecker delta function
  !!<
  function krondelta_int(X,Y)
    integer(long),intent(in)::X,Y
    real(double)::krondelta_int
    
    krondelta_int=0._double
    if(x.EQ.y)krondelta_int=1.0_double
  end function krondelta_int
!---------------------------
  !> \brief identity matrix
  !! \param[in] N size of matrix
  !!<
  function iden(N)
    integer(long),intent(in)::N
    real(double)::iden(N,N)
    integer(long)::i

    iden=0._double
    do i=1,N
       iden(i,i)=1._double
    end do
  end function iden
!---------------------------
  !> \brief trace of a complex matrix
  !! \param[in] A matrix
  !!<
  function trace(A)
    complex(double),intent(in)::A(:,:)
    real(double)::trace
    integer(long)::i,N
    N=size(A,1)

    trace=0._double
    do i=1,N
       trace=trace+sqrt(real(A(i,i)*conjg(A(i,i))))
    end do
  end function trace

  !> \brief Entrywise p-norm of a complex square matrix with p=1
  !! \param[in] A matrix
  !!<
  function Cnorm(A)
    complex(double),intent(in)::A(:,:)
    real(double)::Cnorm
    integer(long)::i,j,N
    N=size(A,1)

    cnorm=0._double
    do i=1,N
       do j=1,N
          cnorm=cnorm+sqrt(real(A(i,j)*conjg(A(i,j))))
       end do
    end do
  end function Cnorm

  !> \brief Entrywise p-norm of a complex matrix
  !! \param[in] A matrix
  !! \param[in] p optional p-norm integer default: p=2 Frobinius normalization 
  !!<
  function Fnorm(A,p)
    complex(double),intent(in)::A(:,:)
    integer(long),intent(in),optional::p
    real(double)::Fnorm,s,sh,si
    integer(long)::i,j,N,M
    N=size(A,1)
    M=size(A,2)

    fnorm=0._double
    s=2._double
    if(present(p))then
       s=real(p,double)
       sh=s/2._double
       si=1._double/s
       do i=1,N
          do j=1,M
             fnorm=fnorm+real(A(i,j)*conjg(A(i,j)))**sh
          end do
       end do
       fnorm=fnorm**si
    else
       do i=1,N
          do j=1,M
             fnorm=fnorm+real(A(i,j)*conjg(A(i,j)))
          end do
       end do
       fnorm=sqrt(fnorm)
    end if
  end function Fnorm

end module math
!===========================================================
