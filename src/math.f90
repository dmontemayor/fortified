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
  public::numerov,solvestate

  interface krondelta
     module procedure krondelta_real
     module procedure krondelta_int
  end interface krondelta

  interface solvestate
module procedure solve_numerov
  end interface solvestate

contains
!---------------------------
  !> \brief State energy for 1D potential
  !! \param[in] N state number, 0 for ground state
  !! \param[inout] E state energy in hartree output
  !! \param[inout] WF optional wavefunction output
  !! \param[in] V array containing potential energy in hartrees on equally spaced grid
  !! \param[in] MASS particle mass in atomic mass units
  !! \param[in] dx space discretization increment in bohr
  !!<
  subroutine solve_numerov(N,E,wf,V,mass,dx)
    integer(long),intent(in)::N
    real(double),intent(inout)::E
    real(double),intent(in)::V(:),mass,dx
    real(double),intent(inout),optional::wf(:)

    integer(long)::i,npt,nnode,nattempt,statesign
    real(double),dimension(size(V))::psi,g,s
    real(double)::E0 !energy origin minval(V)
    real(double)::Eb !lower energy bound (bottom)
    real(double)::Et !upper energy bound (top)
    real(double)::Ep !previous energy to test convergence
    real(double)::error,accuracy,growthfactor

    !constants
    npt=size(V)
    s=0.0_double
    E0=minval(V)
    statesign=2*mod(N,2)-1
    accuracy=dx**4
    growthfactor=1.05

!!$    if(maxval(V).NE.minval(V))then
!!$       !set initial lower energy bound to 1% of potential energy range
!!$       E=1E-2*(maxval(V)-minval(V))
!!$    else
       !set initial lower energy bound to 5% of grnd state particle in a box
       !En=(n*hbar*pi/L)**2/(2m) n=1 for grnd state
       E=(pi/(npt*dx))**2/(2*mass)*.05
!!$    end if



    !Find lower energy bound
    !Assume nnode<=N; increase E until we have the correct number of nodes
    !Eb is penultimate energy setting before nnode=N
    nnode=-1
    nattempt=0
    Eb=E
    do while(nnode.LT.N.and.nattempt.LT.300)
       nattempt=nattempt+1
       !nnode<N so set lowerbound energy
       Eb=E
       !recalculate energy difference
       g=2.0_double*mass*(E+E0-V)
       !calculate numerov solution
       psi=numerov(dx,g,s)
       !count nodes
       nnode=0
       do i=2,npt
          if(psi(i)/psi(i-1).LT.0)nnode=nnode+1
       end do
       !write(*,*)E,nnode,nattempt
       !increase energy if necessary
       if(nnode.LT.N)E=E*growthfactor
    end do
    if(nattempt.GE.300)then
       write(*,*)'solve_numerov error: too many attempts to find lower bound state'
       stop
    end if

    !Find upper energy bound
    !increase E until nnode>N
    nnode=-1
    nattempt=0
    do while(nnode.LE.N.and.nattempt.LT.300)
       nattempt=nattempt+1
       !recalculate energy difference
       g=2.0_double*mass*(E+E0-V)
       !calculate numerov solution
       psi=numerov(dx,g,s)
       !count nodes
       nnode=0
       do i=2,npt
          if(psi(i)/psi(i-1).LT.0)nnode=nnode+1
       end do
       !write(*,*)E,nnode,nattempt
       !increase energy if neccessary
       if(nnode.LE.N)E=E*growthfactor
    end do
    if(nattempt.GE.300)then
       write(*,*)'solve_numerov error: too many attempts to find upper bound state'
       stop
    end if
    !set upper energy bound
    Et=E

    !adjust energy by divide and conquer until
    !magnitude of last point is tiny
    error=huge(dx)
    nattempt=0
    accuracy=epsilon(dx)
    Ep=Eb !set previous energy to lower energy bound
    E=Et  !set energy to upper bound to ensure energy gap
    !loop until wf stops diverging .and. too many attempts .and. energy stops changing
    do while(abs(error).GT.accuracy.and.nattempt.LT.300.and.abs(E-Ep).GT.epsilon(E))
       nattempt=nattempt+1
       !save previous energy
       Ep=E
       !set energy between energy bounds
       E=(Eb+Et)/2.0_double
       !recalculate energy difference
       g=2.0_double*mass*(E+E0-V)
       !calculate numerov solution
       psi=numerov(dx,g,s)
       error=psi(npt)
       !count nodes
       nnode=0
       do i=2,npt
          if(psi(i)/psi(i-1).LT.0)nnode=nnode+1
       end do
       !write(*,*)Eb,E,(Eb+Et)/2.0_double,Et,nnode,nattempt
       !energy too big if number of nodes is greater than N
       if(nnode.GT.N)Et=E
       !energy too small if number of nodes is less than N
       if(nnode.LT.N)Eb=E
       !if number of nodes equals N then
       if(nnode.EQ.N)then
          !energy too small if divergence is negative for even state
          !or positive divergence for odd state
          if(psi(npt)/statesign.LT.accuracy)Eb=E
          !energy too big if divergence is positive for even state
          !or negative divergence for odd state
          if(psi(npt)/statesign.GT.accuracy)Et=E
       end if
    end do
    if(nattempt.GE.300)then
       write(*,*)'solve_numerov error: too many attempts to find state'
       stop
    end if

    !add energy origin to final energy
    E=E+E0

    !normalize wavefunction
    psi=psi/sqrt(sum(psi**2)*dx)

    if(present(wf))then
       if(size(wf).EQ.npt)then
          wf=psi
       else
          write(*,*)'solve_numerov warning: wavefunction array provided wf is not same size as potential V&
               ; program will not change wf'
       end if
    end if


  end subroutine solve_numerov

!---------------------------
  !> \brief solution to diffeq of the form d^2/dx^2 y(x)=s(x)-g(x)y(x)
  !! \param[in] x grid of x values
  !! \param[in] s real function on grid
  !! \param[in] g real function on grid
  !!<
  function numerov (dx,g,s)
    real(double),intent(in)::dx,g(:),s(:)
    real(double)::numerov(size(g))
    real(double)::hh12,hh512,accuracy
    integer(long)::i,npt

    !parameters
    npt=size(g)
    accuracy=dx*dx*dx*dx
    hh12=dx*dx/12.0_double
    hh512=5.0_double*hh12

    !initialize numerov's solution to zero
    numerov=0.0_double

    !set first point very small
    !numerov(2)=epsilon(1.0_double)
    numerov(2)=0.1_double*accuracy!epsilon(1.0_double)

    !propagate solution
    do i=3,npt
       numerov(i)=(2.0_double*numerov(i-1)&
            *(1.0_double-hh512*g(i-1))&
            -numerov(i-2)*(1.0_double+hh12*g(i-2))&
            +hh12*(s(i)+10.0_double*s(i-1)+s(i-2)))&
            /(1.0_double+hh12*g(i))
    end do
    
    !!print solution
    !open(123,file='numerov.out')
    !do i=1,npt
    !   write(123,*)i*dx,numerov(i)
    !end do
    !close(123)


  end function numerov

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
