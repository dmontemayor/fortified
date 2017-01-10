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
  integer(long),parameter::prime(0:167)=(/2, 3, 5, 7, 11, 13, 17, 19, 23, 29&
       , 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103&
       , 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179&
       , 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257&
       , 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347&
       , 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431&
       , 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509&
       , 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607&
       , 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691&
       , 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797&
       , 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883&
       , 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991&
       , 997/)
  integer(long),parameter::twinprime(0:34)=(/3, 5, 11, 17, 29, 41, 59, 71&
       , 101, 107, 137, 149, 179, 191, 197, 227&
       , 239, 269, 281, 311, 347, 419, 431, 461&
       , 521, 569, 599, 617, 641, 659, 809, 821&
       , 827, 857, 881/)
       !,1019,1031,1049,1061,1091&
       !,1151,1229,1277,1289,1301,1319,1427,1451&
       !,1481,1487,1607,1619,1667,1697,1721,1787&
       !,1871,1877,1931,1949,1997,2027,2081,2087&
       !,2111,2129,2141,2237,2267,2309,2339,2381&
       !,2549,2591,2657,2687,2711,2729,2789,2801&
       !,2969,2999,3119,3167,3251,3257,3299,3329&
       !,3359,3371,3389,3461,3467,3527,3539,3557&
       !,3581,3671,3767,3821,3851,3917,3929,4001&
       !,4019,4049,4091,4127,4157,4217,4229,4241&
       !,4259,4271,4337,4421,4481,4517,4547,4637&
       !,4649,4721,4787,4799,4931,4967,5009,5021&
       !,5099,5231,5279,5417,5441,5477,5501,5519&
       !,5639,5651,5657,5741,5849,5867,5879,6089&
       !,6131,6197,6269,6299,6359,6449,6551,6569&
       !,6659,6689,6701,6761,6779,6791,6827,6869&
       !,6947,6959,7127,7211,7307,7331,7349,7457/)

  public::iden,cnorm,fnorm,krondelta,diag
  public::trace,imtrace
  public::numerov,solvestate,HOwf
  public::Hpoly

  !space filling curves
  public::spacefilling_index2coord,spacefilling_coord2index
  public::Hilbertcurve,Rowmajorcurve

  !statistics
  public::mean,variance,covariance,correlation,autocovariance,autocorrelation

  !sorting
  public::sort

  interface sort
     module procedure sort3D_real
  end interface sort
  
  interface trace
     module procedure trace_real
     module procedure trace_complex
  end interface trace
  
  interface Hilbertcurve
     module procedure hilbertcurve_real
     module procedure hilbertcurve_complex
  end interface Hilbertcurve

  interface Rowmajorcurve
     module procedure rowmajorcurve_real
     module procedure rowmajorcurve_complex
  end interface Rowmajorcurve

  interface krondelta
     module procedure krondelta_real
     module procedure krondelta_int
  end interface krondelta

  interface solvestate
     module procedure solve_numerov
  end interface solvestate

contains
  function sort3D_real(A)
    real(double),intent(in)::A(:,:,:)
    integer(long)::sort3D_real(size(A),3)
    logical::mask(size(A,1),size(A,2),size(A,3))
    integer(long)::curr,loc(3)
    sort3D_real=0
    mask=.true.
    curr=0
    do while(curr.LT.size(A))
       curr=curr+1
       loc=minloc(A,mask)
       sort3D_real(curr,:)=loc
       mask(loc(1),loc(2),loc(3))=.false.
    end do
  end function sort3D_real
  !-----------------------
  function mean(A)
    real(double),intent(in)::A(:)
    real(double)::mean,norm
    norm=real(size(A),double)
    mean=sum(A)/norm
  end function mean
  !----------------
  function variance(A)
    real(double),intent(in)::A(:)
    real(double)::variance,norm
    norm=real(size(A),double)
    variance=sum(A**2)/norm-mean(A)**2
  end function variance
  !--------------------
  function covariance(A,B)
    real(double),intent(in)::A(:),B(:)
    real(double)::covariance,comean,norm
    if(size(A).NE.size(B))then
       write(*,*)'covariance error: Array sizes are not the same size'
       stop
    end if
    norm=real(size(A),double)
    comean=mean(A)*mean(B)
    covariance=sum(A*B)/norm-comean
  end function covariance
  !----------------------
  function correlation(A,B)
    real(double),intent(in)::A(:),B(:)
    real(double)::correlation
    correlation=covariance(A,B)/sqrt(variance(A)*variance(B))
  end function correlation
  !-----------------------
  function autocovariance(A,inc)
    real(double),intent(in)::A(:)
    integer(long),intent(in)::inc
    real(double)::B(size(A)-inc)
    real(double)::autocovariance
    integer(long)::N,i
    N=size(A)-inc
    if(N.LT.0)then
       write(*,*)'autocovariance error: increment parameter cannot be greater than size of input array'
       stop
    end if
    if(N.EQ.0)then
       autocovariance=variance(A)
    else
       do i=1,N
          B(i)=A(i+inc)
       end do
       autocovariance=covariance(A(1:N),B(1:N))
    end if
  end function autocovariance
  !-----------------------
  function autocorrelation(A,inc)
    real(double),intent(in)::A(:)
    integer(long),intent(in)::inc
    real(double)::B(size(A)-inc)
    real(double)::autocorrelation
    integer(long)::N,i
    N=size(A)-inc
    if(N.LT.0)then
       write(*,*)'autocorrelation error: increment parameter cannot be greater than size of input array'
       stop
    end if
    if(N.EQ.0)then
       autocorrelation=1.0_double
    else
       do i=1,N
          B(i)=A(i+inc)
       end do
       autocorrelation=correlation(A(1:N),B(1:N))
    end if
  end function autocorrelation
  !---------------------------
  !> \brief Matrix coordinate from space filling curve coordinate
  !! \param[in] N leading dimension of NxM matrix
  !! \param[in] INDEX space filling curve index
  !! \param[in] TYPE optional string determining the type of 
  !!                 space filling curve. When missing or of unkown type
  !!                 a rowmajor curve type will be assumed. 
  !!<
  function spacefilling_index2coord(N,index,type)
    integer::spacefilling_index2coord(2) !matrix coordinate vector [0:N-1,0:M-1] output
    integer(long),intent(in)::N          !matrix leading dimension
    integer(long),intent(in)::index      !space filling curve index [0:N*M-1]
    character(len=*),optional,intent(in)::type
    character(len=label)::nameofcurve
    
    nameofcurve='rowmajor'
    if(present(type))nameofcurve=trim(type)
    
    select case(trim(adjustl(nameofcurve)))
    case('rowmajor')
       spacefilling_index2coord=rowmajor_index2coord(N,index)
    case('hilbert')
       spacefilling_index2coord=hilbert_index2coord(N,index)
    case default
       spacefilling_index2coord=rowmajor_index2coord(N,index)
    end select
    
  contains
    !- - - - - - - - - - - - - - - - - - - -
    !returns Matrix coordinates (i,j) from Row-Major curve index
    function Rowmajor_index2coord(N,index)
      integer(long)::Rowmajor_index2coord(2) !matrix coordinate vector      
      integer(long),intent(in)::N            !matrix leading dimension
      integer(long),intent(in)::index        !Row-Major curve index        
      
      rowmajor_index2coord(1)=index/N
      rowmajor_index2coord(2)=mod(index,N)
      
    end function Rowmajor_index2coord
    !- - - - - - - - - - - - - - - - - - - -
    !returns Square matrix coordinates (i,j) from Hilbert curve index
    function Hilbert_index2coord(N,index)
      integer(long)::Hilbert_index2coord(2)  !matrix coordinate vector      
      integer(long),intent(in)::N            !matrix leading dimension
      integer(long),intent(in)::index        !Hilbert curve index        
      integer(long)::s,rx,ry,x,y,z
      
      !initial value
      x=0
      y=0
      z=index
      s=1
      do while(s.LT.N)
         rx=iand(1,z/2)
         ry=iand(1,ieor(z,rx))
         call Hilbert_rotatequadrant(s,x,y,rx,ry)
         x=x+s*rx
         y=y+s*ry
         z=z/4
         s=s*2
      end do
      
      Hilbert_index2coord(1)=x
      Hilbert_index2coord(2)=y
      
    end function Hilbert_index2coord
    !- - - - - - - - - - - - - - - 
  end function spacefilling_index2coord
  !---------------------------
  !> \brief Space filling curve coordinate from matrix coordinate
  !! \param[in] N leading dimension of NxM matrix
  !! \param[in] COORD NxM Matrix coordinate vector 
  !! \param[in] TYPE optional string determining the type of 
  !!                 space filling curve. When missing or of unkown type
  !!                 a rowmajor curve type will be assumed. 
  !!<
  function spacefilling_coord2index(N,coord,type)
    integer::spacefilling_coord2index       !space filling curve index [0:N*M-1] output
    integer(long),intent(in)::N             !matrix leading dimension
    integer(long),intent(in)::coord(2)      !matrix coordinate vector [0:N-1,0:M-1]
    character(len=*),optional,intent(in)::type
    character(len=label)::nameofcurve

    nameofcurve='rowmajor'
    if(present(type))nameofcurve=trim(type)
    
    select case(trim(adjustl(nameofcurve)))
    case('rowmajor')
       spacefilling_coord2index=rowmajor_coord2index(N,coord)
    case('hilbert')
       spacefilling_coord2index=hilbert_coord2index(N,coord)
    case default
       spacefilling_coord2index=rowmajor_coord2index(N,coord)
    end select
    
  contains
    !- - - - - - - - - - - - - - - - - - - -
    !returns Row-Major curve index from Matrix coordinate vector
    function Rowmajor_coord2index(N,coord)
      integer(long)::Rowmajor_coord2index !Row-Major curve index
      integer(long),intent(in)::N          !matrix leading dimension
      integer(long),intent(in)::coord(2)  !matrix coordinate vector
      
      Rowmajor_coord2index=coord(1)*N+coord(2)
      
    end function Rowmajor_coord2index
    !- - - - - - - - - - - - - - - - - - - -
    !returns Hilbert curve index from Square Matrix coordinate vector
    function Hilbert_coord2index(N,coord)
      integer(long)::Hilbert_coord2index  !Hilbert curve index
      integer(long),intent(in)::N         !square matrix leading dimension
      integer(long),intent(in)::coord(2)  !square matrix coordinate vector
      integer(long)::s,rx,ry,x,y,z
      
      !initial value
      x=coord(1)
      y=coord(2)
      Hilbert_coord2index=0
      s=N
      !write(*,*)N*N,N,i,j
      do while(s.GT.0)
         s=s/2
         rx=iand(x,s)
         if(rx.GT.0)rx=1
         ry=iand(y,s)
         if(ry.GT.0)ry=1
         Hilbert_coord2index=Hilbert_coord2index+s*s*ieor(3*rx,ry)
         call Hilbert_rotatequadrant(N,x,y,rx,ry)
      end do
      !Hilbert_coord2index=Hilbert_coord2index+1
      !write(*,*)Hilbert_coord2index
      !if(j.EQ.10)stop
      
    end function Hilbert_coord2index
    !- - - - - - - - - - - - - - - 
  end function spacefilling_coord2index
  !- - - - - - - - - - - - - - - 
  subroutine Hilbert_rotatequadrant(N,x,y,rx,ry)
    integer(long),intent(in)::N,rx,ry
    integer(long),intent(inout)::x,y
    integer(long)::z
    
    if(ry.EQ.0)then
       if(rx.EQ.1)then
          x=n-1-x
          y=n-1-y
       end if
       !swap corrdinates
       z=x
       x=y
       y=z
    end if
    
  end subroutine Hilbert_rotatequadrant
!---------------------------
  !> \brief Real Rowmajor space filling curve from matrix
  !! \param[in] A real matrix of dimensions NxM
  !!<
  function Rowmajorcurve_real(A)
    real(double),intent(in)::A(:,:)
    real(double)::Rowmajorcurve_real(0:size(A,1)*size(A,2)-1)

    integer(long)::N,M,i,j,l,coord(2)

    !initial value
    Rowmajorcurve_real=0._double

    N=size(A,1)
    M=size(A,2)

    !check number of elements is not to large
    if(N*M.GT.huge(N))then
       write(*,*)'number of elements of input matrix A is too large'
       stop
    end if

    !assign values along rowmajor curve
    do i=0,N-1
       do j=0,N-1
          coord(1)=i
          coord(2)=j
          l=spacefilling_coord2index(N,coord,type='rowmajor')
          !write(*,*)i,j,l
          Rowmajorcurve_real(l)=A(i+1,j+1)
       end do
    end do
  end function Rowmajorcurve_real
!---------------------------
  !> \brief Complex Rowmajor space filling curve from 2D grid
  !! \param[in] A complex matrix of dimensions NxM
  !!<
  function Rowmajorcurve_complex(A)
    complex(double),intent(in)::A(:,:)
    real(double)::Rowmajorcurve_complex(0:size(A,1)*size(A,2)-1)

    integer(long)::N,M,i,j,l,coord(2)

    !initial value
    Rowmajorcurve_complex=0._double

    N=size(A,1)
    M=size(A,2)

    !check number of elements is not to large
    if(N*M.GT.huge(N))then
       write(*,*)'number of elements of input matrix A is too large'
       stop
    end if

    !assign values along rowmajor curve
    do i=0,N-1
       do j=0,N-1
          coord(1)=i
          coord(2)=j
          l=spacefilling_coord2index(N,coord,type='rowmajor')
          !write(*,*)i,j,l
          Rowmajorcurve_complex(l)=A(i+1,j+1)
       end do
    end do
  end function Rowmajorcurve_complex
  !-------------------------------------------------------------
  !> \brief Real Hilbert space filling curve from 2D grid
  !! \param[in] A real square matrix with leading dimension size power of 2
  !!<
  function Hilbertcurve_real(A)
    real(double),intent(in)::A(:,:)
    real(double)::Hilbertcurve_real(0:size(A,1)**2-1)

    integer(long)::N,M,i,j,l,coord(2)

    !initial value
    Hilbertcurve_real=0._double

    !check leading dimension is a power of 2
    N=size(A,1)
    M=1
    do while(M.LT.N)
       M=M*2
    end do
    if(M.NE.N)then
       write(*,*)'Hilbertcurve_real Error: input matrix leading&
            & dimension is not a power of 2'
       stop
    end if

    !check leading dimension is not to large
    if(N*N.GT.huge(N))then
       write(*,*)'leading dimension of input matrix A is too large'
       stop
    end if

    !check sure A is square    
    if(size(A).NE.N*N)then
       write(*,*)'Hilbertcurve_real Error: input matrix is not square'
       return
    end if

    !assign values along hilbert curve
    do i=0,N-1
       do j=0,N-1
          coord(1)=i
          coord(2)=j
          l=spacefilling_coord2index(N,coord,type='hilbert')
          !write(*,*)i,j,l
          Hilbertcurve_real(l)=A(i+1,j+1)
       end do
    end do
  end function Hilbertcurve_real
  !--------------------------------
  !> \brief Complex Hilbert space filling curve from 2D grid
  !! \param[in] A real square matrix with leading dimension size power of 2
  !!<
  function Hilbertcurve_complex(A)
    complex(double),intent(in)::A(:,:)
    real(double)::Hilbertcurve_complex(0:size(A,1)**2-1)

    integer(long)::N,M,i,j,l,coord(2)

    !initial value
    Hilbertcurve_complex=0._double

    !check leading dimension is a power of 2
    N=size(A,1)
    M=1
    do while(M.LT.N)
       M=M*2
    end do
    if(M.NE.N)then
       write(*,*)'Hilbertcurve_complex Error: input matrix leading&
            & dimension is not a power of 2'
       stop
    end if

    !check leading dimension is not to large
    if(N*N.GT.huge(N))then
       write(*,*)'leading dimension of input matrix A is too large'
       stop
    end if

    !check sure A is square    
    if(size(A).NE.N*N)then
       write(*,*)'Hilbertcurve_complex Error: input matrix is not square'
       return
    end if

    !assign values along hilbert curve
    do i=0,N-1
       do j=0,N-1
          coord(1)=i
          coord(2)=j
          l=spacefilling_coord2index(N,coord,type='hilbert')
          !write(*,*)i,j,l
          Hilbertcurve_complex(l)=A(i+1,j+1)
       end do
    end do
  end function Hilbertcurve_complex
  !--------------------------------

  !---------------------------
  !> \brief Harmonic Oscillator wavefunction 
  !! \param[in] N state number, 0 for ground state
  !! \param[in] MASS particle mass in atomic mass units
  !! \param[in] OMEGA angular frequency in atomic units hartree/hbar
  !! \param[in] X array defining values of x dof on a grid
  !! \param[inout] WF normalized wavefunction
  !! \param[inout] E (optional) state energy in hartree
  !!<
  subroutine HOwf(N,x,mass,omega,wf,E)
    integer(long),intent(in)::N
    real(double),intent(in)::x(:),mass,omega
    real(double),intent(inout)::wf(:)
    real(double),intent(inout),optional::E
    
    real(double)::sigma,norm
    integer(long)::i,nfactorial

    sigma=1.0_double/sqrt(mass*omega)
   ! nfactorial=1
   ! if(n.GT.1)then
   !    do i=1,n
   !       nfactorial=nfactorial*i
   !    end do
   ! end if
   ! norm=sqrt(2**n*nfactorial*sigma*sqrt(pi))
    wf=exp(-0.5_double*(x/sigma)**2)*Hpoly(n,x/sigma)
    norm=sqrt(sum(wf*wf))
    wf=wf/norm
   

 if(present(E))E=omega*(N+0.5_double)

  end subroutine HOwf
!---------------------------
  !> \brief Hermite polynomial
  !! \param[in] n nth polynomial to compute
  !! \param[in] x grid of x values
  !!<
  function Hpoly (n,x)
    integer(long),intent(in)::n
    real(double),intent(in)::x(:)
    real(double)::Hpoly(size(x)),Hprev(size(x)),Hcurr(size(x))
    integer(long)::i

    !initialize Hpoly for n=1 and Hprev for n=0
    Hprev=1.0_double
    Hpoly=2.0_double*x

    !Return if n 0 or 1
    if(n.EQ.0)then
       Hpoly=Hprev
       return
    end if
    if(n.EQ.1)return

    !Solve for n>1 recrusively
    do i=2,n
       Hcurr=Hpoly
       Hpoly=2.0_double*(x*Hcurr-(i-1)*Hprev)
       Hprev=Hcurr
    end do
    return
  end function Hpoly

!---------------------------
  !> \brief State energy for 1D potential
  !! \param[in] N state number, 0 for ground state
  !! \param[inout] E state energy in hartree output
  !! \param[inout] WF optional wavefunction output
  !! \param[in] V array containing potential energy in hartrees on equally spaced grid
  !! \param[in] MASS particle mass in atomic mass units
  !! \param[in] dx space discretization increment in bohr
  !!<
  subroutine solve_numerov(N,E,wf,V,mass,dx,dE,growth)
    integer(long),intent(in)::N
    real(double),intent(inout)::E
    real(double),intent(in)::V(:),mass,dx
    real(double),intent(inout),optional::wf(:)
    real(double),intent(in),optional::dE,growth

    integer(long)::i,npt,nnode,statesign
    integer(long)::maxattempt,nattempt
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
    maxattempt=1000
    growthfactor=1.05
    if(present(growth).and.growth.GT.0)growthfactor=1+growth

    !manually set lower energy bound to 0 relative to potential minimum E0
    Eb=0.0

    !set default initial upper energy bound to a tiny number
    E=epsilon(E)

    !Let user set initial upper energy bound
    if(present(dE))E=dE

    !Find upper energy bound
    !increase E until nnode>N
    nnode=-1
    nattempt=0
    do while(nnode.LE.N.and.nattempt.LT.maxattempt)
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
    if(nattempt.GE.maxattempt)then
       write(*,*)'solve_numerov error: too many attempts to find upper bound state'
       stop
    end if
    !set upper energy bound
    Et=E

    !adjust energy by divide and conquer until
    !change between last 2 points is tiny
    error=huge(dx)
    nattempt=0
    accuracy=epsilon(dx)
    !accuracy=dx**4
    Ep=Eb !set previous energy to lower energy bound
    E=Et  !set energy to upper bound to ensure energy gap
    !loop until wf stops diverging .and. too many attempts
    !    .and. energy stops changing
    do while(abs(error).GT.accuracy.and.nattempt.LT.maxattempt&
         .and.abs(E-Ep).GT.epsilon(E))
       nattempt=nattempt+1
       !save previous energy
       Ep=E
       !set energy between energy bounds
       E=(Eb+Et)/2.0_double
       !recalculate energy difference
       g=2.0_double*mass*(E+E0-V)

       !calculate numerov solution
       psi=numerov(dx,g,s)

       !calcualte error - asymptotic slope
       !error=psi(npt)
       error=abs(psi(npt)-psi(npt-1))

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
    if(nattempt.GE.maxattempt)then
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

    !!initialize numerov's solution to zero
    !numerov=0.0_double
    !initialize numerov's solution to accurracy
    numerov=accuracy

    !set first point very small
    !numerov(2)=epsilon(1.0_double)
    !numerov(2)=0.1_double*accuracy
    numerov(2)=2.0_double*accuracy

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
  !> \brief trace of imaginary part of a complex matrix
  !! \param[in] A matrix
  !!<
  function imtrace(A)
    complex(double),intent(in)::A(:,:)
    real(double)::imtrace
    integer(long)::i,N
    N=size(A,1)

    imtrace=0._double
    do i=1,N
       imtrace=imtrace+aimag(A(i,i))
    end do
  end function imtrace
  
!---------------------------
  !> \brief trace of a real matrix
  !! \param[in] A matrix
  !!<
  function trace_real(A)
    real(double),intent(in)::A(:,:)
    real(double)::trace_real
    integer(long)::i,N
    N=size(A,1)

    trace_real=0._double
    do i=1,N
       trace_real=trace_real+A(i,i)
    end do
  end function trace_real

!---------------------------
  !> \brief trace of a complex matrix
  !! \param[in] A matrix
  !!<
  function trace_complex(A)
    complex(double),intent(in)::A(:,:)
    real(double)::trace_complex
    integer(long)::i,N
    N=size(A,1)

    trace_complex=(0._double,0._double)
    do i=1,N
       trace_complex=trace_complex+sqrt(real(A(i,i)*conjg(A(i,i))))
    end do
  end function trace_complex
!--------------------------------
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
