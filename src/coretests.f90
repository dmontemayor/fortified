subroutine coretests
  use type_kinds
  use testing_class
  use omp_lib
  use filemanager
  use math

  implicit none

  integer(short)::ierr
  logical::pass

  real(double)::x,y,S
  complex(double)::z
  real(double),pointer::ptr(:) 
  complex(double),pointer::zptr(:) 

  integer(long),parameter::order=7
  integer(long),parameter::npt=2**order
  real(double)::dx=0.05,Eng
  real(double),dimension(0:npt-1)::psi,grid,array1,array2,phi,chi
  real(double),dimension(0:npt-1,0:npt-1)::rho
  real(double),dimension(3,3)::ElCoupMat
  integer(long)::i,nnode,j,k

  !spacefilling curves
  integer(long)::N
  integer(long)::coord(2)
  real(double),dimension(0:(npt*npt)-1)::Lrho
  
  !Cuckier model constants
  real(double),parameter::hbar=1 !planck const
  real(double),parameter::me=1   !electron rest mass
  real(double),parameter::a0=1   !bhor radius
  real(double),parameter::Eh=1   !hartree
  real(double),parameter::invcm=4.5563E-6    !wavenumber in Eh
  real(double),parameter::angstrom=1.8897*a0 !angstrom in bhor radii
  !Cuckier model parameters
  real(double),parameter::gamma=1.0!ion-electron interaction decay distance
  real(double),parameter::mass=1836*me  !mobile ion mass (H+)
  real(double),parameter::Eb=2000*invcm !double well barrier height
  real(double),parameter::w0=1200*invcm !well frequency in harmonic approx 
  real(double),parameter::D=Eb/(hbar*w0)!Dimensionless barrier height 
  real(double),parameter::a=sqrt(hbar/(mass*w0))!spring cutoff distance



  !-----LOGICAL ASSERT TESTS-----

  write(*,*)'checking logical assert can be called with only a condition.....'
  call assert(.true.)

  write(*,*)'checking logical assert exits gracefully for false condition.....'
  call assert(.false.,iostat=ierr)

  write(*,*)'checking logical assert writes msg for false condition.....'
  write(*,*)'--------------------------------------------------------------'
  call assert(.false.,msg='logical assert is writting you a message.',iostat=ierr)
  write(*,*)'--------------------------------------------------------------'
  write(*,*)'Do you see a message in between the lines above? Answer true or false.'
  read(*,*)pass
  if(.not.pass)then
     write(*,*)'User cannot see a message passed by logical assert.'
     stop
  end if

  write(*,*)'checking logical assert does not write error message for true condition.....'
  call assert(.true.,msg='Error: logical assert writes error message for true condition.')
  
  write(*,*)'checking logical assert iostat returns 0 when cond is true.....' 
  call assert(.true.,iostat=ierr)
  if(ierr.NE.0)then
     write(*,*)'logical assert: iostat option does not return 0 for true condition.'
     stop
  end if

  write(*,*)'checking logical assert iostat returns 1 when cond is false.....' 
  call assert(.false.,iostat=ierr)
  if(ierr.NE.1)then
     write(*,*)'logical assert: iostat option does not return 1 for false condition'
     stop
  end if

  !-----SHORT INT ASSERT TESTS-----
  write(*,*)'checking short int assert can be called with only two short integers.....'
  call assert(1_short,1_short)

  write(*,*)'checking short int assert exits gracefully for not equal condition.....'
  call assert(1_short,2_short,iostat=ierr)

  write(*,*)'checking short int assert writes msg for not equal condition.....'
  write(*,*)'--------------------------------------------------------------'
  call assert(1_short,2_short,msg='short int assert is writting you a message.',iostat=ierr)
  write(*,*)'--------------------------------------------------------------'
  write(*,*)'Do you see a message in between the lines above? Answer true or false.'
  read(*,*)pass
  if(.not.pass)then
     write(*,*)'User cannot see a message passed by short int assert.'
     stop
  end if

  write(*,*)'checking short int assert does not write error message for equal condition.....'
  call assert(1_short,1_short,msg='Error: short int assert writes error message for equal condition.')
  
  write(*,*)'checking short int assert iostat returns 0 for equal condition.....' 
  call assert(1_short,1_short,iostat=ierr)
  if(ierr.NE.0)then
     write(*,*)'short int assert: iostat option does not return 0 for equal condition.'
     stop
  end if

  write(*,*)'checking short int assert iostat returns 1 for not equal condition.....' 
  call assert(1_short,2_short,iostat=ierr)
  if(ierr.NE.1)then
     write(*,*)'short int assert: iostat option does not return 1 for not equal condition'
     stop
  end if

  write(*,*)'checking short int assert returns -1 for huge numbers.....' 
  call assert(Huge(1_short),Huge(1_short),iostat=ierr)
  if(ierr.NE.-1)then
     write(*,*)'short int assert: iostat option does not return -1 for huge numbers'
     stop
  end if

  !-----LONG INT ASSERT TESTS-----
  write(*,*)'checking long int assert can be called with only two long integers.....'
  call assert(1_long,1_long)

  write(*,*)'checking long int assert exits gracefully for not equal condition.....'
  call assert(1_long,2_long,iostat=ierr)

  write(*,*)'checking long int assert writes msg for not equal condition.....'
  write(*,*)'--------------------------------------------------------------'
  call assert(1_long,2_long,msg='long int assert is writting you a message.',iostat=ierr)
  write(*,*)'--------------------------------------------------------------'
  write(*,*)'Do you see a message in between the lines above? Answer true or false.'
  read(*,*)pass
  if(.not.pass)then
     write(*,*)'User cannot see a message passed by long int assert.'
     stop
  end if

  write(*,*)'checking long int assert does not write error message for equal condition.....'
  call assert(1_long,1_long,msg='Error: long int assert writes error message for equal condition.')
  
  write(*,*)'checking long int assert iostat returns 0 for equal condition.....' 
  call assert(1_long,1_long,iostat=ierr)
  if(ierr.NE.0)then
     write(*,*)'long int assert: iostat option does not return 0 for equal condition.'
     stop
  end if

  write(*,*)'checking long int assert iostat returns 1 for not equal condition.....' 
  call assert(1_long,2_long,iostat=ierr)
  if(ierr.NE.1)then
     write(*,*)'long int assert: iostat option does not return 1 for not equal condition'
     stop
  end if

  write(*,*)'checking long int assert returns -1 for huge numbers.....' 
  call assert(Huge(1_long),Huge(1_long),iostat=ierr)
  if(ierr.NE.-1)then
     write(*,*)'long int assert: iostat option does not return -1 for huge numbers'
     stop
  end if

  !-----DOUBLE REAL ASSERT TESTS-----
  write(*,*)'checking double real assert can be called with only two double reals.....'
  call assert(1.0_double,1.0_double)

  write(*,*)'checking double real assert exits gracefully for not equal condition.....'
  call assert(1.0_double,2.0_double,iostat=ierr)

  write(*,*)'checking double real assert writes msg for not equal condition.....'
  write(*,*)'--------------------------------------------------------------'
  call assert(1.0_double,2.0_double,msg='double real assert is writting you a message.',iostat=ierr)
  write(*,*)'--------------------------------------------------------------'
  write(*,*)'Do you see a message in between the lines above? Answer true or false.'
  read(*,*)pass
  if(.not.pass)then
     write(*,*)'User cannot see a message passed by double real assert.'
     stop
  end if

  write(*,*)'checking double real assert does not write error message for equal condition.....'
  call assert(1.0_double,1.0_double,msg='Error: double real assert writes error message for equal condition.')
  
  write(*,*)'checking double real assert iostat returns 0 for equal condition.....' 
  call assert(1.0_double,1.0_double,iostat=ierr)
  if(ierr.NE.0)then
     write(*,*)'double real assert: iostat option does not return 0 for equal condition.'
     stop
  end if

  write(*,*)'checking double real assert iostat returns 1 for not equal condition.....' 
  call assert(1.0_double,2.0_double,iostat=ierr)
  if(ierr.NE.1)then
     write(*,*)'double real assert: iostat option does not return 1 for not equal condition'
     stop
  end if

  write(*,*)'checking double real assert returns -1 for huge numbers.....' 
  call assert(Huge(1.0_double),Huge(1.0_double),iostat=ierr)
  if(ierr.NE.-1)then
     write(*,*)'double real assert: iostat option does not return -1 for huge numbers'
     stop
  end if

  write(*,*)'checking double real assert can be called with tolerance option.....'
  call assert(1.0_double,1.0_double,tol=0.1_double,iostat=ierr) 

  write(*,*)'checking double real assert returns 0 for two vaules within tolerance.....'
  call assert(1.0_double,2.0_double,tol=3.0_double,iostat=ierr)
  if(ierr.NE.0)then
     write(*,*)'double real assert: iostat option does not return 0 for two values within tolerance'
     stop
  end if

  write(*,*)'checking double real assert returns 1 for two vaules not within tolerance.....'
  call assert(1.0_double,2.0_double,tol=0.1_double,iostat=ierr)
  if(ierr.NE.1)then
     write(*,*)'double real assert: iostat option does not return 1 for two values not within tolerance'
     stop
  end if

  !-----FILE ASSERT
  write(*,*)'checking file assert returns 1 for missing file.....'
  call system('rm -f coretest.tmpfile')
  call assert('coretest.tmpfile',iostat=ierr)
  if(ierr.NE.1)then
     write(*,*)'file assert: iostat option does not return 1 for missing file'
     stop
  end if

  write(*,*)'checking file assert returns 0 for present file.....'
  call system('touch coretest.tmpfile')
  call assert('coretest.tmpfile',msg='file check did not return 0 for present file')
  call system('rm -f coretest.tmpfile')


!!$  !--------- OPEN MP ----------
!!$  write(*,*)'checking omp is available.....'
!!$  call assert(omp_get_num_procs().GT.0,msg=&
!!$       'number of available processors is not greater than 0')
!!$  call assert(omp_get_max_threads().GT.0,msg=&
!!$       'number of available threads is not greater than 0')
!!$  call assert(omp_get_max_threads().GE.omp_get_num_procs(),msg=&
!!$       'number of threads is less than number of processors available'&
!!$       ,iostat=ierr)
!!$  if(ierr.NE.0)then
!!$     write(*,*)'try setting the environment variable with'
!!$     write(*,*)'export OMP_NUM_THREADS=',omp_get_num_procs()
!!$  end if


  !--- filemanager ------
  write(*,*)'checking filemanager assigns new units starting at 1000.....'
  call assert(newunit().GE.1000,msg='newunit returns unit below 1000')

  write(*,*)'checking filemanager can check for missing file.....'
  call system('rm -f coretest.tmpfile')
  call assert(check('coretest.tmpfile').EQ.1,msg='file check did not return 1 for missing file')

  write(*,*)'checking filemanager can check for present file.....'
  call system('touch coretest.tmpfile')
  call assert(check('coretest.tmpfile').EQ.0,msg='file check did not return 0 for present file')
  call system('rm -f coretest.tmpfile')

  write(*,*)'checking filemanager can check for present file of known type.....'
  call system('touch coretest.tmpfile')
  call system("echo 'knowntype' >>coretest.tmpfile")
  call assert(check('coretest.tmpfile','knowntype').EQ.0,msg='file check did not return 0 for present file of known type')
  call system('rm -f coretest.tmpfile')

  write(*,*)'checking filemanager returns -1 for present file of wrong type.....'
  call system('touch coretest.tmpfile')
  call system("echo 'knowntype' >>coretest.tmpfile")
  call assert(check('coretest.tmpfile','unknowntype').EQ.-1,msg='file check did not return -1 for present file of wrong type')
  call system('rm -f coretest.tmpfile')

  
  !-----SHORT INT CHECK TESTS-----
  write(*,*)'test check short int returns 0 for the integer 1.....' 
  if(check(1_short).NE.0)then
     write(*,*)'check short int does not return 0 for the integer 1'
     stop
  end if
  write(*,*)'test check short int returns -1 for huge numbers.....' 
  if(check(Huge(1_short)).NE.-1)then
     write(*,*)'check short int does not return -1 for huge numbers'
     stop
  end if
  write(*,*)'test check short int returns -1 for negative huge numbers.....' 
  if(check(-Huge(1_short)).NE.-1)then
     write(*,*)'check short int does not return -1 for negative huge numbers'
     stop
  end if
  !-----LONG INT CHECK TESTS-----
  write(*,*)'test check long int returns 0 for the integer 1.....' 
  if(check(1_long).NE.0)then
     write(*,*)'check long int does not return 0 for the integer 1'
     stop
  end if
  write(*,*)'test check long int returns -1 for huge numbers.....' 
  if(check(Huge(1_long)).NE.-1)then
     write(*,*)'check long int does not return -1 for huge numbers'
     stop
  end if
  write(*,*)'test check long int returns -1 for negative huge numbers.....' 
  if(check(-Huge(1_long)).NE.-1)then
     write(*,*)'check long int does not return -1 for negative huge numbers'
     stop
  end if

  !-----DOUBLE REAL CHECK TESTS-----
  write(*,*)'test check double real returns 0 for the value 1.0.....' 
  if(check(1.0_double).NE.0)then
     write(*,*)'check double real does not return 0 for the value 1.0'
     stop
  end if
  write(*,*)'test check double real returns -1 for huge numbers.....' 
  if(check(Huge(1.0_double)).NE.-1)then
     write(*,*)'check double real does not return -1 for huge numbers'
     stop
  end if
  write(*,*)'test check double real returns -1 for negative huge numbers.....' 
  if(check(-Huge(1.0_double)).NE.-1)then
     write(*,*)'check double real does not return -1 for negative huge numbers'
     stop
  end if
  write(*,*)'test check double real returns -1 for infinity.....' 
  x=0._double
  x=1/x
  if(check(x).NE.-1)then
     write(*,*)'check double real does not return -1 for infinity'
     stop
  end if
 
  !-----DOUBLE REAL POINTER CHECK TESTS-----
  write(*,*)'test check double real pointer returns 0 for associated pointer.....'
  if(associated(ptr))nullify(ptr)
  allocate(ptr(0:1))
  if(check(ptr).NE.0)then
     write(*,*)'check double real pointer does not return 0 for associated pointer'
     stop
  end if
  write(*,*)'test check double real pointer returns -1 for unassociated pointer.....'
  if(associated(ptr))nullify(ptr)
  if(check(ptr).NE.-1)then
     write(*,*)'check double real pointer does not return -1 for unassociated pointer'
     stop
  end if
  write(*,*)'test check double real pointer returns -1 for poor behaved elements.....'
  if(associated(ptr))nullify(ptr)
  allocate(ptr(0:1))
  x=0._double
  ptr(:)=1/x
  if(check(ptr).NE.-1)then
     write(*,*)'check double real pointer does not return -1 for poor behaved elements'
     stop
  end if
  
  !-----DOUBLE COMPLEX POINTER CHECK TESTS-----
  write(*,*)'test check double complex pointer returns 0 for associated pointer.....'
  if(associated(zptr))nullify(zptr)
  allocate(zptr(0:1))
  if(check(zptr).NE.0)then
     write(*,*)'check double complex pointer does not return 0 for associated pointer'
     stop
  end if
  write(*,*)'test check double complex pointer returns -1 for unassociated pointer.....'
  if(associated(zptr))nullify(zptr)
  if(check(zptr).NE.-1)then
     write(*,*)'check double complex pointer does not return -1 for unassociated pointer'
     stop
  end if
  write(*,*)'test check double complex pointer returns -1 for poor behaved elements.....'
  if(associated(zptr))nullify(zptr)
  allocate(zptr(0:1))
  z=(0._double,0._double)
  zptr(:)=(0.,1.)/z
  if(check(zptr).NE.-1)then
     write(*,*)'check double complex pointer does not return -1 for poor behaved elements'
     stop
  end if
  
  !check numerov's method
  !--------------------------------
  !harmonic oscillator constants
  !hbar=1
  !omega=1
  !mass=1
  !En=hbar*omega(n+1/2)
  !V(x)=1/2*mass*omega^2*x^2
  !--------------------------------
  write(*,*)'test numerov solution to ground state harmonic oscillator obeys proper boundary conditions.....'
  !set xgrid and calculated array1 and array2
  do i=0,npt-1
     grid(i)=(i-npt/2)*dx !xgrid 
  end do
  array1=1.0_double-grid*grid !energy diferrence =2*mass(E-V(x))/hbar=2E-2V
  array2=0.0_double           !independent function =0
  psi=numerov(dx,array1,array2)
  !if(abs(psi(0)).GT.epsilon(1.0_double))then
  !   write(*,*)'numerov solution does not equal zero for first point'
  !   stop
  !end if
  x=abs(psi(1)-psi(0))
  if(x.GT.dx**4)then
     write(*,*)'point diff=',x,'accuracy=',dx**4
     write(*,*)'numerov solution starting point difference is greater than dx**4'
     stop
  end if
  x=abs(psi(npt-1)-psi(npt-2))
  if(x.GT.dx**4)then
     write(*,*)'point diff=',x,'accuracy=',dx**4
     write(*,*)'numerov solution final point difference is greater than dx**4'
     stop
  end if

  write(*,*)'test numerov solution to first excited state harmonic oscillator has one node.....'
  !set xgrid and calculated array1 and array2
  do i=0,npt-1
     grid(i)=(i-npt/2)*dx !xgrid 
  end do
  array1=3.0_double-grid*grid !energy diferrence =2*mass(E-V(x))/hbar=2E-2V
  array2=0.0_double           !independent function =0
  psi=numerov(dx,array1,array2)
  !count number of nodes
  nnode=0
  do i=2,npt-1
     if(psi(i)/psi(i-1).LT.0)nnode=nnode+1
  end do
  if(nnode.NE.1)then
     write(*,*)'numerov solution does not have only one node'
     stop
  end if

  !test solvestate
  write(*,*)'test solvestate returns particle in a box ground state energy.....'
  !--------------------------------
  !Particle in a bax constants
  !hbar=1
  !mass=1
  !L=npt*dx=128*.05=6.4
  !En=(n*hbar*pi/L)**2/(2m) n=1 for grnd state
  !V(x)=0=array1
  !--------------------------------
  !set xgrid and calculated array1 and array2
  do i=0,npt-1
     grid(i)=(i-npt/2)*dx !xgrid 
  end do
  array1=0._double !potential energy
  call solvestate(N=0,E=Eng,V=array1,mass=1.0_double,dx=dx)
  x=0.5_double*(pi/(npt*dx))**2
  x=(Eng-x)/x !relative error
  if(abs(x).GT.2E-2)then
     write(*,*)'result=',Eng,'Theo. solution=',x,'Rel. Error=',abs((Eng-x)/x)
     write(*,*)'solvestate energy relative error is greater than .02 for grnd state particle in a box'
     stop
  end if

  write(*,*)'test solvestate returns normalized particle in a box 1st excited state wave function.....'
  !--------------------------------
  !Particle in a bax constants
  !hbar=1
  !mass=1
  !L=npt*dx=128*.05=6.4
  !En=(n*hbar*pi/L)**2/(2m) n=1 for grnd state
  !V(x)=0=array1
  !--------------------------------
  !set xgrid and calculated array1 and array2
  do i=0,npt-1
     grid(i)=(i-npt/2)*dx !xgrid 
  end do
  array1=0._double !potential energy
  call solvestate(N=1,E=Eng,wf=psi,V=array1,mass=1.0_double,dx=dx)
  call solvestate(N=1,E=Eng,wf=phi,V=array1,mass=1.0_double,dx=dx)
  x=0.0_double
  do i=0,npt-1
     x=x+psi(i)*phi(i)
  end do
  x=x*dx
  !write(*,*)x
  if(abs(x-1.0_double).GT.2E-2)then
     write(*,*)'solvestate normalization is greater than .02 for 1th excited state particle in a box'
     stop
  end if

  write(*,*)'test solvestate returns orthogonal particle in a box 1st and 5th excited state wavefunctions.....'
  !--------------------------------
  !Particle in a bax constants
  !hbar=1
  !mass=1
  !L=npt*dx=128*.05=6.4
  !En=(n*hbar*pi/L)**2/(2m) n=1 for grnd state
  !V(x)=0=array1
  !--------------------------------
  !set xgrid and calculated array1 and array2
  do i=0,npt-1
     grid(i)=(i-npt/2)*dx !xgrid 
  end do
  array1=0._double !potential energy
  call solvestate(N=1,E=Eng,wf=psi,V=array1,mass=1.0_double,dx=dx)
  call solvestate(N=5,E=Eng,wf=phi,V=array1,mass=1.0_double,dx=dx)
  x=0.0_double
  do i=0,npt-1
     x=x+psi(i)*phi(i)
  end do
  x=x*dx
  !write(*,*)x
  if(abs(x).GT.2E-2)then
     write(*,*)'solvestate orthogonalization for particle in a box 1st and 5th&
          excited state wavefunctions is greater than 2%'
     stop
  end if

  write(*,*)'test solvestate returns Cukier proton sates.....'
  !--------------------------------
  !Cukier potential
  !hbar=1
  !me=1
  !a0=1
  !invcm=4.5563E-6
  !angstrom=1.8897*a0
  !mass=1836*me
  !Eb=2000*invcm
  !w0=1200*invcm
  !D=Eb/(hbar*w0)
  !a=sqrt(hbar/(mass*w0))
  !L=npt*dx=128*.018=2.304
  !gamma=1.0_double
  !--------------------------------
  dx=0.018*angstrom
  !set xgrid and calculated array1(proton potential) and array2
  do i=0,npt-1
     grid(i)=(i-npt/2)*dx !xgrid
     x=grid(i)/a
     y=x*x/4
     array1(i)=y*(1-y/(4*D))   !potential energy
  end do
  array1=Eb-(hbar*w0)*array1

  !calculate inverted potential
  do i=0,npt-1
     array2(npt-1-i)=array1(i) !inverted potential energy
  end do
  !calculate first excited states
  call solvestate(N=1,E=x,wf=psi,V=array1,mass=mass,dx=dx)
  call solvestate(N=1,E=y,wf=phi,V=array2,mass=mass,dx=dx)
  !normalize wavefunctions
  psi=psi/sqrt(sum(psi*psi))
  phi=-phi/sqrt(sum(phi*phi))

  !write(*,*)'Diagnostic: write Cukier potential with first and&
  !& inverted state solutions.'
  !open(123,file='Cukier-invertedstate.dat')
  !do i=0,npt-1
  !   write(123,*)grid(i),array1(i)&
  !        ,psi(i)*(hbar*w0),phi(i)*(hbar*w0)
  !end do
  !close(123)
  !write(*,*)'solvestate returns same first excited wf when grid is&
  !& inverted for Cukier potential'

  !calculate overlap
  S=0.0_double
  do i=0,npt-1
     S=S+psi(i)*phi(npt-1-i)
  end do
  S=sqrt(S)
  !write(*,*)S
  call assert(S.GT.0.98,msg='first excited wf not same (<2% error) as with inverted grid.')


  !!diagnostic write potential and first 3 state solutions
  !call solvestate(N=0,E=x,wf=psi,V=array1,mass=mass,dx=dx)
  !call solvestate(N=1,E=y,wf=phi,V=array1,mass=mass,dx=dx)
  !call solvestate(N=2,E=Eng,wf=chi,V=array1,mass=mass,dx=dx)
  !open(123,file='CukierModel.dat')
  !do i=0,npt-1
  !   write(123,*)grid(i)/angstrom,array1(i)&
  !        ,psi(i)*(hbar*w0)+x,phi(i)*(hbar*w0)+y&
  !        ,chi(i)*(hbar*w0)+Eng,.5*(psi(i)*(hbar*w0)+phi(i)*(hbar*w0))
  !end do
  !close(123)

  !!diagnostic calculate electronic coupling matrix
  !call solvestate(N=0,E=x,wf=phi,V=array1,mass=mass,dx=dx)
  !call solvestate(N=1,E=y,wf=chi,V=array1,mass=mass,dx=dx)
  !psi=.5*(phi+chi)
  !phi=.5*(phi-chi)
  !call solvestate(N=2,E=Eng,wf=chi,V=array1,mass=mass,dx=dx)
  !ElCoupMat(1,1)=sum(psi*exp(-gamma*abs(grid))*psi)
  !ElCoupMat(1,2)=sum(psi*exp(-gamma*abs(grid))*phi)
  !ElCoupMat(1,3)=sum(psi*exp(-gamma*abs(grid))*chi)
  !ElCoupMat(2,1)= ElCoupMat(1,2)
  !ElCoupMat(2,2)=sum(phi*exp(-gamma*abs(grid))*phi)
  !ElCoupMat(2,3)=sum(phi*exp(-gamma*abs(grid))*chi)
  !ElCoupMat(3,1)= ElCoupMat(1,3)
  !ElCoupMat(3,2)= ElCoupMat(2,3)
  !ElCoupMat(3,3)=sum(chi*exp(-gamma*abs(grid))*chi)
  !open(123,file='CukierCouplingMat.dat')
  !do i=1,3
  !   write(123,*)(ElCoupMat(i,j),j=1,3)
  !   write(*,*)(ElCoupMat(i,j),j=1,3)
  !end do
  !close(123)

  write(*,*)'test solvestate returns particle in a box 1st excited state energy with error<2% .....'
  !--------------------------------
  !Particle in a bax constants
  !hbar=1
  !mass=1
  !L=npt*dx=128*.05=6.4
  !En=(n*hbar*pi/L)**2/(2m) n=2 for first excited state
  !V(x)=0=array1
  !--------------------------------
  !set xgrid and calculated array1 and array2
  do i=0,npt-1
     grid(i)=(i-npt/2)*dx !xgrid 
  end do
  array1=0._double !potential energy
  call solvestate(N=1,E=Eng,V=array1,mass=1.0_double,dx=dx)
  x=0.5_double*(2*pi/(npt*dx))**2
  !write(*,*)Eng,x,abs((Eng-x)/x)
  x=(Eng-x)/x !relative error
  if(abs(x).GT.2E-2)then
     write(*,*)'solvestate energy relative error is greater than .02&
          &for 1t exited state particle in a box'
     stop
  end if

  write(*,*)'test solvestate can be called with energy increment option .....'
  !--------------------------------
  !Particle in a bax constants
  !hbar=1
  !mass=1
  !L=npt*dx=128*.05=6.4
  !En=(n*hbar*pi/L)**2/(2m) n=2 for first excited state
  !V(x)=0=array1
  !--------------------------------
  !set xgrid and calculated array1 and array2
  do i=0,npt-1
     grid(i)=(i-npt/2)*dx !xgrid 
  end do
  array1=0._double !potential energy
  x=0.1
  call solvestate(N=1,E=Eng,V=array1,mass=1.0_double,dx=dx,dE=x)

  write(*,*)'test solvestate can be called with growth factor option .....'
  !--------------------------------
  !Particle in a bax constants
  !hbar=1
  !mass=1
  !L=npt*dx=128*.05=6.4
  !En=(n*hbar*pi/L)**2/(2m) n=2 for first excited state
  !V(x)=0=array1
  !--------------------------------
  !set xgrid and calculated array1 and array2
  do i=0,npt-1
     grid(i)=(i-npt/2)*dx !xgrid 
  end do
  array1=0._double !potential energy
  x=0.1
  call solvestate(N=1,E=Eng,V=array1,mass=1.0_double,dx=dx,growth=x)

  write(*,*)'test solvestate returns correct ground state harmonic oscillator energy.....'
  !--------------------------------
  !harmonic oscillator constants
  !hbar=1
  !omega=1
  !mass=1
  !L=npt*dx=128*.05=6.4
  !En=hbar*omega(n+1/2)
  !V(x)=1/2*mass*omega^2*x^2
  !--------------------------------
  !set xgrid and calculated array1 and array2
  dx=0.05_double
  do i=0,npt-1
     grid(i)=(i-npt/2)*dx !xgrid 
  end do
  array1=0.5_double*grid**2 !potential energy
  call solvestate(N=0,E=Eng,wf=psi,V=array1,mass=1.0_double,dx=dx)
  x=0.5_double
  x=(Eng-x)/x !relative error
  if(abs(x).GT.2E-2)then
     write(*,*)'Calc. Eng=',Eng,'Theo. Eng=',x,'Rel. error=',abs((Eng-x)/x)
     write(*,*)'solvestate energy relative error is greater than .02&
          & for ground state harmonic oscillator'

     write(*,*)'Diagnostic: write hamonic oscillator potential with&
     & calculated ground state.'
     open(123,file='solvestate-HOgroundstate.dat')
     do i=0,npt-1
        write(123,*)grid(i),array1(i)&
             ,psi(i)+Eng
     end do
     close(123)

     stop
  end if

  write(*,*)'test solvestate returns correct ground state energy for shifted harmonic oscillator.....'
  !--------------------------------
  !harmonic oscillator constants
  !hbar=1
  !omega=1
  !mass=1
  !L=npt*dx=128*.05=6.4
  !En=hbar*omega(n+1/2)
  !V(x)=1/2*mass*omega^2*x^2-1
  !--------------------------------
  !set xgrid and calculated array1 and array2
  dx=0.05_double
  do i=0,npt-1
     grid(i)=(i-npt/2)*dx !xgrid 
  end do
  array1=0.5_double*grid**2-1.0_double !potential energy
  call solvestate(N=0,E=Eng,wf=psi,V=array1,mass=1.0_double,dx=dx)
  x=-0.5_double
  x=(Eng-x)/x !relative error
  if(abs(x).GT.2E-2)then
     write(*,*)'Calc. Eng=',Eng,'Theo. Eng=',x,'Rel. error=',abs((Eng-x)/x)
     write(*,*)'solvestate energy relative error is greater than .02&
          & for ground state harmonic oscillator'

     write(*,*)'Diagnostic: write shifted hamonic oscillator potential with&
     & calculated ground state.'
     open(123,file='solvestate-shiftedHOgroundstate.dat')
     do i=0,npt-1
        write(123,*)grid(i),array1(i)&
             ,psi(i)+Eng
     end do
     close(123)

     stop
  end if



  write(*,*)'test Hpoly returns correct hermite polynomial defintion.....'
  !set xgrid and calculated array1 and array2
  dx=0.05_double
  do i=0,npt-1
     grid(i)=(i-npt/2)*dx !xgrid 
  end do
  ! Calculate 4th Hermite polynomial H_4(x)=16x^4-48x^2+12
  array1=Hpoly(4,grid)
  array2=16*grid**4-48*grid**2+12
  !calcualte rmsd x
  x=sqrt(sum((array1-array2)**2)/real(npt))
  if(x.GT.2E-2)then
     write(*,*)'Hpoly(n=4) RMSD=',x
     write(*,*)'RMSD error is greater than .02 on xrange [-3.4:3.4]&
          & with discritization of 128 grid points'
     open(123,file='Hpolytest.dat')
     do i=0,npt-1
        write(123,*)grid(i),array1(i),array2(i)
     end do
     close(123)
     stop
  end if

  write(*,*)'test HOwf returns correct harmonic oscillator wavefunction.....'
  !--------------------------------
  !harmonic oscillator constants
  !hbar=1
  !omega=1
  !mass=1
  !L=npt*dx=128*.05=6.4
  !En=hbar*omega(n+1/2)
  !V(x)=1/2*mass*omega^2*x^2-1
  !--------------------------------
  !set xgrid and calculated array1 and array2
  dx=0.05_double
  do i=0,npt-1
     grid(i)=(i-npt/2)*dx !xgrid 
  end do
  ! Calculate 4th excited state wavefunction 
  !psi_n(x)=Norm*exp(-0.5*(x/sigma)**2)*Hpoly(n,x/sigma)
  !Norm=1/sqrt((sqrt(pi)*sigma)*(2**n*factorial(n)))=1/sqrt(sqrt(pi)*384)
  !sigma=sqrt(hbar/(mass*omega))=1
  !Y=sqrt(sqrt(pi)*384)
  psi=exp(-0.5_double*grid**2)*Hpoly(4,grid)!/Y
  Y=sqrt(sum(psi*psi))
  psi=psi/Y
  x=1.0_double
  call HOwf(4,grid,x,x,chi)
  !calcualte rmsd x
  x=sqrt(sum((chi-psi)**2)/real(npt))
  if(x.GT.2E-2)then
     write(*,*)'HOwf(n=4) RMSD=',x
     write(*,*)'RMSD error is greater than .02 on xrange [-3.4:3.4]&
          & with discritization of 128 grid points'
     open(123,file='HOwf.dat')
     do i=0,npt-1
        write(123,*)grid(i),chi(i),psi(i)
     end do
     close(123)
     stop
  end if

  write(*,*)'test HOwf returns 11th excited  wavefunction.....'
  !--------------------------------
  !harmonic oscillator constants
  !hbar=1
  !omega=1
  !mass=1
  !L=npt*dx=128*.07=8.96
  !En=hbar*omega(n+1/2)
  !V(x)=1/2*mass*omega^2*x^2-1
  !--------------------------------
  !set xgrid and calculated array1 and array2
  dx=0.07_double
  do i=0,npt-1
     grid(i)=(i-npt/2)*dx !xgrid 
  end do
  ! Calculate 4th excited state wavefunction 
  !psi_n(x)=Norm*exp(-0.5*(x/sigma)**2)*Hpoly(n,x/sigma)
  !Norm=1/sqrt((sqrt(pi)*sigma)*(2**n*factorial(n)))=arb
  !sigma=sqrt(hbar/(mass*omega))=1
  psi=exp(-0.5_double*grid**2)*Hpoly(11,grid)
  Y=sqrt(sum(psi*psi))
  psi=psi/Y
  x=1.0_double
  call HOwf(11,grid,x,x,chi)
  !calcualte rmsd x
  x=sqrt(sum((chi-psi)**2)/real(npt))
  if(x.GT.2E-2)then
     write(*,*)'HOwf(n=4) RMSD=',x
     write(*,*)'RMSD error is greater than .02 on xrange [-8.96:8.96]&
          & with discritization of 128 grid points'
     open(123,file='HOwf.dat')
     do i=0,npt-1
        write(123,*)grid(i),chi(i),psi(i)
     end do
     close(123)
     stop
  end if

  !test hilbert curve
  write(*,*)'test that hilbert curve is max at point npt/2 for ground state HO density'
  !set xgrid
  dx=0.05_double
  do i=0,npt-1
     grid(i)=(i-npt/2)*dx !xgrid 
  end do
  x=1.0_double
  !create HO ground state density
  call HOwf(0,grid,x,x,psi)
  do i=0,npt-1
     do j=0,npt-1
        rho(i,j)=psi(i)*psi(j)
     end do
  end do
  !transform density rho into linear hilbert curve
  Lrho=hilbertcurve(rho)
  !assert max of Lrho is at point npt**2/2
  !write(*,*) maxloc(Lrho)-1,size(Lrho),npt*npt/2
  if(maxloc(Lrho,1)-1.NE.npt*npt/2)then
     write(*,*)'maxloc of Lrho', maxloc(Lrho)-1,'is not npt*npt/2',npt*npt/2
     stop
  end if

!!$  write(*,*)'Diagnostic: write denisty along hilbert curve.'
!!$  !create HO |0><1| density
!!$  psi=grid
!!$  chi=-grid
!!$  !call HOwf(0,grid,x,x,psi)
!!$  !call HOwf(0,grid,x,x,chi)
!!$  do i=0,npt-1
!!$     do j=0,npt-1
!!$        rho(i,j)=psi(i)*chi(j)
!!$     end do
!!$  end do
!!$  !transform density rho into linear hilbert curve
!!$  Lrho=hilbertcurve(rho)
!!$  open(123,file='rho.dat')
!!$  do i=0,npt-1
!!$     do j=0,npt-1
!!$        write(123,*)i,j,rho(i,j)
!!$     end do
!!$     write(123,*)
!!$  end do
!!$  close(123)
!!$  open(123,file='Lrho.dat')
!!$  do i=0,npt*npt-1
!!$     write(123,*)i,Lrho(i)
!!$  end do
!!$  close(123)
!!$  call system('gnuplot Lrho.plt')

!!$  write(*,*)'test undoHilbertcurvemap returns proper matrix coordinates.'
!!$  k=Hilbertcurvemap(64,12,34) !input coordinate (12,34) for 64x64 matrix
!!$  call undoHilbertcurvemap(64,k,i,j)
!!$  if(i.NE.12.or.j.NE.34)then
!!$     write(*,*)'undoHilbertcurvemap does not return proper matrix coordinates.'
!!$     stop 
!!$  end if

  write(*,*)'test Hilbert curve index returns proper matrix coordinates.'
  !input coordinate (12,34) for 64x64 matrix  coord(1)=12
  N=64
  coord(1)=12
  coord(2)=34
  k=spacefilling_coord2index(N,coord,type='hilbert')
  coord=spacefilling_index2coord(N,k,type='hilbert')
  if(coord(1).NE.12.or.coord(2).NE.34)then
     write(*,*)'Hilbert curve mapping does not return proper matrix coordinates.'
     stop 
  end if

  write(*,*)'test Row-Major curve index returns proper matrix coordinates.'
  !input coordinate (12,34) for 64x64 matrix  coord(1)=12
  N=64
  coord(1)=12
  coord(2)=34
  k=spacefilling_coord2index(N,coord,type='rowmajor')
  coord=spacefilling_index2coord(N,k,type='rowmajor')
  if(coord(1).NE.12.or.coord(2).NE.34)then
     write(*,*)'Hilbert curve mapping does not return proper matrix coordinates.'
     stop 
  end if


  !test string module
  
  !test csv_io
  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CSV_IO_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the CSV_IO library.'
  
  call csv_io_test01 ( )
  call csv_io_test02 ( )
  !
  !  Terminate.
  !
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CSV_IO_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )
  





  write(*,*) 'ALL CORE TESTS PASSED!'
end subroutine coretests
!-----------------------

subroutine csv_io_test01 ( )

    !*****************************************************************************80
    !
    !! CSV_IO_TEST01 writes a variety of data items to a CSV file.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    29 November 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    implicit none

    integer ( kind = 4 ), parameter :: n = 10

    real ( kind = 4 ), dimension ( n ) :: column1 = (/ &
         1.1E+00, 2.2E+00, 3.3E+00, 4.4E+00,  5.5E+00, &
         6.6E+00, 7.7E+00, 8.8E+00, 9.9E+00, 10.10E+00 /)
    integer   ( kind = 4 ), dimension ( n ) :: column2 = (/ &
         1, 2, 3, 4, 5, 6, 7, 8, 9, 10 /)
    character ( len = 10 ), dimension ( n ) :: column3 = (/ &
         'one       ', 'two       ', 'three     ', 'four      ', 'five      ', &
         'six       ', 'seven     ', 'eight     ', 'nine      ', 'ten       ' /)
    real ( kind = 8 ), dimension ( n ) :: column4 = (/ &
         1.1D+00, 2.2D+00, 3.3D+00, 4.4D+00,  5.5D+00, &
         6.6D+00, 7.7D+00, 8.8D+00, 9.9D+00, 10.10D+00 /)
    character ( len = 80 ) :: csv_file_name = 'csv_4col_5row.csv'
    integer   ( kind = 4 ) csv_file_unit
    integer   ( kind = 4 ) i
    character ( len = 80 ) record

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSV_IO_TEST01'
    write ( *, '(a)' ) '  Test routines to create a CSV file.'

    call csv_file_open_write ( csv_file_name, csv_file_unit )
    !
    !  Write the header line.
    !
    record = ' '
    call csv_record_append_s ( 'R4', record )
    call csv_record_append_s ( 'I4', record )
    call csv_record_append_s ( 'S', record )
    call csv_record_append_s ( 'R8', record )

    call csv_file_header_write ( csv_file_name, csv_file_unit, record )

    do i = 1, n

       record = ' '

       call csv_record_append_r4 ( column1(i), record )
       call csv_record_append_i4 ( column2(i), record )
       call csv_record_append_s  ( column3(i), record )
       call csv_record_append_r8 ( column4(i), record )

       call csv_file_record_write (  csv_file_name, csv_file_unit, record )

    end do

    call csv_file_close_write ( csv_file_name, csv_file_unit )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Data written to "' // trim ( csv_file_name ) //'".'

    return
  end subroutine csv_io_test01
  subroutine csv_io_test02 ( )

    !*****************************************************************************80
    !
    !! CSV_IO_TEST02 reads a CSV file.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    29 November 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    implicit none

    character ( len = 80 ) :: csv_file_name = 'csv_4col_5row.csv'
    integer   ( kind = 4 ) csv_file_status
    integer   ( kind = 4 ) csv_file_unit
    integer   ( kind = 4 ) csv_record_status
    integer   ( kind = 4 ) i
    integer   ( kind = 4 ) line_num
    character ( len = 120 ) record
    integer   ( kind = 4 ) value_count

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSV_IO_TEST02'
    write ( *, '(a)' ) '  Read data from CSV file created by CSV_IO_TEST01.'

    call csv_file_line_count ( csv_file_name, line_num )

    write ( *, '(a,i8,a)' ) '  File contains ', line_num, ' lines.'

    call csv_file_open_read ( csv_file_name, csv_file_unit )

    do i = 1, line_num
       read ( csv_file_unit, '(a)', iostat = csv_file_status ) record
       write ( *, '(a)' ) i, trim ( record )
       call csv_value_count ( record, csv_record_status, value_count )
       write ( *, * ) i, value_count
    end do

    call csv_file_close_read ( csv_file_name, csv_file_unit )

    return
  end subroutine csv_io_test02
