program gencubicdataset
  implicit none
  !Generates a labeled data set.
  !Output is written in plain text as two column space separated file.
  !Independent data {x} is in the first column, dependent classifiers {y} are in second column.
  !The generating function is y=Ax**3+B, where A is a constant parameter 
  !and B is a gaussian random number of variance s**2.
  !The data set consitis of 500 data pairs. x is sampled uniformly from (-5:5).
  !Successfull fitting of the data should recover A with error s.

  real*8,parameter::A=2.0
  real*8,parameter::s=10.0
  integer,parameter::nsample=500
  integer:: i
  real*8::x,z1,z2,B,pi

  real*8::Ebar,Evar

  Ebar=0.
  Evar=0.
  pi=acos(-1D0)
  
  open(111,file='cubicdataset.dat')
  do i=1,nsample
     call random_number(x)
     x=(x-.5)*10
     
     call random_number(z1)
     call random_number(z2)
     B=sqrt(-2D0*log(z1))*cos(2D0*pi*z2)
     z2=B*s
     
     write(111,*)x,A*x**3+z2
     Ebar=Ebar+z2
     Evar=Evar+z2*z2
     
  end do
  Ebar=Ebar/real(nsample)
  Evar=Evar/real(nsample)-Ebar**2
  
  write(*,*)'Error=',sqrt(Evar)
  write(111,*)'#cubic data set: y=Ax**3+gran*s, A=2.0, s=10.0, Specific Error=',sqrt(Evar)
  close(111)

  write(*,*)'Data set generated in cubicdataset.dat'
end program gencubicdataset
