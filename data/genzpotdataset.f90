program genzpotdataset
  implicit none
  !Generates a 2D data set with prob distribution according to a 2D Z potential.
  !Output is written in plain text as two column space separated file.
  !Points are sampled from window -2<x<2 and -1<y<1 and accepted with
  !probability according to Z potential probability of the form
  !P(x,y)=exp(-3*(y-x*(x-1)*(x+1))**2)
  !The data set consitis of 500 data pairs.

  integer,parameter::nsample=500
  integer:: i
  real*8::x,y,z,P
  open(111,file='zpotdataset.dat')
  i=0
  do while(i.LT.nsample)
     call random_number(x)
     call random_number(y)
     call random_number(z)
     x=(x-0.5)*4
     y=(y-0.5)*2
     P=exp(-3*(y-x*(x-1)*(x+1))**2)
     if(P.GE.z)then
        i=i+1
        write(111,*)x,y
     end if
  end do
  write(111,*)'#zpot data set: with probaility dist P(x,y)=exp(-3*(y-x*(x-1)*(x+1))**2)'
  close(111)

  write(*,*)'Data set generated in zpotdataset.dat'
end program genzpotdataset
