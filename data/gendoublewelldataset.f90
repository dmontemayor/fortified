program gendoublewelldataset
  implicit none
  !Generates a 2D data set with prob distribution according to a 2D double well.
  !probability of finding point is the sum of two gaussians of unit variance
  !centered at points (1,2) and (-1,-2).
  !Output is written in plain text as two column space separated file.
  !Points are sampled from window -10<x<10 and -10<y<10 and accepted with
  !probability according to double gaussian probability of the form
  !P(x,y)=exp(-1/2*((x-1)**2+(y-2)**2))+exp(-1/2*((x+1)**2+(y+2)**2))
  !The data set consitis of 500 data pairs.

  real*8,parameter::x1=1.0,x2=-1.0
  real*8,parameter::y1=1.0,y2=-1.0
  integer,parameter::nsample=500
  integer:: i
  real*8::x,y,z,P,N
  open(111,file='doublewelldataset.dat')
  N=1+exp(-.5*((1+1)**2+(2+2)**2))
  i=0
  do while(i.LT.nsample)
     call random_number(x)
     call random_number(y)
     call random_number(z)
     x=(x-0.5)*20
     y=(y-0.5)*20
     P=exp(-.5*((x-1)**2+(y-2)**2))+exp(-.5*((x+1)**2+(y+2)**2))
     if(P.GE.z*N)then
        i=i+1
        write(111,*)x,y
     end if
  end do
  write(111,*)'#doublewell data set: with probaility dist P(x,y)=exp(-1/2*((x-1)**2+(y-2)**2))+exp(-1/2*((x+1)**2+(y+2)**2))'
  close(111)

  write(*,*)'Data set generated in doublewelldataset.dat'
end program gendoublewelldataset
