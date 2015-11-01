program main
  implicit none

  integer,parameter::nsample=500
  integer::i
  integer::a,b,c

  open(123,file='xor.dat')
  do i=1,nsample
     a=floor(rand()*2)
     b=floor(rand()*2)
     c=0
     if(a.NE.b)c=1
     write(123,*)a,b,c
  end do
  close(123)

end program main
