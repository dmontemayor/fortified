program main
  implicit none

  integer,parameter::nsample=500
  integer::i
  integer::a,b,c

  !initialize state
  a=0
  b=0
  c=0

  open(123,file='xor.seqdat')
  do i=1,nsample
     b=a
     a=1
     if(rand().GE.0.5)a=0
     c=0
     if(a.NE.b)c=1
     write(123,*)a,c
  end do
  close(123)

end program main
