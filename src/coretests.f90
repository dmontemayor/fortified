subroutine coretests
  use type_kinds
  use testing_class
  use omp_lib

  implicit none

  integer(short)::ierr
  logical::pass
  real(single)::dum_sing
  real(double)::dum_doub

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


  !--------- OPEN MP ----------
  write(*,*)'checking omp is available.....'
  call assert(omp_get_num_procs().GT.0,msg=&
       'number of available processors is not greater than 0')
  call assert(omp_get_max_threads().GT.0,msg=&
       'number of available threads is not greater than 0')
  call assert(omp_get_max_threads().GE.omp_get_num_procs(),msg=&
       'number of threads is less than number of processors available'&
       ,iostat=ierr)
  if(ierr.NE.0)then
     write(*,*)'try setting the environment variable with'
     write(*,*)'export OMP_NUM_THREADS=',omp_get_num_procs()
  end if



  
!!$  !-----SINGLE REAL ASSERT TESTS-----
!!$  write(*,*)'checking single real assert can be called with only two single reals.....'
!!$  call assert(1_single,1_single)
!!$
!!$  write(*,*)'checking single real assert exits gracefully for not equal condition.....'
!!$  call assert(1_single,2_single,iostat=ierr)
!!$
!!$  write(*,*)'checking single real assert writes msg for not equal condition.....'
!!$  write(*,*)'--------------------------------------------------------------'
!!$  call assert(1_single,2_single,msg='single real assert is writting you a message.',iostat=ierr)
!!$  write(*,*)'--------------------------------------------------------------'
!!$  write(*,*)'Do you see a message in between the lines above? Answer true or false.'
!!$  read(*,*)pass
!!$  if(.not.pass)then
!!$     write(*,*)'User cannot see a message passed by single real assert.'
!!$     stop
!!$  end if
!!$
!!$  write(*,*)'checking single real assert does not write error message for equal condition.....'
!!$  call assert(1_single,1_single,msg='Error: single real assert writes error message for equal condition.')
!!$  
!!$  write(*,*)'checking single real assert iostat returns 0 for equal condition.....' 
!!$  call assert(1_single,1_single,iostat=ierr)
!!$  if(ierr.NE.0)then
!!$     write(*,*)'single real assert: iostat option does not return 0 for equal condition.'
!!$     stop
!!$  end if
!!$
!!$  write(*,*)'checking single real assert iostat returns 1 for not equal condition.....' 
!!$  call assert(1_single,2_single,iostat=ierr)
!!$  if(ierr.NE.1)then
!!$     write(*,*)'single real assert: iostat option does not return 1 for not equal condition'
!!$     stop
!!$  end if
!!$
!!$  write(*,*)'checking single real assert returns -1 for huge numbers.....' 
!!$  call assert(Huge(1_single),Huge(1_single),iostat=ierr)
!!$  if(ierr.NE.-1)then
!!$     write(*,*)'single real assert: iostat option does not return -1 for huge numbers'
!!$     stop
!!$  end if



!!$  write(*,*)'checking single real assert returns 0 when two not equal numbers are within tolerance.....'
!!$  call assert(1_single,2_single,1_single,iostat=ierr)
!!$  if(ierr.NE.0)then
!!$     write(*,*)'single real assert: does not return 0 when two not equal numbers are within tolerance'
!!$     stop
!!$  end if

!!$  write(*,*)'checking single real assert returns 1 when two not equal numbers are not within tolerance.....'
!!$  call assert(-1_single,1_single,1_single,iostat=ierr)
!!$  if(ierr.NE.1)then
!!$     write(*,*)'single real assert: does not return 1 when two not equal numbers are not within tolerance'
!!$     stop
!!$  end if


  write(*,*) 'ALL CORE TESTS PASSED!'
end subroutine coretests
!-----------------------
