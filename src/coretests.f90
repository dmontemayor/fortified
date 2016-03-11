subroutine coretests
  use type_kinds
  use testing_class
  use omp_lib
  use filemanager
  
  implicit none

  integer(short)::ierr
  logical::pass

  real(double)::x
  real(double),pointer::ptr(:) 
  
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
  ptr(:)=0._double
  ptr(:)=1/ptr
  if(check(ptr).NE.-1)then
     write(*,*)'check double real pointer does not return -1 for poor behaved elements'
     stop
  end if
  
  write(*,*) 'ALL CORE TESTS PASSED!'
end subroutine coretests
!-----------------------
