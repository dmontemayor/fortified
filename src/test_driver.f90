program test
  implicit none

  character(len=50)::cmd
  integer::len,ierr


  call get_command_argument(1,cmd,len,ierr)
  if(ierr.GT.0)write(*,*)'Test name cannot be read properly! (try: coretest help) Program will stop.'
  if(ierr.LT.0)write(*,*)trim(cmd)//' Test name was truncated! (try: coretest help) Program will stop.'
  if(ierr.NE.0)stop

  select case(trim(cmd))
  case ('help')
     write(*,*)'Summary: Unit testing for core functions and class methods.'
     write(*,*)'Usage: test testname'
     write(*,*)'available testname:'
     write(*,*)
     write(*,*)'help'
     write(*,*)'all'
     write(*,*)'core'
     write(*,*)
  case('all')
     call coretests
  case('core')
     call coretests
  case default
     write(*,*)'Unknown test (try: test help)'
  end select

  
  write(*,*) 'ALL TESTS PASSED!'
  
end program test
