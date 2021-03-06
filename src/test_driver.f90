program test
  use layer_class
  implicit none

  character(len=50)::cmd
  integer::len,ierr


  call get_command_argument(1,cmd,len,ierr)
  if(ierr.GT.0)write(*,*)'Test name cannot be read properly! (try: test help) Program will stop.'
  if(ierr.LT.0)write(*,*)trim(cmd)//' Test name was truncated! (try: test help) Program will stop.'
  if(ierr.NE.0)stop

  select case(trim(cmd))
  case ('help')
     write(*,*)'Summary: Unit testing for class methods.'
     write(*,*)'Usage: test classname'
     write(*,*)'available classname:'
     write(*,*)
     write(*,*)'help'
     write(*,*)'all'
     write(*,*)'core'
     write(*,*)'template'
     write(*,*)'layer'
     write(*,*)
  case('all')
     call coretests
     call template_test
     call layer_test
  case('core')
     call coretests
  case('template')
     call template_test
  case('layer')
     call layer_test
  case default
     write(*,*)'Unknown test (try: test help)'
  end select

  
  write(*,*) 'ALL TESTS PASSED!'
  
end program test
