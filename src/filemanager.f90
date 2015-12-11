module filemanager
!> \brief
!! IO utility
!! \details 
!! Controls Fortran unit assignments for read and write statements and checks
!! if files exist.
!< --------------------------------------------------------------
  use type_kinds
  implicit none
  
  private
  public:: newunit,check

  interface check
     module procedure checkfile
  end interface

contains
  !----------------------------
  integer(long) function newunit()
    logical::usedunit  
    newunit=1000
    usedunit=.true.
    do while (usedunit)
       newunit=newunit+1
       inquire(newunit,opened=usedunit)
    end do
    return
  end function newunit
  !-----------------------------
  integer(short) function checkfile(file)
    character*(*),intent(in)::file
    logical::there
    
    checkfile=0
    inquire(file=file,exist=there)
    if(.not.there)then
       checkfile=1
    end if
    return
  end function checkfile

end module filemanager
