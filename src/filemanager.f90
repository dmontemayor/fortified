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
  integer(short) function checkfile(file,type)
    character*(*),intent(in)::file
    character*(*),optional,intent(in)::type
    character(len=label)::header
    integer(long)::unit
    logical::there
    
    inquire(file=file,exist=there)
    if(there)then
       checkfile=0
       if(present(type))then
          unit=newunit()
          open(unit,file=file)
          read(unit,*)header
          close(unit)
          if(trim(header).NE.type)checkfile=-1
       end if
    else
       checkfile=1
    end if
    return
  end function checkfile

end module filemanager
