!>\brief
!! wavelet class
!!\details
!<------------------------------------------------------------------------
module wavelet_class
  use type_kinds
  implicit none
  private

  public::wavelet, wavelet_test
  public::check, make, kill, backup, update, reset, status, describe

  type wavelet
     logical::initialized=.false.
     character(len=label)::name='wavelet'
     
!!$     type(primitive)::primitive
     !**********     Enter your derived type's attributes here     **********!


     !***********************************************************************!
     !=======================================================================!
     !***   Here are some example attributes your derived type may have   ***!
     ! type(primitive)::primitive2                     !an extra primitive   !
     ! integer(short)::label                           !a short integer      !
     ! integer(long)::ndim                             !a long integer       !
     ! real(double)::var                               !a real variable      !
     ! complex(double)::zed                            !a complex variable   !
     ! real(double),dimension(:,:),pointer::matrix     !a real matrix        !
     ! complex(double),dimension(:,:),pointer::Zmatrix !a complex matrix     !
     !***********************************************************************!
  end type wavelet

  !> Creates the wavelet object.
  interface make
     module procedure wavelet_init
  end interface

  !> Destroys the wavelet object.
  interface kill
     module procedure wavelet_kill
  end interface

  !> Returns current state of the wavelet object.
  interface status
     module procedure wavelet_status
  end interface

  !> Returns a plain text description of the wavelet object.
  interface describe
     module procedure wavelet_describe
  end interface
  
  !> Backups the current state of the wavelet object.
  interface backup
     module procedure wavelet_backup
  end interface

  !> Recaluclates the wavelet object.
  interface update
     module procedure wavelet_update
  end interface

  !> Reinitializes the wavelet object.
  interface reset
     module procedure wavelet_reset
  end interface

  !> Checks that the wavelet object.
  interface check
     module procedure wavelet_check
  end interface

contains
  !======================================================================
  !> \brief Retruns a description of wavelet as a string.
  !> \param[in] this is the wavelet object.
  !======================================================================
  character(len=comment) function wavelet_describe(this)
    type(wavelet),intent(in)::this
    character(len=5)::FMT='(A)'

    write(wavelet_describe,FMT)'The wavelet is the basis of any quantum&
         & mechanical object. It has attributes of mass and charge that&
         & are delocalized on a grid in N dimensional space.'
   
  end function wavelet_describe

  !======================================================================
  !> \brief Creates and initializes the wavelet object.
  !> \param this is the wavelet object to be initialized.
  !> \param[in] file is an optional string containing the name of a previously backuped wavelet file.
!!$  !> \remark If no input file is provided the user must manually initialize THIS using stout.
  !=====================================================================
  subroutine wavelet_init(this,file)!,param)
    use filemanager
    use testing_class
    type(wavelet),intent(inout)::this
    character*(*),intent(in),optional::file
    integer(long)::unit
    logical::fileisopen=.false.
    character(len=label)::header
    !integer(long),intent(in),optional::param
    !integer(long)::i

    !set scalar parameters
    if(present(file))then 
       
       !check input file    
       inquire(file=file,opened=fileisopen,number=unit)
       if(unit.LT.0)unit=newunit()
       if(.not.fileisopen)open(unit,file=file)
    
       !check if file is of type wavelet
       read(unit,*)header
       call assert(trim(header).EQ.'wavelet',msg='wavelet_init: bad input file header in file'//file)
       
       !read scalar parameters
       !read(unit,*)this%XXX
    else
       !set default scalar parameters
       !this%XXX=YYY
       !if(present(param))this%XXX=param
    end if

    !allocate dynamic arrays
    !if(associated(this%XXX))nullify(this%XXX)
    !allocate(this%XXX(0:thi%YYY-1))

    !set dynamic arrays
    if(present(file))then 
       !read dynamic arrays
       !read(unit,*)(this%XXX(i),i=0,this%YYY-1)
    else
       !set default array values
       !this%XXX=YYY
    end if
    
    !finished reading all attributes - now close store file
    if(.not.fileisopen)close(unit)

    !declare initialization complete
    this%initialized=.true.

  end subroutine wavelet_init

  !======================================================================
  !> \brief Destroys the wavelet object.
  !> \param this is the wavelet object to be destroyed.
  !====================================================================
  subroutine wavelet_kill(this)
    type(wavelet),intent(inout)::this
 
!!$    call note('Begin wavelet_kill.')
!!$
!!$    !kill the primitive
!!$    call kill(this%primitive)
!!$
!!$
!!$    !*******************       Nullify all pointers    **********************!
!!$
!!$ 
!!$
!!$
!!$
!!$    !************************************************************************!
!!$    !========================================================================!
!!$    !******        Example - cleanup matrix attribute 'matrix'       ********!
!!$    !                                                                        !
!!$    ! if(associated(this%matrix))nullify(this%matrix)                        !
!!$    !                                                                        !
!!$    !************************************************************************!
!!$
!!$
    !un-initialized wavelet object
    this%initialized=.false.

  end subroutine wavelet_kill

  !======================================================================
  !> \brief Computes the current state of wavelet object.
  !> \param this is the wavelet  object to be updated.
  !======================================================================
  subroutine wavelet_update(this)
    type(wavelet),intent(inout)::this

!!$    call Note('Begin wavelet_update.')
!!$    
!!$    !Primitives usually dont get updated
!!$    !    call update(this%primitive)    
!!$
!!$
!!$    !******   Recompute any attribute values that might have evolved   ******!
!!$
!!$
!!$
!!$
!!$
!!$    !************************************************************************!
!!$    !========================================================================!
!!$    !*****    Example - attribute 'var' is always equall to the trace   *****!
!!$    !                   of the primitive's denisity                          !
!!$    !                                                                        !
!!$    ! this%var=0._double                                                     !
!!$    ! do istate=1,this%primitive%nstate                                             !
!!$    !    this%var=this%var+this%primitive%den(istate,istate)                        !
!!$    ! end do                                                                 !
!!$    !                                                                        !
!!$    !************************************************************************!
!!$
!!$    !Usually leave out final check before we exit
!!$    !if(check(this).EQ.1)call stop('wavelet_update: failed final check!')

  end subroutine wavelet_update

  !======================================================================
  !> \brief Re-initiallizes the wavelet object.
  !> \param this is the wavelet  object to be re-initialized.
  !======================================================================
  subroutine wavelet_reset(this)
    type(wavelet),intent(inout)::this
!!$    integer(long)::istate,jstate
!!$
!!$    call Note('Begin wavelet_reset.')
!!$
!!$    !reset the primitive
!!$    call reset(this%primitive)
!!$
!!$    !****  Reset any attributes to satisfy re-initiation conditions   ****! 
!!$
!!$
!!$
!!$
!!$
!!$    !************************************************************************!
!!$    !========================================================================!
!!$    !********    Example - attribute 'var' is always initially a     ********!
!!$    !                      Gaussian random number                            !
!!$    !                                                                        !
!!$    ! this%var=gran()                                                        !
!!$    !                                                                        !
!!$    !************************************************************************!
!!$
!!$    !update now that we have changed some values and do a final check
!!$    call update(this)
!!$    if(check(this).EQ.1)call stop('wavelet_reset: failed final check!')

  end subroutine wavelet_reset

  !======================================================================
  !> \brief Backups the current state of the wavelet object to file.
  !> \param[in] this is the wavelet  object to be updated.
  !> \param[in] file is a string containing the location of the backup file.
  !======================================================================
  subroutine wavelet_backup(this,file)
    use filemanager
    type(wavelet),intent(in)::this
    character*(*),intent(in)::file
    integer(short)::unit
    logical::fileisopen
    integer(long)::i,j
    
    !check input file
    inquire(file=file,opened=fileisopen,number=unit)
    if(unit.LT.0)unit=newunit()
    if(.not.fileisopen)open(unit,file=file)

!!$    logical::usedunit      
!!$
!!$    call note('Begin wavelet_backup.')
!!$    call Note('input file= '//file)
!!$    if(check(this).NE.0)then
!!$       call warn('wavelet_backup: failed check.','not saving object.')
!!$    else
!!$
!!$       !assign a unique unit label
!!$       unit=newunit()
!!$
!!$       !open backup file
!!$       open(unit,file=file)
!!$
       !always write the data type on the first line
       write(unit,*)'wavelet'
!!$
!!$       !backup the primitive
!!$       call backup(this%primitive,file//'.primitive')
!!$
!!$       !write the location of the primitive
!!$       write(unit,*)quote(file//'.primitive')
!!$
!!$       !******      Backup below all the derived type's attributes       ******!
!!$       !******         in the order the MAKE command reads them         ******!
!!$
!!$
!!$
!!$
!!$
!!$       !*********************************************************************!
!!$       !=====================================================================!
!!$       !******      Example - Backup an attribute called 'var  '    ***********!
!!$       ! write(unit,*)this%var                                               !
!!$       !*********************************************************************!
!!$       !=====================================================================!
!!$       !***  Example - Backup an NxM matrix attribute called 'matrix'  ********!
!!$       ! write(unit,*)((this%matrix(i,j),j=1,M),i=1,N)                       !
!!$       !*********************************************************************!
!!$
!!$
!!$       !finished saving all attributes - now close backup file
       close(unit)
!!$    end if
  end subroutine wavelet_backup

  !======================================================================
  !> \brief Retrun the current state of wavelet as a string.
  !> \param[in] this is the wavelet object.
  !> \param[in] msg is an optional string message to annotate the status.
  !======================================================================
  character(len=line) function wavelet_status(this,msg)
    type(wavelet),intent(in)::this
    character*(*),intent(in),optional::msg
    character(len=5)::FMT='(A)'

    write(wavelet_status,FMT)'wavelet'
   
  end function wavelet_status

 !======================================================================
  !> \brief Checks the wavelet object.
  !> \param[in] this is the wavelet object to be checked.
  !> \return 0 if all checks pass or exit at first failed check and returm non zero.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function wavelet_check(this)
    use testing_class
    type(wavelet),intent(in)::this

    !initiate with no problems found 
    wavelet_check=0

    !check that object is initialized
    call assert(this%initialized,msg='wavelet_check: wavelet object not initialized.',iostat=wavelet_check)
    if(wavelet_check.NE.0)return

    !check the primitive
    !if(check(this%primitive))call stop('wavelet_check: failed primitive check!')


    !********    Check all attributes are within acceptable values    *******!





    !***********************************************************************!
    !=======================================================================!
    !**********     Example - check an integer attribute 'ndim'    *********!
    !                                                                       !
    ! !check if integer 'ndim' is NAN (not a number)                        !
    ! if(this%ndim.NE.this%ndim)then                                        !
    !    call Warn('wavelet_check: ndim not a number.')                    !
    !    wavelet_check=1                                                   !
    !    return                                                             !
    ! end if                                                                !
    !                                                                       !
    ! !check if 'ndim' is too big to fit in its memory                      !
    ! if(abs(this%ndim).GE.huge(this%ndim))then                             !
    !    call Warn('wavelet_check: ndim is too big.')                      !
    !    wavelet_check=1                                                   !
    !    return                                                             !
    ! end if                                                                !
    !                                                                       !
    ! !add a constrain that says 'ndim' can only be positive                !
    ! if(this%ndim.LE.0)then                                                !
    !    call Warn('wavelet_check: ndim not a positive integer.')          !
    !    wavelet_check=1                                                   !
    !    return                                                             !
    ! end if                                                                !
    !                                                                       !
    !***********************************************************************!
    !=======================================================================!
    !**********    Example - check a real number attribute 'var'   *********!
    !                                                                       !
    ! !check if 'var' is not a number                                       !
    ! if(this%var.NE.this%var)then                                          !
    !    call Warn('wavelet_check: var is not a number.')                  !
    !    wavelet_check=1                                                   !
    !    return                                                             !
    ! end if                                                                !
    !                                                                       !
    ! !check if 'var' is too big to fit in its memory                       !
    ! if(abs(this%var).GE.huge(this%var))then                               !
    !    call Warn('wavelet_check: var is too big.')                       !
    !    wavelet_check=1                                                   !
    !   return                                                              !
    ! end if                                                                !
    !                                                                       !
    ! !add a constrain that says 'var' can not be zero:                     !
    ! !      'var' can not be smaller than the smallest computable value    !
    ! if(abs(this%var).LE.epsilon(this%var))then                            !
    !    call Warn('wavelet_check: var is too small.')                     !
    !    wavelet_check=1                                                   !
    !    return                                                             !
    ! end if                                                                !
    !                                                                       !
    !***********************************************************************!
    !=======================================================================!
    !********* Example - check an NxM matrix attribute 'matrix' ************!
    !                                                                       !
    ! !check that 'matrix' points to something                              !
    ! if(.not.associated(this%matrix))then                                  !
    !    call Warn('wavelet_check: matrix memory not associated.')         !
    !    wavelet_check=1                                                   !
    !    return                                                             !
    ! end if                                                                !
    !                                                                       !
    ! !check that 'matrix' has the right dimensions                         !
    ! if(size(this%matrix).NE.N*M)then                                      !
    !    call Warn('wavelet_check: number of matrix elements not = N*M.')  !
    !    wavelet_check=1                                                   !
    !    return                                                             !
    ! end if                                                                !
    !                                                                       !
    ! !check for NAN values in the matrix                                   !
    ! if(any(this%matrix.NE.this%matrix))then                               !
    !    call Warn('wavelet_check: matirx has NAN values.')                !
    !    wavelet_check=1                                                   !
    !    return                                                             !
    ! end if                                                                !
    !                                                                       !
    ! !check if any matrix element values are too big for thier memory      !
    ! if(any(abs(this%matirx).GT.huge(this%matirx)))then                    !
    !    call Warn('wavelet_check: matrix has huge values.')               !
    !    mappingH_check=1                                                   !
    !    wavelet_check=1                                                   !
    !    return                                                             !
    ! end if                                                                !
    !                                                                       !
    !***********************************************************************!

  end function wavelet_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the wavelet methods.
  !> \param[in] this is the wavelet object whose methods will be excercised.
  !> \return Nothing if all tests pass or 1 and a stop for the first failed test.
  !> \remark Will stop after first failed check.
  !======================================================================
  subroutine wavelet_test
    use testing_class
    use filemanager
    type(wavelet)::this
    character(len=label)::string
    integer(long)::unit
    
    !verify wavelet is compatible with current version
    include 'verification'

    !----- additional make tests -----
    write(*,*)'test make sets correct default values'
    call make(this)
    !call assert(this%XXX.EQ.YYY,msg='wavelet default XXX is not YYY')
    call kill(this)

    
    !----- additional kill tests -----
    write(*,*)'test kill cleans up dynamic memory and pointers'
    !call assert(.not.associated(this%PPP),msg='wavelet pointer PPP remains associated after killed.')


    !----- additional status tests -----


    !----- additional backup tests -----
    write(*,*)'test attributes are stored properly stored in backup file'
    call make(this)
    !manually set wavelet attributes to non-default values
    call system('rm -f wavelet.tmpfile')
    call backup(this,file='wavelet.tmpfile')
    call kill(this)
    call make(this,file='wavelet.tmpfile')
    !assert non default attribute values are conserved
    call kill(this)    

    
    !----- additional update tests -----

    
    !----- additional reset tests -----
    

    write(*,*)'ALL wavelet TESTS PASSED!'
  end subroutine wavelet_test
  !-----------------------------------------

end module wavelet_class

