!>\brief
!! template class
!!\details
!! Use this template to help you start building a new object and methods.
!! Do not remove the \a make, \a kill, \a status, \a backup, \a update, \a reset, \a check, and \a test methods.
!! These methods must be present for the class to integrate properly.
!! You can leave these methods blank if you do not want those particular funcitionalities, although, this is not advised.
!! Follow the examples commented in the source code to define the various attributes of your class.
!! DO NOT OVERWRITE THIS FILE. Make a copy and title it \a myclass.f90 where 'myclass' is your choice
!! one-word title for your new class.
!! Finally, Replace instances of the string 'template' in this template with the same
!! one-word title 'myclass' you chose.
!<------------------------------------------------------------------------
module template_class
  use type_kinds
  implicit none
  private

  public::template, template_test
  public::check, make, kill, backup, update, reset, status, describe

  type template
     logical::initialized=.false.
     character(len=label)::name='template'
     
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
  end type template

  !> Creates the template object.
  interface make
     module procedure template_init
  end interface

  !> Destroys the template object.
  interface kill
     module procedure template_kill
  end interface

  !> Returns current state of the template object.
  interface status
     module procedure template_status
  end interface

  !> Returns a plain text description of the template object.
  interface describe
     module procedure template_describe
  end interface
  
  !> Backups the current state of the template object.
  interface backup
     module procedure template_backup
  end interface

  !> Recaluclates the template object.
  interface update
     module procedure template_update
  end interface

  !> Reinitializes the template object.
  interface reset
     module procedure template_reset
  end interface

  !> Checks that the template object.
  interface check
     module procedure template_check
  end interface

contains
  !======================================================================
  !> \brief Retruns a description of template as a string.
  !> \param[in] this is the template object.
  !======================================================================
  character(len=comment) function template_describe(this)
    type(template),intent(in)::this
    character(len=5)::FMT='(A)'

    write(template_describe,FMT)'No description for template has been provided.'
   
  end function template_describe

  !======================================================================
  !> \brief Creates and initializes the template object.
  !> \param this is the template object to be initialized.
  !> \param[in] file is an optional string containing the name of a previously backupd template file.
!!$  !> \remark If no input file is provided the user must manually initialize THIS using stout.
  !=====================================================================
  subroutine template_init(this,file)!,param)
    use filemanager
    use testing_class
    type(template),intent(inout)::this
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
    
       !check if file is of type template
       read(unit,*)header
       call assert(trim(header).EQ.'template',msg='template_init: bad input file header in file'//file)
       
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

  end subroutine template_init

  !======================================================================
  !> \brief Destroys the template object.
  !> \param this is the template object to be destroyed.
  !====================================================================
  subroutine template_kill(this)
    type(template),intent(inout)::this
 
!!$    call note('Begin template_kill.')
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
    !un-initialized template object
    this%initialized=.false.

  end subroutine template_kill

  !======================================================================
  !> \brief Computes the current state of template object.
  !> \param this is the template  object to be updated.
  !======================================================================
  subroutine template_update(this)
    type(template),intent(inout)::this

!!$    call Note('Begin template_update.')
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
!!$    !if(check(this).EQ.1)call stop('template_update: failed final check!')

  end subroutine template_update

  !======================================================================
  !> \brief Re-initiallizes the template object.
  !> \param this is the template  object to be re-initialized.
  !======================================================================
  subroutine template_reset(this)
    type(template),intent(inout)::this
!!$    integer(long)::istate,jstate
!!$
!!$    call Note('Begin template_reset.')
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
!!$    if(check(this).EQ.1)call stop('template_reset: failed final check!')

  end subroutine template_reset

  !======================================================================
  !> \brief Backups the current state of the template object to file.
  !> \param[in] this is the template  object to be updated.
  !> \param[in] file is a string containing the location of the backup file.
  !======================================================================
  subroutine template_backup(this,file)
    use filemanager
    type(template),intent(in)::this
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
!!$    call note('Begin template_backup.')
!!$    call Note('input file= '//file)
!!$    if(check(this).NE.0)then
!!$       call warn('template_backup: failed check.','not saving object.')
!!$    else
!!$
!!$       !assign a unique unit label
!!$       unit=newunit()
!!$
!!$       !open backup file
!!$       open(unit,file=file)
!!$
       !always write the data type on the first line
       write(unit,*)'template'
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
  end subroutine template_backup

  !======================================================================
  !> \brief Retrun the current state of template as a string.
  !> \param[in] this is the template object.
  !> \param[in] msg is an optional string message to annotate the status.
  !======================================================================
  character(len=line) function template_status(this,msg)
    type(template),intent(in)::this
    character*(*),intent(in),optional::msg
    character(len=5)::FMT='(A)'

    write(template_status,FMT)'template'
   
  end function template_status

 !======================================================================
  !> \brief Checks the template object.
  !> \param[in] this is the template object to be checked.
  !> \return 0 if all checks pass or exit at first failed check and returm non zero.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function template_check(this)
    use testing_class
    type(template),intent(in)::this

    !initiate with no problems found 
    template_check=0

    !check that object is initialized
    call assert(this%initialized,msg='template_check: template object not initialized.',iostat=template_check)
    if(template_check.NE.0)return

    !check the primitive
    !if(check(this%primitive))call stop('template_check: failed primitive check!')


    !********    Check all attributes are within acceptable values    *******!





    !***********************************************************************!
    !=======================================================================!
    !**********     Example - check an integer attribute 'ndim'    *********!
    !                                                                       !
    ! !check if integer 'ndim' is NAN (not a number)                        !
    ! if(this%ndim.NE.this%ndim)then                                        !
    !    call Warn('template_check: ndim not a number.')                    !
    !    template_check=1                                                   !
    !    return                                                             !
    ! end if                                                                !
    !                                                                       !
    ! !check if 'ndim' is too big to fit in its memory                      !
    ! if(abs(this%ndim).GE.huge(this%ndim))then                             !
    !    call Warn('template_check: ndim is too big.')                      !
    !    template_check=1                                                   !
    !    return                                                             !
    ! end if                                                                !
    !                                                                       !
    ! !add a constrain that says 'ndim' can only be positive                !
    ! if(this%ndim.LE.0)then                                                !
    !    call Warn('template_check: ndim not a positive integer.')          !
    !    template_check=1                                                   !
    !    return                                                             !
    ! end if                                                                !
    !                                                                       !
    !***********************************************************************!
    !=======================================================================!
    !**********    Example - check a real number attribute 'var'   *********!
    !                                                                       !
    ! !check if 'var' is not a number                                       !
    ! if(this%var.NE.this%var)then                                          !
    !    call Warn('template_check: var is not a number.')                  !
    !    template_check=1                                                   !
    !    return                                                             !
    ! end if                                                                !
    !                                                                       !
    ! !check if 'var' is too big to fit in its memory                       !
    ! if(abs(this%var).GE.huge(this%var))then                               !
    !    call Warn('template_check: var is too big.')                       !
    !    template_check=1                                                   !
    !   return                                                              !
    ! end if                                                                !
    !                                                                       !
    ! !add a constrain that says 'var' can not be zero:                     !
    ! !      'var' can not be smaller than the smallest computable value    !
    ! if(abs(this%var).LE.epsilon(this%var))then                            !
    !    call Warn('template_check: var is too small.')                     !
    !    template_check=1                                                   !
    !    return                                                             !
    ! end if                                                                !
    !                                                                       !
    !***********************************************************************!
    !=======================================================================!
    !********* Example - check an NxM matrix attribute 'matrix' ************!
    !                                                                       !
    ! !check that 'matrix' points to something                              !
    ! if(.not.associated(this%matrix))then                                  !
    !    call Warn('template_check: matrix memory not associated.')         !
    !    template_check=1                                                   !
    !    return                                                             !
    ! end if                                                                !
    !                                                                       !
    ! !check that 'matrix' has the right dimensions                         !
    ! if(size(this%matrix).NE.N*M)then                                      !
    !    call Warn('template_check: number of matrix elements not = N*M.')  !
    !    template_check=1                                                   !
    !    return                                                             !
    ! end if                                                                !
    !                                                                       !
    ! !check for NAN values in the matrix                                   !
    ! if(any(this%matrix.NE.this%matrix))then                               !
    !    call Warn('template_check: matirx has NAN values.')                !
    !    template_check=1                                                   !
    !    return                                                             !
    ! end if                                                                !
    !                                                                       !
    ! !check if any matrix element values are too big for thier memory      !
    ! if(any(abs(this%matirx).GT.huge(this%matirx)))then                    !
    !    call Warn('template_check: matrix has huge values.')               !
    !    mappingH_check=1                                                   !
    !    template_check=1                                                   !
    !    return                                                             !
    ! end if                                                                !
    !                                                                       !
    !***********************************************************************!

  end function template_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the template methods.
  !> \param[in] this is the template object whose methods will be excercised.
  !> \return Nothing if all tests pass or 1 and a stop for the first failed test.
  !> \remark Will stop after first failed check.
  !======================================================================
  subroutine template_test
    use testing_class
    use filemanager
    type(template)::this
    character(len=label)::string
    integer(long)::unit
    
    !verify template is compatible with current version
    include 'verification'

    !----- additional make tests -----
    write(*,*)'test make sets correct default values'
    call make(this)
    !call assert(this%XXX.EQ.YYY,msg='template default XXX is not YYY')
    call kill(this)

    
    !----- additional kill tests -----
    write(*,*)'test kill cleans up dynamic memory and pointers'
    !call assert(.not.associated(this%PPP),msg='template pointer PPP remains associated after killed.')


    !----- additional status tests -----


    !----- additional backup tests -----
    write(*,*)'test attributes are stored properly stored in backup file'
    call make(this)
    !manually set template attributes to non-default values
    call system('rm -f template.tmpfile')
    call backup(this,file='template.tmpfile')
    call kill(this)
    call make(this,file='template.tmpfile')
    !assert non default attribute values are conserved
    call kill(this)    

    
    !----- additional update tests -----

    
    !----- additional reset tests -----
    

    write(*,*)'ALL template TESTS PASSED!'
  end subroutine template_test
  !-----------------------------------------

end module template_class

