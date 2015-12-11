!>\brief
!! template class
!!\details
!! Use this template to help you start building a new object and methods.
!! Do not remove the \a make, \a kill, \a display, \a store, \a update, \a reset, \a check, and \a test methods.
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
  public::make, kill, display, store, update, reset, check

  type template
     logical::initialized=.false.
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

  !> Displays the current state of the template object.
  interface display
     module procedure template_display
  end interface

  !> Stores the current state of the template object.
  interface store
     module procedure template_store
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
  !> \brief Creates and initializes the template object.
  !> \param this is the template object to be initialized.
  !> \param[in] file is an optional string containing the name of a previously stored template file.
!!$  !> \remark If no input file is provided the user must manually initialize THIS using stout.
  !=====================================================================
  subroutine template_init(this,file)
    type(template),intent(inout)::this
    character*(*),intent(in),optional::file
!!$    character(len=title)::filetype
!!$    character(len=path)::infile
!!$
!!$    integer::unit
!!$    logical::usedefault,usedunit
!!$


    !declare initialization complete
    this%initialized=.true.




!!$    call Note('Begin template_init.')
!!$
!!$    !check if input file is present and valid
!!$    if(present(file))then
!!$
!!$       !check if file is there
!!$       if(check(file).EQ.1)call stop('template_init: cannot find input file '//file)
!!$
!!$       !assign a unique unit label
!!$       unit=newunit()
!!$
!!$       !open the file
!!$       open(unit,file=file)
!!$
!!$       !read the file type - should always be on the first line
!!$       read(unit,*)filetype
!!$       filetype=adjustl(filetype)
!!$
!!$       !check if input file is the right kind of file
!!$       if(trim(filetype).NE.'template')& 
!!$            call Stop('template_init: input file is not valid.')
!!$
!!$    end if
!!$    
!!$    ! prepare primative
!!$    if(present(file))then
!!$       read(unit,*)infile
!!$       infile=adjustl(infile)
!!$       call make(this%primitive,trim(infile))
!!$    else
!!$       call make(this%primitive)
!!$    end if
!!$
!!$    !******     Initiate all your derived type's attributes below      ******!
!!$    !************************************************************************!
!!$
!!$
!!$
!!$
!!$    !************************************************************************!
!!$    !========================================================================!
!!$    !********           Example - Setup attribute 'var'            **********!
!!$    !                                                                        !
!!$    ! if(present(file))then                                                  !
!!$    !    read(unit,*)this%var          !Read from input file if present      !
!!$    ! else                                                                   !
!!$    !    write(*,*)'Please enter VAR.' !manually enter information           !
!!$    !    read(*,*)this%var                                                   !
!!$    ! end if                                                                 !
!!$    !                                                                        !
!!$    !************************************************************************!
!!$    !========================================================================!
!!$    !********   Example - Setup a matrix attribute called 'matrix'   ********!
!!$    !                                                                        !
!!$    ! if(associated(this%matrix))nullify(this%matrix) !cleanup memory        !
!!$    ! allocate(this%matirx(N,M))                      !create an NxM matrix  !
!!$    !                                                                        !
!!$    ! !Now fill that matrix with data from the input file or manually        !
!!$    ! if(present(file))then                                                  !
!!$    !    read(unit,*)((this%matrix(i,j),j=1,M),i=1,N) !read M loop first     !
!!$    ! else                                                                   !
!!$    !    write(*,*)'Enter value of matrix element:'                          !
!!$    !    do i=1,N                                                            !
!!$    !       do j=1,M                                                         !
!!$    !          write(*,*)'index=',i,j      !prompt user the matrix index     !
!!$    !          read(*,*)this%matrix(i,j)   !read data manually from keyboard !
!!$    !       end do                                                           !
!!$    !    end do                                                              !
!!$    ! end if                                                                 !
!!$    !                                                                        !
!!$    !************************************************************************!
!!$
!!$
!!$    !finished reading data - now close input file
!!$    if(present(file))close(unit)
!!$
!!$    !declare that initialization is complete
!!$    this%initialized=.true.
!!$
!!$    !update (or reset if brand new) this derived type
!!$    if(present(file))then
!!$       call update(this)
!!$    else
!!$       call reset(this)
!!$    end if
!!$
!!$    !do a final check before we exit
!!$    if(check(file).EQ.1)call stop('template_init: failed final check!')

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
  !> \brief Stores the current state of the template object to file.
  !> \param[in] this is the template  object to be updated.
  !> \param[in] file is a string containing the location of the store file.
  !======================================================================
  subroutine template_store(this,file)
    use filemanager

    type(template),intent(in)::this
    character*(*),intent(in)::file
    integer(short)::unit

!!$    logical::usedunit      
!!$
!!$    call note('Begin template_store.')
!!$    call Note('input file= '//file)
!!$    if(check(this).NE.0)then
!!$       call warn('template_store: failed check.','not saving object.')
!!$    else
!!$
       !assign a unique unit label
       unit=newunit()

       !open store file
       open(unit,file=file)

       !always write the data type on the first line
       write(unit,*)'template'
!!$
!!$       !store the primitive
!!$       call store(this%primitive,file//'.primitive')
!!$
!!$       !write the location of the primitive
!!$       write(unit,*)quote(file//'.primitive')
!!$
!!$       !******      Store below all the derived type's attributes       ******!
!!$       !******         in the order the MAKE command reads them         ******!
!!$
!!$
!!$
!!$
!!$
!!$       !*********************************************************************!
!!$       !=====================================================================!
!!$       !******      Example - Store an attribute called 'var  '    ***********!
!!$       ! write(unit,*)this%var                                               !
!!$       !*********************************************************************!
!!$       !=====================================================================!
!!$       !***  Example - Store an NxM matrix attribute called 'matrix'  ********!
!!$       ! write(unit,*)((this%matrix(i,j),j=1,M),i=1,N)                       !
!!$       !*********************************************************************!
!!$
!!$
!!$       !finished saving all attributes - now close store file
       close(unit)
!!$    end if
  end subroutine template_store


  !======================================================================
  !> \brief Displays the template object.
  !> \param[in] this is the template object.
  !> \param[in] msg is an optional string message to preface the displayed object.
  !======================================================================
  subroutine template_display(this,msg)
    type(template),intent(in)::this
    character*(*),intent(in),optional::msg

!!$    call Note('Begin template_display.')
!!$
!!$    if(check(this).NE.0)then
!!$       call warn('template_display: failed check','displaying nothing.')
!!$       return
!!$    end if
!!$
!!$    write(Dunit,*)'____________________________________________________'
!!$    write(Dunit,*)'-------------------   template   -------------------'
!!$    if(present(msg))write(Dunit,*)msg
!!$    write(Dunit,*)'____________________________________________________'
!!$    
!!$    
!!$    !****    Display the derived type's attributes here if you want   ****!
!!$    
!!$    
!!$    
!!$    !*********************************************************************!
!!$    !=====================================================================!
!!$    !*******        Example display scalar attribute 'var'      **********!
!!$    !                                                                     !
!!$    ! write(Dunit,*)'VAR=',this%var                                       !
!!$    !                                                                     !
!!$    !*********************************************************************!
!!$    !=====================================================================!
!!$    !*******    Example display NxM matrix attribute 'matrix'   **********!
!!$    !                                                                     !
!!$    ! write(Dunit,*)'MATRIX=',this%matrix   !a simple example             !
!!$    !                                                                     !
!!$    ! write(Dunit,*)'MATRIX='               !a better example:            !
!!$    ! do i=1,N                          !write each row on a new line     !
!!$    !    write(Dunit,*)i,(this%matrix(i,j),j=1,M)                         !
!!$    ! end do                                                              !
!!$    !*********************************************************************!
!!$
!!$    call display(this%primitive,msg='template primitive')
!!$    write(Dunit,*)'===================================================='

  end subroutine template_display

  !======================================================================
  !> \brief Checks the template object.
  !> \param[in] this is the template object to be checked.
  !> \return Nothing if all checks pass or 1 and a warn for the first failed check.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function template_check(this)
    use testing_class
    type(template),intent(in)::this

    !initiate with no problems found 
    template_check=0
    !call Note('Checking template.')

    !check that object is initialized
    call assert(this%initialized,msg='template_check: template object not initialized.',iostat=template_check)
    if(template_check.NE.0)return

!!$    !check the primitive
!!$    if(check(this%primitive))call stop('template_check: failed primitive check!')
!!$
!!$
!!$    !********    Check all attributes are within acceptable values    *******!
!!$
!!$
!!$
!!$
!!$
!!$    !**************************************************************************************!
!!$    !======================================================================================!
!!$    !**********     Example - check an integer attribute 'ndim'    ************************!
!!$    !                                                                                      !
!!$    ! !check if integer 'ndim' is NAN (not a number)                                       !
!!$    ! if(this%ndim.NE.this%ndim)then                                                       !
!!$    !    call Warn('template_check: ndim not a number.')                                   !
!!$    !    template_check=1                                                                  !
!!$    !    return                                                                            !
!!$    ! end if                                                                               !
!!$    !                                                                                      !
!!$    ! !check if 'ndim' is too big to fit in its memory                                     !
!!$    ! if(abs(this%ndim).GE.huge(this%ndim))then                                            !
!!$    !    call Warn('template_check: ndim is too big.')                                     !
!!$    !    template_check=1                                                                  !
!!$    !    return                                                                            !
!!$    ! end if                                                                               !
!!$    !                                                                                      !
!!$    ! !add a constrain that says 'ndim' can only be positive                               !
!!$    ! if(this%ndim.LE.0)then                                                               !
!!$    !    call Warn('template_check: ndim not a positive integer.')                         !
!!$    !    template_check=1                                                                  !
!!$    !    return                                                                            !
!!$    ! end if                                                                               !
!!$    !                                                                                      !
!!$    !**************************************************************************************!
!!$    !======================================================================================!
!!$    !**********    Example - check a real number attribute 'var'   ************************!
!!$    !                                                                                      !
!!$    ! !check if 'var' is not a number                                                      !
!!$    ! if(this%var.NE.this%var)then                                                         !
!!$    !    call Warn('template_check: var is not a number.')                                 !
!!$    !    template_check=1                                                                  !
!!$    !    return                                                                            !
!!$    ! end if                                                                               !
!!$    !                                                                                      !
!!$    ! !check if 'var' is too big to fit in its memory                                      !
!!$    ! if(abs(this%var).GE.huge(this%var))then                                              !
!!$    !    call Warn('template_check: var is too big.')                                      !
!!$    !    template_check=1                                                                  !
!!$    !   return                                                                             !
!!$    ! end if                                                                               !
!!$    !                                                                                      !
!!$    ! !add a constrain that says 'var' can not be zero:                                    !
!!$    ! !      'var' can not be smaller than the smallest computable value                   !
!!$    ! if(abs(this%var).LE.epsilon(this%var))then                                           !
!!$    !    call Warn('template_check: var is too small.')                                    !
!!$    !    template_check=1                                                                  !
!!$    !    return                                                                            !
!!$    ! end if                                                                               !
!!$    !                                                                                      !
!!$    !**************************************************************************************!
!!$    !======================================================================================!
!!$    !*********          Example - check an NxM matrix attribute 'matrix'        ***********!
!!$    !                                                                                      !
!!$    ! !check that 'matrix' points to something                                             !
!!$    ! if(.not.associated(this%matrix))then                                                 !
!!$    !    call Warn('template_check: matrix memory not associated.')                        !
!!$    !    template_check=1                                                                  !
!!$    !    return                                                                            !
!!$    ! end if                                                                               !
!!$    !                                                                                      !
!!$    ! !check that 'matrix' has the right dimensions                                        !
!!$    ! if(size(this%matrix).NE.N*M)then                                                     !
!!$    !    call Warn('template_check: number of matrix elements not = N*M.')                 !
!!$    !    template_check=1                                                                  !
!!$    !    return                                                                            !
!!$    ! end if                                                                               !
!!$    !                                                                                      !
!!$    ! !check for NAN values in the matrix                                                  !
!!$    ! if(any(this%matrix.NE.this%matrix))then                                              !
!!$    !    call Warn('template_check: matirx has NAN values.')                               !
!!$    !    template_check=1                                                                  !
!!$    !    return                                                                            !
!!$    ! end if                                                                               !
!!$    !                                                                                      !
!!$    ! !check if any matrix element values are too big for thier memory                     !
!!$    ! if(any(abs(this%matirx).GT.huge(this%matirx)))then                                   !
!!$    !    call Warn('template_check: matrix has huge values.')                              !
!!$    !    mappingH_check=1                                                                  !
!!$    !    template_check=1                                                                  !
!!$    !    return                                                                            !
!!$    ! end if                                                                               !
!!$    !                                                                                      !
!!$    !**************************************************************************************!

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
    type(template)::this
    character(len=label)::string
    
    write(*,*)'test template can be checked prior to being made.'
    call assert(check(this).EQ.1,msg='check template object does not return 1 prior to make.')

    write(*,*)'test template can be created.'
    call make(this)
    call assert(check(this).EQ.0,msg='template object was not created properly.')
    write(*,*)'test template kill method sets initiallization flag to false.'
    call kill(this)
    call assert(.not.this%initialized,msg='template object remains initialized after killed.')

    write(*,*)'test template can be stored'
    call make(this)
    call system('rm -f testtemplate.tmpfile')
    call store(this,file='testtemplate.tmpfile')
    call assert('testtemplate.tmpfile',msg='template store file was not created.')
    call system('rm -f testtemplate.tmpfile')
    call kill(this)

    write(*,*)'test template can be created with store file'
    call make(this)
    call store(this,file='testtemplate.tmpfile')
    call kill(this)
    call make(this,file='testtemplate.tmpfile')
    call assert(check(this).EQ.0,msg='template object was not created properly from storefile.')
    call system('rm -f testtemplate.tmpfile')
    
    write(*,*)'test savefile begins with template object name in first line'
    call make(this)
    call system('rm -f testtemplate.tmpfile')
    call store(this,file='testtemplate.tmpfile')
    open(123,file='testtemplate.tmpfile')
    read(123,*)string
    call assert(trim(string).EQ.'template',msg='save file does not have template label on first line.')
    call system('rm -f testtemplate.tmpfile')
    call kill(this)
    
    write(*,*)'test template can be updated'
    call make(this)
    call update(this)
    call assert(check(this).EQ.0,msg='template object was not updated properly.')
    call kill(this)
    
    write(*,*)'test template can be resetted'
    call make(this)
    call reset(this)
    call assert(check(this).EQ.0,msg='template object was not resetted properly.')
    call kill(this)
    

    write(*,*)'ALL template TESTS PASSED!'
  end subroutine template_test
  !-----------------------------------------

end module template_class

