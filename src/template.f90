!!!!!!!! silcrow ATTRIBUTE LIST !!!!!!!!
!!! use silcrow wizard to pre-generate
!!! attribute source, checks, and tests.
!!! usage:
!!!
!§ ATTRIBUTE TYPE MEMORY RANK DEPENDENCIES
!!!
!!! '!§' identifies the silcrow directive
!!! 'ATTRIBUTE' is the attribute name
!!! 'TYPE' can be any of the known type_kinds
!!!     or classes.
!!! 'MEMORY' can be 'static' or 'dynamic'.
!!! 'RANK' is an integer number of dimensions
!!! 'DEPENDENCIES' is a list parameters of
!!!     length RANK ued to define the size.
!!!     Can be an integer for static memory
!!!     or another static attribute of type
!!!     integer for dynamic memory.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>\brief
!! template class
!!\details
!! DO NOT OVERWRITE THIS FILE.
!! Make a copy and title it \a myclass.f90
!! where 'myclass' is your choice one-word title for your new class.
!! Then replace instances of the string 'template' with the same
!! one-word title 'myclass' you chose before you begin building your class.
!!
!! Use this template to help you start building a new object.
!! Do not remove the \a make, \a kill, \a status, \a backup, \a update,
!! \a reset, \a check, and \a test methods.
!! These methods must be present for the class to integrate properly.
!! Augment these methods to ensure the new object behaves properly.
!! These methods will operate on the object attributes which need to be
!! defined in the template type.
!! Commented examples are provided to help build the object class.
!!
!! A good start is to compose a description of the class and enter it
!! in the description method. Throuought development ensure this
!! description acurately represents the object as it will be the primary
!! resource used by others to understand the object and how to use it.
!!
!! Development should be test driven. Motivated by the object description,
!! create unit tests in the test method that challenge the other object
!! methods to behave properly. Write the tests first, ensure they fail,
!! then modify the class to correct test failure. Don't modify or remove
!! previous tests.
!!
!! The check method does not challenge the object to behave properly,
!! rather it mearly checks that all of the objects attributes are within
!! acceptable parameters. A unit test may ensure the check function returns
!! an error after purposely assigning a bad attribute value but the check
!! fucntion itself is not a unit test. Write checks after the unit test.
!<------------------------------------------------------------------------
module template_class
  use type_kinds
  implicit none
  private

  public::template, template_test
  public::describe, check, make, kill, backup, update, reset, status

  type template
     logical::initialized=.false.
     character(len=label)::name='template'

     !***************      Enter template attributes here     ***************!


     !***********************************************************************!
     !=======================================================================!
     !***   Here are some example attributes your derived type may have   ***!
     ! type(object)::that                              !an other type        !
     ! integer(short)::ierr                            !a short integer      !
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

  !> Returns the current state of the template object.
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

  !> Checks the template object.
  interface check
     module procedure template_check
  end interface

contains
  !======================================================================
  !> \brief Retruns a description of template as a string.
  !> \param[in] THIS is the template object.
  !======================================================================
  character(len=comment) function template_describe(this)
    type(template),intent(in)::this
    character(len=5)::FMT='(A)'

    write(template_describe,FMT)'No description for template has been provided.'

  end function template_describe

  !======================================================================
  !> \brief Creates and initializes template.
  !! \param THIS is the template object.
  !! \param[in] FILE is an optional string containing the name of a
  !! previously backuped template file.
  !=====================================================================
  subroutine template_init(this,file)
    use filemanager
    use testing_class
    type(template),intent(inout)::this
    character*(*),intent(in),optional::file
    integer(long)::unit
    logical::fileisopen=.false.
    character(len=label)::header
    character(len=path)::infile

    !initialize all sub-objects
    !***      example     ***
    !call make(this%object)
    !************************

    !check input file
    if(present(file))then

       !check input file
       inquire(file=file,opened=fileisopen,number=unit)
       if(unit.LT.0)unit=newunit()
       if(.not.fileisopen)open(unit,file=file)

       !check if file is of type template
       read(unit,*)header
       call assert(trim(header).EQ.this%name,msg='template_init: bad input file header in file'//file)

       !read static parameters
       !***             example            ***
       !***read scalar parameters from file***
       !read(unit,*)this%NNN
       !**************************************

       !use reset to manage dynamic memory
       !, reset sub-objects, and set calculated variable seeds
       call reset(this,state=1)

       !READ dynamic array values
       !***      example     ***
       !read(unit,*)(this%PPP(i),i=0,this%NNN-1)
       !************************

       !READ sub-objects
       !***      example     ***
       !read(unit,*)infile
       !call make(this%object,file=trim(infile))
       !************************

       !finished reading all attributes - now close backup file
       if(.not.fileisopen)close(unit)

       !declare initialization complete
       this%initialized=.true.
    else
       !Set static parameters to default settings
       !***       example      ***
       !***set scalar parameter***
       !this%NNN=123
       !**************************

       !Use reset to make a default object
       call reset(this,state=1)
    end if


    !finished making now update object
    call update(this)

  end subroutine template_init

  !======================================================================
  !> \brief Destroys the template object.
  !> \param THIS is the template object to be destroyed.
  !> \remarks kill is simply the reset method passed with a null flag
  !====================================================================
  subroutine template_kill(this)
    type(template),intent(inout)::this

    call reset(this,0)

  end subroutine template_kill

  !======================================================================
  !> \brief Computes the current state of template object.
  !> \param THIS is the template  object to be updated.
  !======================================================================
  subroutine template_update(this)
    type(template),intent(inout)::this

    !Recompute calculated variables that might have evolved

    !***************************      Example     *****************************
    !** attribute 'var' is always equall to the trace of the denisity matrix **
    ! this%var=0._double
    ! do istate=1,this%object%nstate
    !    this%var=this%var+this%den(istate,istate)
    ! end do
    !**************************************************************************

  end subroutine template_update

  !======================================================================
  !> \brief Re-initiallizes the template object.
  !> \param THIS is the template  object to be re-initialized.
  !> \param STATE is an optional integer:
  !>        when 0, will create a null state by deallocating all dynamic
  !>        memory and returning the object to an un-initiallized state;
  !>        when not 0, will return the object to the default settings;
  !>        when not present, object will reset based on current scalar
  !>        parameters.
  !======================================================================
  subroutine template_reset(this,state)
    type(template),intent(inout)::this
    integer(long),intent(in),optional::STATE

    if(present(state))then
      !un-initialize template object until properly reset
      this%initialized=.false.
       if(state.EQ.0)then
          !nullify all pointers
          !******        Example - cleanup pointer attribute 'PPP'       ****
          !if(associated(this%PPP))nullify(this%PPP)
          !******************************************************************

          !kill all sub-objects
          !**** example **********
          !call kill(this%object)
          !***********************

          !set all scalar parameters to error values
          !**** example **********
          !this%nstate=-1
          !***********************

       else
         if(checkparam(this).EQ.0) then
           !reallocate dynamic memory
           !***  Example - cleanup pointer attribute 'PPP'     ***
           !***            then reallocate memory              ***
           !if(associated(this%PPP))nullify(this%PPP)
           !allocate(this%PPP(0:this%npt-1))
           !******************************************************

           !Set default dynamic memory values
           !***  Example - set values in pointer 'PPP' to zero ***
           !this%PPP(:)=0.0_double
           !******************************************************

           !overwrite sub-object default static parameters
           !***       example      ***
           !***set object static parameter***
           !this%object%NNN=123
           !**************************

           !reset all sub objects to correct any memory issues
           !***      example     ***
           !call reset(this%object,state=1)
           !************************

           !overwrite sub-object default dynamic array
           !***       example      ***
           !***set object pointer array values***
           !this%object%PPP(:)=XXX
           !**************************

           !declare initialization complete
           this%initialized=.true.
         end if
       end if
    end if

    !reset object based on current static parameters
    if(this%initialized)then

       !Sample calculated variable seeds
       !***  Example - attribute 'var' samples a Gaussian random number
       ! this%var=gran()

       !Resample sub-objects
       !**** example **********
       !call reset(this%object)
       !***********************

       !update now that object is fully reset and not in null state
       call update(this)

    end if


  end subroutine template_reset

  !======================================================================
  !> \brief Backups the current state of the template object to file.
  !> \param[in] THIS is the template  object to be updated.
  !> \param[in] FILE is a string containing the location of the backup file.
  !======================================================================
  subroutine template_backup(this,file)
    use filemanager
    use string
    use testing_class
    type(template),intent(in)::this
    character*(*),intent(in)::file
    integer(short)::unit
    logical::fileisopen
    integer(long)::i,j

    !check input file
    inquire(file=file,opened=fileisopen,number=unit)
    if(unit.LT.0)unit=newunit()
    if(.not.fileisopen)open(unit,file=file)

    !check template object
    call assert(check(this).EQ.0,msg='template object does not pass check.')

    !always write the data type on the first line
    write(unit,*)'template'

    !******      Backup below all the derived type's attributes       ****
    !******         in the order the MAKE method reads them           ****


    !First, Scalar parameters
    !******          Example - Backup a scalar parameter            ******
    ! write(unit,*)this%var
    !*********************************************************************


    !Second, Dynamic arrays
    !***       Example - Backup an NxM matrix                          ***
    ! write(unit,*)((this%matrix(i,j),j=1,M),i=1,N)
    !*********************************************************************


    !Last,objects
    !******              Example - Backup an object            ***********
    ! call backup(this%object,file//'.object')
    ! write(unit,*)quote(file//'.object')!write the object location
    !*********************************************************************


    !finished writing all attributes - now close backup file
    close(unit)
  end subroutine template_backup

  !======================================================================
  !> \brief Retrun the current state of template as a string.
  !> \param[in] THIS is the template object.
  !> \param[in] MSG is an optional string message to annotate the status.
  !======================================================================
  character(len=line) function template_status(this,msg)
    type(template),intent(in)::this
    character*(*),intent(in),optional::msg
    character(len=5)::FMT='(A)'

    !Edit the status prompt to suit your needs
    write(template_status,FMT)'template status is currently not available'

  end function template_status
  !======================================================================
  !> \brief Checks the template object parameters are good.
  !> \remark Will stop program after first failed check.
  !======================================================================
  integer(short)function checkparam(this)
    use testing_class
    type(template),intent(in)::this

    !initiate with no problems found
    checkparam=0

    !***   Example - check an integer parameter 'ndim' is well behaved   ***
    !call assert(check(this%ndim).EQ.0&
    !     ,msg='checkparam: ndim failed check',iostat=checkparam)
    !if(checkparam.NE.0)return
    !***********************************************************************

    !*** Example - add a constrain that says 'ndim' can only be positive ***
    !call assert(this%ndim.GT.0&
    !     ,msg='checkparam: ndim is not positive',iostat=checkparam)
    !if(checkparam.NE.0)return
    !***********************************************************************

    !***  Example - check a real valued parameter 'var' is well behaved  ***
    !call assert(check(this%var).EQ.0&
    !     ,msg='checkparam: var failed check',iostat=checkparam)
    !if(checkparam.NE.0)return
    !***********************************************************************

    !***  Example - add a constrain that says 'var' can not be zero     ***
    !call assert(abs(this%var).GT.epsilon(this%var)&
    !     ,msg='checkparam: var is tiny',iostat=checkparam)
    !if(checkparam.NE.0)return
    !***********************************************************************
  end function checkparam
  !======================================================================
  !> \brief Checks the template object.
  !> \param[in] THIS is the template object to be checked.
  !> \return 0 if all checks pass or exit at first failed check and returm non zero.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function template_check(this)
    use testing_class
    type(template),intent(in)::this

    !initiate with no problems found
    template_check=0

    !check that object is initialized
    call assert(this%initialized&
         ,msg='template_check: template object not initialized.'&
         ,iostat=template_check)
    if(template_check.NE.0)return

    !check that object has correct name
    call assert(this%name.EQ.'template'&
         ,msg='template_check: template name is not set.'&
         ,iostat=template_check)
    if(template_check.NE.0)return

    !check all parameters are within acceptable values
    call assert(checkparam(this).EQ.0&
         ,msg='template_check: unacceptable parameters found!'&
         ,iostat=template_check)
    if(template_check.NE.0)return

    !Check all sub-objects
    !**********   Example - check an object attribute 'that'  *********
    !call assert(check(this%that).EQ.0&
    !     ,msg='template_check: that sub-object failed check'&
    !     ,iostat=template_check)
    !if(template_check.NE.0)return
    !***********************************************************************

    !Check dynamic attributes
    !***  Example - check a real valued pointer attribute 'matrix'       ***
    !***            is well behaved                                      ***
    !call assert(check(this%matrix).EQ.0&
    !     ,msg='template_check: matrix failed check',iostat=template_check)
    !if(template_check.NE.0)return
    !***********************************************************************

    !********* Example - check an NxM matrix has right dimensions **********
    !call assert(size(this%matrix).EQ.N*M&
    !     ,msg='template_check: number of matrix elements not = N*M.'&
    !     ,iostat=template_check)
    !if(template_check.NE.0)return
    !***********************************************************************

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

    !======Unit testing=====
    ! Testing of individual components. Typically done by
    ! the developer and not by testers as it requires
    ! detailed knowledge of the internal program design.


    !======Functional testing======
    ! This type of testing igores the internal parts
    ! and focuses on the output as per requirements.


    !=====End-to-end testing=====
    ! Tests that mimics real-world use with physically
    ! reasonable inputs or well established examples.


    !=====Load testing=====
    ! Test performance behavior under various loads.
    ! Examples include calculation speed as a function
    ! of parallel processors used, or as a funtion of
    ! system size.


    !================== consider the following tests for =====================
    !=========================================================================

    !==================       static parameters             ==================

    !!write(*,*)'test kill template breaks static parameter NNN.....'
    !!call make(this) !make template
    !!call kill(this)
    !!call assert(check(this%NNN).NE.0,&
    !!     msg='kill template does not break static parameter NNN')
    !! check method should be extended to optionally check individual
    !! attributes.
    !!Usage: check(this,element='NNN')
    !!check method
    !!    if(present(element))
    !!      select case(element)
    !!        case (element='NNN')
    !!             ...
    !!        case (element='XXX')
    !!             ...
    !!        default
    !!             throw error unkown element
    !!      end select
    !!    else
    !!      run all checks
    !!    end if

    !write(*,*)'test static parameter NNN is stored properly in backup file'
    !call make(this) !make template
    !this%NNN=MMM    !manually set static parameter to non-default value
    !call system('rm -f template.tmpfile*') !remove any previous backup file(s)
    !call backup(this,file='template.tmpfile') !create backup file
    !call kill(this) !destroy template
    !call make(this,file='template.tmpfile') !make new template from backup file
    !call system('rm -f template.tmpfile*') !remove backup file(s)
    !!assert non default parameter is conserved
    !call assert(this%NNN.EQ.MMM,&
    !     msg='template static parameter NNN is not stored properly')
    !call kill(this)    !destroy template to clean up memory

    !write(*,*)'test make sets correct default value for static parameter NNN'
    !call make(this) !make template
    !call assert(this%NNN.EQ.MMM,&
    !     msg='template default static parameter NNN is not MMM')
    !call kill(this)

    !write(*,*)'test edge case for static parameter NNN breaks template'
    !call make(this) !make template
    !this%NNN=MMM
    !call assert(check(this).NE.0,&
    !     msg='edge case value MMM for static parameter NNN &
    !     &does not break template')
    !call kill(this)

    !==================    dynamic arrays and pointers      ==================

    !write(*,*)'test make allocates memory for dyanamic pointer array PPP'
    !call make(this)
    !call assert(associated(this%PPP),msg='template dynamic pointer array &
    !     &PPP remains associated after killed.')
    !call kill(this)

    !write(*,*)'test kill deallocates memory for dynamic pointer array PPP'
    !call make(this)
    !call kill(this)
    !call assert(.not.associated(this%PPP),msg='template dynamic pointer array &
    !     &PPP remains associated after killed.')

    !write(*,*)'test dynamic pointer array PPP is saved properly in backup file'
    !call make(this) !make template
    !this%PPP=YYY    !manually set dynamic pointer array to non-default value
    !call system('rm -f template.tmpfile*') !remove any previous backup file(s)
    !call backup(this,file='template.tmpfile') !create backup file
    !call kill(this) !destroy template
    !call make(this,file='template.tmpfile') !make new template from backup file
    !call system('rm -f template.tmpfile*') !remove backup file(s)
    !!assert non default parameter is conserved
    !call assert(all(this%PPP.EQ.YYY),&
    !     msg='template dynamic pointer array PPP is not stored properly')
    !call kill(this)    !destroy template to clean up memory

    !write(*,*)'test make allocates dynamic pointer array PPP of default size'
    !call make(this) !make template
    !!assert dynamic pointer array has proper size for all dimensions
    !call assert(size(this%PPP,J).EQ.N),&
    !     msg='template dynamic pointer array PPP is not of size N for &
    !     &dimension J')
    !call assert(size(this%PPP,I).EQ.N),&
    !     msg='template dynamic pointer array PPP is not of size N for &
    !     &dimension I')
    !call kill(this)    !destroy template to clean up memory

    !write(*,*)'test dynamic pointer array PPP can be resized by adjusting &
    !     & static parameter NNN then reseting with state=1'
    !call make(this) !make template
    !this%NNN=MMM !adjust static parameter NNN to non-default value
    !call reset(this,state=1) !reset template to reallocate dynamic memory
    !!assert dynamic pointer array size has changed properly
    !call assert(size(this%PPP,J).EQ.MMM),&
    !     msg='template dynamic pointer array PPP did not change size upon &
    !     &reset with state=1')
    !call kill(this)    !destroy template to clean up memory

    !write(*,*)'test edge case for dynamic pointer array PPP breaks template'
    !call make(this) !make template
    !this%PPP=YYY !set edge case value
    !assert edge case value breaks template
    !call assert(check(this).NE.0,&
    !     msg='edge case value YYY for dynamic pointer array PPP &
    !     &does not break template')
    !call kill(this)    !destroy template to clean up memory

    !==================        calculated variables         ==================

    !write(*,*)'test make sets default value for calculated variable XXX'
    !call make(this) !make template
    !call assert(this%XXX.EQ.YYY,&
    !     msg='template default calculated variable XXX is not YYY')
    !call kill(this)

    !write(*,*)'test edge case for calculated variable XXX breaks template'
    !call make(this) !make template
    !this%XXX=YYY
    !call assert(check(this).NE.0,&
    !     msg='edge case value YYY for calculated variable XXX &
    !     &does not break template')
    !call kill(this)
    !========================================================================





    write(*,*)'ALL template TESTS PASSED!'
  end subroutine template_test
  !-----------------------------------------

end module template_class
