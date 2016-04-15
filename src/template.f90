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

  !> Checks that the template object.
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
       !read(unit,*)this%XXX
       !**************************************
       
       !use reset to manage dynamic memory, reset sub-objects, and set random parameters
       call reset(this,state=1)
       
       !READ dynamic array values
       !***      example     ***
       !read(unit,*)(this%PPP(i),i=0,this%XXX-1)
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
       !this%XXX=123
       !**************************
       
       !Use reset to make a default object
       call reset(this,state=1)
    end if

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

    !Recompute any attribute values that might have evolved

    !****************************      Example     ******************************
    !*** attribute 'var' is always equall to the trace of the denisity matrix *** 
    ! this%var=0._double
    ! do istate=1,this%object%nstate
    !    this%var=this%var+this%den(istate,istate)
    ! end do
    !****************************************************************************

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
       if(state.EQ.0)then
          !nullify all pointers
          !******        Example - cleanup pointer attribute 'PPP'       ****
          !if(associated(this%PPP))nullify(this%PPP)
          !******************************************************************
          
          !kill all sub-objects
          !**** example **********
          !call kill(this%object)
          !***********************
          
          !un-initialized metiu object
          this%initialized=.false.
       else
          !allocate dynamic memory
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
          !this%object%XXX=123
          !**************************

          !reset all sub objects to correct any memory issues
          !***      example     ***
          !call reset(this%object,state=1)
          !************************

          !overwrite sub-object default dynamic parameters
          !***       example      ***
          !***set object pointer array values***
          !this%object%PPP(:)=XXX
          !**************************

          !declare initialization complete
          this%initialized=.true.
          
       end if
    end if
    
    !reset object based on current static parameters
    if(this%initialized)then

       !Sample Random parameters
       !***  Example - attribute 'var' samples a Gaussian random number
       ! this%var=gran()
       
!!$       !Reallocate all dynamic memory
!!$       !***  Example - cleanup pointer attribute 'PPP'     ***
!!$       !***            then reallocate memory              ***
!!$       !if(associated(this%PPP))nullify(this%PPP)
!!$       !allocate(this%PPP(0:this%npt-1))
!!$       !******************************************************
!!$       
!!$       !Reset dynamic memory values
!!$       !***  Example - set values in pointer 'PPP' to zero ***
!!$       !this%PPP(:)=0.0_double
!!$       !******************************************************

       !Resample sub-objects
       !**** example **********
       !call reset(this%object)
       !***********************

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


    !First, Scalar attributes
    !******          Example - Backup a scalar attribute            ******
    ! write(unit,*)this%var
    !*********************************************************************


    !Second, Dynamic attributes
    !***       Example - Backup an NxM matrix attribute                ***
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

    !Check all attributes are within acceptable values


    !**********   Example - check an object attribute 'that'  *********
    !call assert(check(this%that).EQ.0&
    !     ,msg='template_check: that sub-object failed check'&
    !     ,iostat=template_check)
    !if(template_check.NE.0)return
    !***********************************************************************
    
    !***   Example - check an integer attribute 'ndim' is well behaved   ***
    !call assert(check(this%ndim).EQ.0&
    !     ,msg='template_check: ndim failed check',iostat=template_check)
    !if(template_check.NE.0)return
    !***********************************************************************
 
    !*** Example - add a constrain that says 'ndim' can only be positive ***
    !call assert(this%ndim.GT.0&
    !     ,msg='template_check: ndim is not positive',iostat=template_check)
    !if(template_check.NE.0)return
    !***********************************************************************

    !***  Example - check a real valued attribute 'var' is well behaved  ***
    !call assert(check(this%var).EQ.0&
    !     ,msg='template_check: var failed check',iostat=template_check)
    !if(template_check.NE.0)return
    !***********************************************************************

    !***  Example - add a constrain that says 'var' can not be zero     ***
    !call assert(abs(this%var).GT.epsilon(this%var)&
    !     ,msg='template_check: var is tiny',iostat=template_check)
    !if(template_check.NE.0)return
    !***********************************************************************

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

    !================== consider the following tests ========================
    !-----  make tests -----
    !***          example          ****
    !write(*,*)'test make sets correct default values'
    !call make(this)
    !call assert(this%XXX.EQ.YYY,msg='template default XXX is not YYY')
    !call kill(this)
    !**********************************

    !----- kill tests -----
    !***          example          ****
    !write(*,*)'test kill cleans up dynamic memory and pointers'
    !call make(this)
    !call kill(this)
    !call assert(.not.associated(this%PPP),msg='template pointer PPP remains associated after killed.')
    !**********************************

    !----- backup tests -----
    !***          example          ****
    !write(*,*)'test attributes are stored properly in backup file'
    !call make(this)
    !this%var=XXX!First, manually set template attributes to non-default values
    !call system('rm -f template.tmpfile')
    !call backup(this,file='template.tmpfile')
    !call kill(this)
    !call make(this,file='template.tmpfile')
    !call system('rm -f template.tmpfile')
    !call assert(this%var.EQ.XXX)!Then, assert non default attribute values are conserved
    !call kill(this)    
    !**********************************

    !----- status tests -----

    !----- update tests -----

    !----- reset tests -----
    
    !----- fail cases -----

    !========================================================================



    write(*,*)'ALL template TESTS PASSED!'
  end subroutine template_test
  !-----------------------------------------

end module template_class

