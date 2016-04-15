!>\brief
!! sbham class
!!\details
!! DO NOT OVERWRITE THIS FILE. 
!! Make a copy and title it \a myclass.f90
!! where 'myclass' is your choice one-word title for your new class.
!! Then replace instances of the string 'sbham' with the same
!! one-word title 'myclass' you chose before you begin building your class.
!!
!! Use this sbham to help you start building a new object.
!! Do not remove the \a make, \a kill, \a status, \a backup, \a update,
!! \a reset, \a check, and \a test methods.
!! These methods must be present for the class to integrate properly.
!! Augment these methods to ensure the new object behaves properly.
!! These methods will operate on the object attributes which need to be
!! defined in the sbham type.
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
module sbham_class
  use type_kinds
  implicit none
  private

  public::sbham, sbham_test
  public::describe, check, make, kill, backup, update, reset, status

  type sbham
     logical::initialized=.false.
     character(len=label)::name='sbham'
     
     !***************      Enter sbham attributes here     ***************!


     !***********************************************************************!
     !=======================================================================!
     !***   Here are some example attributes your derived type may have   ***!
     ! type(primitive)::primitive                      !an other type        !
     ! integer(short)::ierr                            !a short integer      !
     ! integer(long)::ndim                             !a long integer       !
     ! real(double)::var                               !a real variable      !
     ! complex(double)::zed                            !a complex variable   !
     ! real(double),dimension(:,:),pointer::matrix     !a real matrix        !
     ! complex(double),dimension(:,:),pointer::Zmatrix !a complex matrix     !
     !***********************************************************************!
  end type sbham

  !> Creates the sbham object.
  interface make
     module procedure sbham_init
  end interface

  !> Destroys the sbham object.
  interface kill
     module procedure sbham_kill
  end interface

  !> Returns current state of the sbham object.
  interface status
     module procedure sbham_status
  end interface

  !> Returns a plain text description of the sbham object.
  interface describe
     module procedure sbham_describe
  end interface
  
  !> Returns the current state of the sbham object.
  interface backup
     module procedure sbham_backup
  end interface

  !> Recaluclates the sbham object.
  interface update
     module procedure sbham_update
  end interface

  !> Reinitializes the sbham object.
  interface reset
     module procedure sbham_reset
  end interface

  !> Checks that the sbham object.
  interface check
     module procedure sbham_check
  end interface

contains
  !======================================================================
  !> \brief Retruns a description of sbham as a string.
  !> \param[in] THIS is the sbham object.
  !======================================================================
  character(len=comment) function sbham_describe(this)
    type(sbham),intent(in)::this
    character(len=5)::FMT='(A)'

    write(sbham_describe,FMT)'No description for sbham has been provided.'
   
  end function sbham_describe

  !======================================================================
  !> \brief Creates and initializes sbham.
  !! \param THIS is the sbham object.
  !! \param[in] FILE is an optional string containing the name of a
  !! previously backuped sbham file.
  !=====================================================================
  subroutine sbham_init(this,file)
    use filemanager
    use testing_class
    type(sbham),intent(inout)::this
    character*(*),intent(in),optional::file
    integer(long)::unit
    logical::fileisopen=.false.
    character(len=label)::header

    !assign scalar parameters
    if(present(file))then 
       
       !check input file
       inquire(file=file,opened=fileisopen,number=unit)
       if(unit.LT.0)unit=newunit()
       if(.not.fileisopen)open(unit,file=file)
    
       !check if file is of type sbham
       read(unit,*)header
       call assert(trim(header).EQ.this%name,msg='sbham_init: bad input file header in file'//file)
       
       !read scalar parameters from file
       !***      example     ***
       !read(unit,*)this%XXX
       !************************
    else
       !set default scalar parameters
       !***      example     ***
       !this%XXX=YYY
       !************************
    end if

    !allocate dynamic arrays from scalar parameters
    !***      example     ***
    !if(associated(this%PPP))nullify(this%PPP)
    !allocate(this%PPP(0:this%XXX-1))
    !************************

    !assign dynamic array values
    if(present(file))then 
       !read dynamic arrays from file
       !***      example     ***
       !read(unit,*)(this%PPP(i),i=0,this%XXX-1)
       !************************
    else
       !set default array values
       !***      example     ***
       !this%PPP=YYY
       !************************
    end if
    
    !finished reading all attributes - now close backup file
    if(.not.fileisopen)close(unit)

    !declare initialization complete
    this%initialized=.true.

  end subroutine sbham_init

  !======================================================================
  !> \brief Destroys the sbham object.
  !> \param THIS is the sbham object to be destroyed.
  !====================================================================
  subroutine sbham_kill(this)
    type(sbham),intent(inout)::this
 
    !nullify all pointers
    !******        Example - cleanup pointer attribute 'PPP'       ****
    !if(associated(this%PPP))nullify(this%PPP)
    !******************************************************************

    !un-initialized sbham object
    this%initialized=.false.

  end subroutine sbham_kill

  !======================================================================
  !> \brief Computes the current state of sbham object.
  !> \param THIS is the sbham  object to be updated.
  !======================================================================
  subroutine sbham_update(this)
    type(sbham),intent(inout)::this

    !Recompute any attribute values that might have evolved

    !****************************      Example     ******************************
    !*** attribute 'var' is always equall to the trace of the denisity matrix *** 
    ! this%var=0._double
    ! do istate=1,this%primitive%nstate
    !    this%var=this%var+this%den(istate,istate)
    ! end do
    !****************************************************************************

  end subroutine sbham_update

  !======================================================================
  !> \brief Re-initiallizes the sbham object.
  !> \param THIS is the sbham  object to be re-initialized.
  !======================================================================
  subroutine sbham_reset(this)
    type(sbham),intent(inout)::this

    !Reset any attributes to satisfy re-initiation conditions

    !**********  Example - attribute 'var' is always initially a     ********
    !**********            Gaussian random number                    ******** 
    ! this%var=gran()
    !************************************************************************

  end subroutine sbham_reset

  !======================================================================
  !> \brief Backups the current state of the sbham object to file.
  !> \param[in] THIS is the sbham  object to be updated.
  !> \param[in] FILE is a string containing the location of the backup file.
  !======================================================================
  subroutine sbham_backup(this,file)
    use filemanager
    type(sbham),intent(in)::this
    character*(*),intent(in)::file
    integer(short)::unit
    logical::fileisopen
    integer(long)::i,j
    
    !check input file
    inquire(file=file,opened=fileisopen,number=unit)
    if(unit.LT.0)unit=newunit()
    if(.not.fileisopen)open(unit,file=file)
    
    !always write the data type on the first line
    write(unit,*)'sbham'

    !******      Backup below all the derived type's attributes       ****
    !******         in the order the MAKE method reads them           ****




    !******              Example - Backup an object            ***********
    ! call backup(this%primitive,file//'.primitive')
    ! write(unit,*)quote(file//'.primitive')!write the object location
    !*********************************************************************


    !******          Example - Backup a scalar attribute            ******
    ! write(unit,*)this%var
    !*********************************************************************


    !***       Example - Backup an NxM matrix attribute                ***
    ! write(unit,*)((this%matrix(i,j),j=1,M),i=1,N)
    !*********************************************************************
    
    
    !finished writing all attributes - now close backup file
    close(unit)  
  end subroutine sbham_backup
  
  !======================================================================
  !> \brief Retrun the current state of sbham as a string.
  !> \param[in] THIS is the sbham object.
  !> \param[in] MSG is an optional string message to annotate the status.
  !======================================================================
  character(len=line) function sbham_status(this,msg)
    type(sbham),intent(in)::this
    character*(*),intent(in),optional::msg
    character(len=5)::FMT='(A)'
    
    !Edit the status prompt to suit your needs
    write(sbham_status,FMT)'sbham status is currently not available'
    
  end function sbham_status

 !======================================================================
  !> \brief Checks the sbham object.
  !> \param[in] THIS is the sbham object to be checked.
  !> \return 0 if all checks pass or exit at first failed check and returm non zero.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function sbham_check(this)
    use testing_class
    type(sbham),intent(in)::this

    !initiate with no problems found 
    sbham_check=0

    !check that object is initialized
    call assert(this%initialized&
         ,msg='sbham_check: sbham object not initialized.'&
         ,iostat=sbham_check)
    if(sbham_check.NE.0)return

    !check that object has correct name
    call assert(this%name.EQ.'sbham'&
         ,msg='sbham_check: sbham name is not set.'&
         ,iostat=sbham_check)
    if(sbham_check.NE.0)return

    !Check all attributes are within acceptable values


    !**********   Example - check an object attribute 'primitive'  *********
    !call assert(check(this%primitive).EQ.0&
    !     ,msg='sbham_check: primitive sub-object failed check'&
    !     ,iostat=sbham_check)
    !if(sbham_check.NE.0)return
    !***********************************************************************
    
    !***   Example - check an integer attribute 'ndim' is well behaved   ***
    !call assert(check(this%ndim).EQ.0&
    !     ,msg='sbham_check: ndim failed check',iostat=sbham_check)
    !if(sbham_check.NE.0)return
    !***********************************************************************
 
    !*** Example - add a constrain that says 'ndim' can only be positive ***
    !call assert(this%ndim.GT.0&
    !     ,msg='sbham_check: ndim is not positive',iostat=sbham_check)
    !if(sbham_check.NE.0)return
    !***********************************************************************

    !***  Example - check a real valued attribute 'var' is well behaved  ***
    !call assert(check(this%var).EQ.0&
    !     ,msg='sbham_check: var failed check',iostat=sbham_check)
    !if(sbham_check.NE.0)return
    !***********************************************************************

    !***  Example - add a constrain that says 'var' can not be zero     ***
    !call assert(abs(this%var).GT.epsilon(this%var)&
    !     ,msg='sbham_check: var is tiny',iostat=sbham_check)
    !if(sbham_check.NE.0)return
    !***********************************************************************

    !***  Example - check a real valued pointer attribute 'matrix'       ***
    !***            is well behaved                                      ***
    !call assert(check(this%matrix).EQ.0&
    !     ,msg='sbham_check: matrix failed check',iostat=sbham_check)
    !if(sbham_check.NE.0)return
    !***********************************************************************

    !********* Example - check an NxM matrix has right dimensions **********
    !call assert(size(this%matrix).EQ.N*M&
    !     ,msg='sbham_check: number of matrix elements not = N*M.'&
    !     ,iostat=sbham_check)
    !if(sbham_check.NE.0)return
    !***********************************************************************

  end function sbham_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the sbham methods.
  !> \param[in] this is the sbham object whose methods will be excercised.
  !> \return Nothing if all tests pass or 1 and a stop for the first failed test.
  !> \remark Will stop after first failed check.
  !======================================================================
  subroutine sbham_test
    use testing_class
    use filemanager
    type(sbham)::this
    character(len=label)::string
    integer(long)::unit
    
    !verify sbham is compatible with current version
    include 'verification'

    !================== consider the following tests ========================
    !-----  make tests -----
    !***          example          ****
    !write(*,*)'test make sets correct default values'
    !call make(this)
    !call assert(this%XXX.EQ.YYY,msg='sbham default XXX is not YYY')
    !call kill(this)
    !**********************************

    !----- kill tests -----
    !***          example          ****
    !write(*,*)'test kill cleans up dynamic memory and pointers'
    !call make(this)
    !call kill(this)
    !call assert(.not.associated(this%PPP),msg='sbham pointer PPP remains associated after killed.')
    !**********************************

    !----- backup tests -----
    !***          example          ****
    !write(*,*)'test attributes are stored properly stored in backup file'
    !call make(this)
    !this%var=XXX!First, manually set sbham attributes to non-default values
    !call system('rm -f sbham.tmpfile')
    !call backup(this,file='sbham.tmpfile')
    !call kill(this)
    !call make(this,file='sbham.tmpfile')
    !call assert(this%var.EQ.XXX)!Then, assert non default attribute values are conserved
    !call kill(this)    
    !**********************************

    !----- status tests -----

    !----- update tests -----

    !----- reset tests -----
    
    !----- fail cases -----

    !========================================================================



    write(*,*)'ALL sbham TESTS PASSED!'
  end subroutine sbham_test
  !-----------------------------------------

end module sbham_class

