!>\brief
!! wf class
!!\details
!! Quantum object prototype specifying multiple states, in multiple dofs
!! in delta function represention  
!<------------------------------------------------------------------------
module wf_class
  use type_kinds
  implicit none
  private

  public::wf, wf_test
  public::describe, check, make, kill, backup, update, reset, status

  type wf
     logical::initialized=.false.
     character(len=label)::name='wf'
     
     integer(long)::nstate
     integer(long)::ndof
     integer(long)::npt
     !density elements in delta rep last index for position and weight
     complex(double),dimension(:,:,:,:,:),pointer::psi

     !***************      Enter wf attributes here     ***************!


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
  end type wf

  !> Creates the wf object.
  interface make
     module procedure wf_init
  end interface

  !> Destroys the wf object.
  interface kill
     module procedure wf_kill
  end interface

  !> Returns current state of the wf object.
  interface status
     module procedure wf_status
  end interface

  !> Returns a plain text description of the wf object.
  interface describe
     module procedure wf_describe
  end interface
  
  !> Returns the current state of the wf object.
  interface backup
     module procedure wf_backup
  end interface

  !> Recaluclates the wf object.
  interface update
     module procedure wf_update
  end interface

  !> Reinitializes the wf object.
  interface reset
     module procedure wf_reset
  end interface

  !> Checks that the wf object.
  interface check
     module procedure wf_check
  end interface

contains
  !======================================================================
  !> \brief Retruns a description of wf as a string.
  !> \param[in] THIS is the wf object.
  !======================================================================
  character(len=comment) function wf_describe(this)
    type(wf),intent(in)::this
    character(len=5)::FMT='(A)'

    write(wf_describe,FMT)'No description for wf has been provided.'
   
  end function wf_describe

  !======================================================================
  !> \brief Creates and initializes wf.
  !! \param THIS is the wf object.
  !! \param[in] FILE is an optional string containing the name of a
  !! previously backuped wf file.
  !=====================================================================
  subroutine wf_init(this,file)
    use filemanager
    use testing_class
    type(wf),intent(inout)::this
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
       
       !check if file is of type wf
       read(unit,*)header
       call assert(trim(header).EQ.this%name,msg='wf_init: bad input file header in file'//file)
       
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

  end subroutine wf_init

  !======================================================================
  !> \brief Destroys the wf object.
  !> \param THIS is the wf object to be destroyed.
  !> \remarks kill is simply the reset method passed with a null flag 
  !====================================================================
  subroutine wf_kill(this)
    type(wf),intent(inout)::this
 
    call reset(this,0)

  end subroutine wf_kill

  !======================================================================
  !> \brief Computes the current state of wf object.
  !> \param THIS is the wf  object to be updated.
  !======================================================================
  subroutine wf_update(this)
    type(wf),intent(inout)::this

    !Recompute calculated variables that might have evolved

    !***************************      Example     *****************************
    !** attribute 'var' is always equall to the trace of the denisity matrix ** 
    ! this%var=0._double
    ! do istate=1,this%object%nstate
    !    this%var=this%var+this%den(istate,istate)
    ! end do
    !**************************************************************************

  end subroutine wf_update

  !======================================================================
  !> \brief Re-initiallizes the wf object.
  !> \param THIS is the wf  object to be re-initialized.
  !> \param STATE is an optional integer:
  !>        when 0, will create a null state by deallocating all dynamic
  !>        memory and returning the object to an un-initiallized state;
  !>        when not 0, will return the object to the default settings;
  !>        when not present, object will reset based on current scalar
  !>        parameters.
  !======================================================================
  subroutine wf_reset(this,state)
    type(wf),intent(inout)::this
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
          
          !set all scalar parameters to error values
          !**** example **********
          !this%nstate=-1
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
    

  end subroutine wf_reset

  !======================================================================
  !> \brief Backups the current state of the wf object to file.
  !> \param[in] THIS is the wf  object to be updated.
  !> \param[in] FILE is a string containing the location of the backup file.
  !======================================================================
  subroutine wf_backup(this,file)
    use filemanager
    use string
    use testing_class
    type(wf),intent(in)::this
    character*(*),intent(in)::file
    integer(short)::unit
    logical::fileisopen
    integer(long)::i,j
    
    !check input file
    inquire(file=file,opened=fileisopen,number=unit)
    if(unit.LT.0)unit=newunit()
    if(.not.fileisopen)open(unit,file=file)
    
    !check wf object
    call assert(check(this).EQ.0,msg='wf object does not pass check.')

    !always write the data type on the first line
    write(unit,*)'wf'

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
  end subroutine wf_backup
  
  !======================================================================
  !> \brief Retrun the current state of wf as a string.
  !> \param[in] THIS is the wf object.
  !> \param[in] MSG is an optional string message to annotate the status.
  !======================================================================
  character(len=line) function wf_status(this,msg)
    type(wf),intent(in)::this
    character*(*),intent(in),optional::msg
    character(len=5)::FMT='(A)'
    
    !Edit the status prompt to suit your needs
    write(wf_status,FMT)'wf status is currently not available'
    
  end function wf_status

 !======================================================================
  !> \brief Checks the wf object.
  !> \param[in] THIS is the wf object to be checked.
  !> \return 0 if all checks pass or exit at first failed check and returm non zero.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function wf_check(this)
    use testing_class
    type(wf),intent(in)::this



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !initiate with no problems found 
    wf_check=0

    !check that object is initialized
    call assert(this%initialized&
         ,msg='wf_check: wf object not initialized.'&
         ,iostat=wf_check)
    if(wf_check.NE.0)return

    !check that object has correct name
    call assert(this%name.EQ.'wf'&
         ,msg='wf_check: wf name is not set.'&
         ,iostat=wf_check)
    if(wf_check.NE.0)return

    !Check all attributes are within acceptable values


    !**********   Example - check an object attribute 'that'  *********
    !call assert(check(this%that).EQ.0&
    !     ,msg='wf_check: that sub-object failed check'&
    !     ,iostat=wf_check)
    !if(wf_check.NE.0)return
    !***********************************************************************
    
    !***   Example - check an integer attribute 'ndim' is well behaved   ***
    !call assert(check(this%ndim).EQ.0&
    !     ,msg='wf_check: ndim failed check',iostat=wf_check)
    !if(wf_check.NE.0)return
    !***********************************************************************
 
    !*** Example - add a constrain that says 'ndim' can only be positive ***
    !call assert(this%ndim.GT.0&
    !     ,msg='wf_check: ndim is not positive',iostat=wf_check)
    !if(wf_check.NE.0)return
    !***********************************************************************

    !***  Example - check a real valued attribute 'var' is well behaved  ***
    !call assert(check(this%var).EQ.0&
    !     ,msg='wf_check: var failed check',iostat=wf_check)
    !if(wf_check.NE.0)return
    !***********************************************************************

    !***  Example - add a constrain that says 'var' can not be zero     ***
    !call assert(abs(this%var).GT.epsilon(this%var)&
    !     ,msg='wf_check: var is tiny',iostat=wf_check)
    !if(wf_check.NE.0)return
    !***********************************************************************

    !***  Example - check a real valued pointer attribute 'matrix'       ***
    !***            is well behaved                                      ***
    !call assert(check(this%matrix).EQ.0&
    !     ,msg='wf_check: matrix failed check',iostat=wf_check)
    !if(wf_check.NE.0)return
    !***********************************************************************

    !********* Example - check an NxM matrix has right dimensions **********
    !call assert(size(this%matrix).EQ.N*M&
    !     ,msg='wf_check: number of matrix elements not = N*M.'&
    !     ,iostat=wf_check)
    !if(wf_check.NE.0)return
    !***********************************************************************

  end function wf_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the wf methods.
  !> \param[in] this is the wf object whose methods will be excercised.
  !> \return Nothing if all tests pass or 1 and a stop for the first failed test.
  !> \remark Will stop after first failed check.
  !======================================================================
  subroutine wf_test
    use testing_class
    use filemanager
    type(wf)::this
    character(len=label)::string
    integer(long)::unit
    
    !verify wf is compatible with current version
    include 'verification'

    write(*,*)'test static parameter nstate is stored properly in backup file'
    call make(this) !make wf
    this%nstate=3    !manually set static parameter to non-default value
    call system('rm -f wf.tmpfile*') !remove any previous backup file(s)
    call backup(this,file='wf.tmpfile') !create backup file
    call kill(this) !destroy wf
    call make(this,file='wf.tmpfile') !make new wf from backup file
    call system('rm -f wf.tmpfile*') !remove backup file(s)
    !assert non default parameter is conserved
    call assert(this%nstate.EQ.3,&
         msg='wf static parameter nstate is not stored properly')
    call kill(this)    !destroy wf to clean up memory

    write(*,*)'test make sets correct default value for static parameter nstate'
    call make(this) !make wf
    call assert(this%nstate.EQ.1,&
         msg='wf default static parameter nstate is not 1')
    call kill(this)

    write(*,*)'test edge case for static parameter nstate breaks wf'
    call make(this) !make wf
    this%nstate=-1
    call assert(check(this).NE.0,&
         msg='edge case value -1 for static parameter nstate &
         &does not break wf')
    call kill(this)

    write(*,*)'test static parameter npt is stored properly in backup file'
    call make(this) !make wf
    this%npt=3    !manually set static parameter to non-default value
    call system('rm -f wf.tmpfile*') !remove any previous backup file(s)
    call backup(this,file='wf.tmpfile') !create backup file
    call kill(this) !destroy wf
    call make(this,file='wf.tmpfile') !make new wf from backup file
    call system('rm -f wf.tmpfile*') !remove backup file(s)
    !assert non default parameter is conserved
    call assert(this%npt.EQ.3,&
         msg='wf static parameter npt is not stored properly')
    call kill(this)    !destroy wf to clean up memory

    write(*,*)'test make sets correct default value for static parameter npt'
    call make(this) !make wf
    call assert(this%npt.EQ.1,&
         msg='wf default static parameter npt is not 1')
    call kill(this)

    write(*,*)'test edge case for static parameter npt breaks wf'
    call make(this) !make wf
    this%npt=-1
    call assert(check(this).NE.0,&
         msg='edge case value -1 for static parameter npt &
         &does not break wf')
    call kill(this)

    write(*,*)'test static parameter ndof is stored properly in backup file'
    call make(this) !make wf
    this%ndof=3    !manually set static parameter to non-default value
    call system('rm -f wf.tmpfile*') !remove any previous backup file(s)
    call backup(this,file='wf.tmpfile') !create backup file
    call kill(this) !destroy wf
    call make(this,file='wf.tmpfile') !make new wf from backup file
    call system('rm -f wf.tmpfile*') !remove backup file(s)
    !assert non default parameter is conserved
    call assert(this%ndof.EQ.3,&
         msg='wf static parameter ndof is not stored properly')
    call kill(this)    !destroy wf to clean up memory

    write(*,*)'test make sets correct default value for static parameter ndof'
    call make(this) !make wf
    call assert(this%ndof.EQ.0,&
         msg='wf default static parameter ndof is not 0')
    call kill(this)

    write(*,*)'test edge case for static parameter ndof breaks wf'
    call make(this) !make wf
    this%ndof=-1
    call assert(check(this).NE.0,&
         msg='edge case value -1 for static parameter ndof &
         &does not break wf')
    call kill(this)


    !================== consider the following tests for =====================
    !=========================================================================

    !==================       static parameters             ==================

    !!write(*,*)'test kill wf breaks static parameter NNN'
    !!call make(this) !make wf
    !!call kill(this)
    !!call assert(check(this%NNN).NE.0,&
    !!     msg='kill wf does not break static parameter NNN')
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
    !call make(this) !make wf
    !this%NNN=MMM    !manually set static parameter to non-default value
    !call system('rm -f wf.tmpfile*') !remove any previous backup file(s)
    !call backup(this,file='wf.tmpfile') !create backup file
    !call kill(this) !destroy wf
    !call make(this,file='wf.tmpfile') !make new wf from backup file
    !call system('rm -f wf.tmpfile*') !remove backup file(s)
    !!assert non default parameter is conserved
    !call assert(this%NNN.EQ.MMM,&
    !     msg='wf static parameter NNN is not stored properly')
    !call kill(this)    !destroy wf to clean up memory

    !write(*,*)'test make sets correct default value for static parameter NNN'
    !call make(this) !make wf
    !call assert(this%NNN.EQ.MMM,&
    !     msg='wf default static parameter NNN is not MMM')
    !call kill(this)

    !write(*,*)'test edge case for static parameter NNN breaks wf'
    !call make(this) !make wf
    !this%NNN=MMM
    !call assert(check(this).NE.0,&
    !     msg='edge case value MMM for static parameter NNN &
    !     &does not break wf')
    !call kill(this)

    !==================    dynamic arrays and pointers      ==================

    !write(*,*)'test make allocates memory for dyanamic pointer array PPP'
    !call make(this)
    !call assert(associated(this%PPP),msg='wf dynamic pointer array &
    !     &PPP remains associated after killed.')
    !call kill(this)

    !write(*,*)'test kill deallocates memory for dynamic pointer array PPP'
    !call make(this)
    !call kill(this)
    !call assert(.not.associated(this%PPP),msg='wf dynamic pointer array &
    !     &PPP remains associated after killed.')

    !write(*,*)'test dynamic pointer array PPP is saved properly in backup file'
    !call make(this) !make wf
    !this%PPP=YYY    !manually set dynamic pointer array to non-default value
    !call system('rm -f wf.tmpfile*') !remove any previous backup file(s)
    !call backup(this,file='wf.tmpfile') !create backup file
    !call kill(this) !destroy wf
    !call make(this,file='wf.tmpfile') !make new wf from backup file
    !call system('rm -f wf.tmpfile*') !remove backup file(s)
    !!assert non default parameter is conserved
    !call assert(all(this%PPP.EQ.YYY),&
    !     msg='wf dynamic pointer array PPP is not stored properly')
    !call kill(this)    !destroy wf to clean up memory

    !write(*,*)'test make allocates dynamic pointer array PPP of default size'
    !call make(this) !make wf
    !!assert dynamic pointer array has proper size for all dimensions
    !call assert(size(this%PPP,J).EQ.N),&
    !     msg='wf dynamic pointer array PPP is not of size N for &
    !     &dimension J')
    !call assert(size(this%PPP,I).EQ.N),&
    !     msg='wf dynamic pointer array PPP is not of size N for &
    !     &dimension I')
    !call kill(this)    !destroy wf to clean up memory
    
    !write(*,*)'test dynamic pointer array PPP can be resized by adjusting &
    !     & static parameter NNN then reseting with state=1'
    !call make(this) !make wf
    !this%NNN=MMM !adjust static parameter NNN to non-default value
    !call reset(this,state=1) !reset wf to reallocate dynamic memory
    !!assert dynamic pointer array size has changed properly
    !call assert(size(this%PPP,J).EQ.MMM),&
    !     msg='wf dynamic pointer array PPP did not change size upon &
    !     &reset with state=1')
    !call kill(this)    !destroy wf to clean up memory

    !write(*,*)'test edge case for dynamic pointer array PPP breaks wf'
    !call make(this) !make wf
    !this%PPP=YYY !set edge case value
    !assert edge case value breaks wf
    !call assert(check(this).NE.0,&
    !     msg='edge case value YYY for dynamic pointer array PPP &
    !     &does not break wf')
    !call kill(this)    !destroy wf to clean up memory

    !==================        calculated variables         ==================

    !write(*,*)'test make sets default value for calculated variable XXX'
    !call make(this) !make wf
    !call assert(this%XXX.EQ.YYY,&
    !     msg='wf default calculated variable XXX is not YYY')
    !call kill(this)

    !write(*,*)'test edge case for calculated variable XXX breaks wf'
    !call make(this) !make wf
    !this%XXX=YYY
    !call assert(check(this).NE.0,&
    !     msg='edge case value YYY for calculated variable XXX &
    !     &does not break wf')
    !call kill(this)
    !========================================================================





    write(*,*)'ALL wf TESTS PASSED!'
  end subroutine wf_test
  !-----------------------------------------

end module wf_class

