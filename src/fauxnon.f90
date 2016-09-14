!>\brief
!! fauxnon class
!!\details
!! A classical collection of non interacting harmonic oscillators used to
!! describe a dissipative heat bath in terms of a spectral density.
!<------------------------------------------------------------------------
module fauxnon_class
  use type_kinds
  implicit none
  private

  public::fauxnon, fauxnon_test
  public::describe, check, make, kill, backup, update, reset, status

  type fauxnon
     logical::initialized=.false.
     character(len=label)::name='fauxnon'
     
     !***************      Enter fauxnon attributes here     ***************!
     !input:
     !spectral density: list of nmode weighted unit mass frequencies
     !temperature: bath temperature
     !statistics: Bose-Einstein, Fermi-Dirac, or Maxwell-Boltzmann 

     !variables:
     !bath mode coordinates: q_k
     !bath mode momenta: p_k

     !calculate:
     !internal forces: f_k from q_k
     !correlation function: FT of specden 

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
  end type fauxnon

  !> Creates the fauxnon object.
  interface make
     module procedure fauxnon_init
  end interface

  !> Destroys the fauxnon object.
  interface kill
     module procedure fauxnon_kill
  end interface

  !> Returns current state of the fauxnon object.
  interface status
     module procedure fauxnon_status
  end interface

  !> Returns a plain text description of the fauxnon object.
  interface describe
     module procedure fauxnon_describe
  end interface
  
  !> Returns the current state of the fauxnon object.
  interface backup
     module procedure fauxnon_backup
  end interface

  !> Recaluclates the fauxnon object.
  interface update
     module procedure fauxnon_update
  end interface

  !> Reinitializes the fauxnon object.
  interface reset
     module procedure fauxnon_reset
  end interface

  !> Checks that the fauxnon object.
  interface check
     module procedure fauxnon_check
  end interface

contains
  !======================================================================
  !> \brief Retruns a description of fauxnon as a string.
  !> \param[in] THIS is the fauxnon object.
  !======================================================================
  character(len=comment) function fauxnon_describe(this)
    type(fauxnon),intent(in)::this
    character(len=5)::FMT='(A)'

    write(fauxnon_describe,FMT)'No description for fauxnon has been provided.'
   
  end function fauxnon_describe

  !======================================================================
  !> \brief Creates and initializes fauxnon.
  !! \param THIS is the fauxnon object.
  !! \param[in] FILE is an optional string containing the name of a
  !! previously backuped fauxnon file.
  !=====================================================================
  subroutine fauxnon_init(this,file)
    use filemanager
    use testing_class
    type(fauxnon),intent(inout)::this
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
       
       !check if file is of type fauxnon
       read(unit,*)header
       call assert(trim(header).EQ.this%name,msg='fauxnon_init: bad input file header in file'//file)
       
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

  end subroutine fauxnon_init

  !======================================================================
  !> \brief Destroys the fauxnon object.
  !> \param THIS is the fauxnon object to be destroyed.
  !> \remarks kill is simply the reset method passed with a null flag 
  !====================================================================
  subroutine fauxnon_kill(this)
    type(fauxnon),intent(inout)::this
 
    call reset(this,0)

  end subroutine fauxnon_kill

  !======================================================================
  !> \brief Computes the current state of fauxnon object.
  !> \param THIS is the fauxnon  object to be updated.
  !======================================================================
  subroutine fauxnon_update(this)
    type(fauxnon),intent(inout)::this

    !Recompute any attribute values that might have evolved

    !****************************      Example     ******************************
    !*** attribute 'var' is always equall to the trace of the denisity matrix *** 
    ! this%var=0._double
    ! do istate=1,this%object%nstate
    !    this%var=this%var+this%den(istate,istate)
    ! end do
    !****************************************************************************

  end subroutine fauxnon_update

  !======================================================================
  !> \brief Re-initiallizes the fauxnon object.
  !> \param THIS is the fauxnon  object to be re-initialized.
  !> \param STATE is an optional integer:
  !>        when 0, will create a null state by deallocating all dynamic
  !>        memory and returning the object to an un-initiallized state;
  !>        when not 0, will return the object to the default settings;
  !>        when not present, object will reset based on current scalar
  !>        parameters.
  !======================================================================
  subroutine fauxnon_reset(this,state)
    type(fauxnon),intent(inout)::this
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
          
          !set all static parameters to error values
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
    

  end subroutine fauxnon_reset

  !======================================================================
  !> \brief Backups the current state of the fauxnon object to file.
  !> \param[in] THIS is the fauxnon  object to be updated.
  !> \param[in] FILE is a string containing the location of the backup file.
  !======================================================================
  subroutine fauxnon_backup(this,file)
    use filemanager
    use string
    use testing_class
    type(fauxnon),intent(in)::this
    character*(*),intent(in)::file
    integer(short)::unit
    logical::fileisopen
    integer(long)::i,j
    
    !check input file
    inquire(file=file,opened=fileisopen,number=unit)
    if(unit.LT.0)unit=newunit()
    if(.not.fileisopen)open(unit,file=file)
    
    !check fauxnon object
    call assert(check(this).EQ.0,msg='fauxnon object does not pass check.')

    !always write the data type on the first line
    write(unit,*)'fauxnon'

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
  end subroutine fauxnon_backup
  
  !======================================================================
  !> \brief Retrun the current state of fauxnon as a string.
  !> \param[in] THIS is the fauxnon object.
  !> \param[in] MSG is an optional string message to annotate the status.
  !======================================================================
  character(len=line) function fauxnon_status(this,msg)
    type(fauxnon),intent(in)::this
    character*(*),intent(in),optional::msg
    character(len=5)::FMT='(A)'
    
    !Edit the status prompt to suit your needs
    write(fauxnon_status,FMT)'fauxnon status is currently not available'
    
  end function fauxnon_status

 !======================================================================
  !> \brief Checks the fauxnon object.
  !> \param[in] THIS is the fauxnon object to be checked.
  !> \return 0 if all checks pass or exit at first failed check and returm non zero.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function fauxnon_check(this)
    use testing_class
    type(fauxnon),intent(in)::this

    !initiate with no problems found 
    fauxnon_check=0

    !check that object is initialized
    call assert(this%initialized&
         ,msg='fauxnon_check: fauxnon object not initialized.'&
         ,iostat=fauxnon_check)
    if(fauxnon_check.NE.0)return

    !check that object has correct name
    call assert(this%name.EQ.'fauxnon'&
         ,msg='fauxnon_check: fauxnon name is not set.'&
         ,iostat=fauxnon_check)
    if(fauxnon_check.NE.0)return

    !Check all attributes are within acceptable values


    !**********   Example - check an object attribute 'that'  *********
    !call assert(check(this%that).EQ.0&
    !     ,msg='fauxnon_check: that sub-object failed check'&
    !     ,iostat=fauxnon_check)
    !if(fauxnon_check.NE.0)return
    !***********************************************************************
    
    !***   Example - check an integer attribute 'ndim' is well behaved   ***
    !call assert(check(this%ndim).EQ.0&
    !     ,msg='fauxnon_check: ndim failed check',iostat=fauxnon_check)
    !if(fauxnon_check.NE.0)return
    !***********************************************************************
 
    !*** Example - add a constrain that says 'ndim' can only be positive ***
    !call assert(this%ndim.GT.0&
    !     ,msg='fauxnon_check: ndim is not positive',iostat=fauxnon_check)
    !if(fauxnon_check.NE.0)return
    !***********************************************************************

    !***  Example - check a real valued attribute 'var' is well behaved  ***
    !call assert(check(this%var).EQ.0&
    !     ,msg='fauxnon_check: var failed check',iostat=fauxnon_check)
    !if(fauxnon_check.NE.0)return
    !***********************************************************************

    !***  Example - add a constrain that says 'var' can not be zero     ***
    !call assert(abs(this%var).GT.epsilon(this%var)&
    !     ,msg='fauxnon_check: var is tiny',iostat=fauxnon_check)
    !if(fauxnon_check.NE.0)return
    !***********************************************************************

    !***  Example - check a real valued pointer attribute 'matrix'       ***
    !***            is well behaved                                      ***
    !call assert(check(this%matrix).EQ.0&
    !     ,msg='fauxnon_check: matrix failed check',iostat=fauxnon_check)
    !if(fauxnon_check.NE.0)return
    !***********************************************************************

    !********* Example - check an NxM matrix has right dimensions **********
    !call assert(size(this%matrix).EQ.N*M&
    !     ,msg='fauxnon_check: number of matrix elements not = N*M.'&
    !     ,iostat=fauxnon_check)
    !if(fauxnon_check.NE.0)return
    !***********************************************************************

  end function fauxnon_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the fauxnon methods.
  !> \param[in] this is the fauxnon object whose methods will be excercised.
  !> \return Nothing if all tests pass or 1 and a stop for the first failed test.
  !> \remark Will stop after first failed check.
  !======================================================================
  subroutine fauxnon_test
    use testing_class
    use filemanager
    type(fauxnon)::this
    character(len=label)::string
    integer(long)::unit
    
    !verify fauxnon is compatible with current version
    include 'verification'

    !================== consider the following tests ========================
    !-----  make tests -----
    !***          example          ****
    !write(*,*)'test make sets correct default values'
    !call make(this)
    !call assert(this%XXX.EQ.YYY,msg='fauxnon default XXX is not YYY')
    !call kill(this)
    !**********************************

    !----- kill tests -----
    !***          example          ****
    !write(*,*)'test kill cleans up dynamic memory and pointers'
    !call make(this)
    !call kill(this)
    !call assert(.not.associated(this%PPP),msg='fauxnon pointer PPP remains associated after killed.')
    !**********************************

    !----- backup tests -----
    !***          example          ****
    !write(*,*)'test attributes are stored properly in backup file'
    !call make(this)
    !this%var=XXX!First, manually set fauxnon attributes to non-default values
    !call system('rm -f fauxnon.tmpfile')
    !call backup(this,file='fauxnon.tmpfile')
    !call kill(this)
    !call make(this,file='fauxnon.tmpfile')
    !call system('rm -f fauxnon.tmpfile')
    !call assert(this%var.EQ.XXX)!Then, assert non default attribute values are conserved
    !call kill(this)    
    !**********************************

    !----- status tests -----

    !----- update tests -----

    !----- reset tests -----
    
    !----- fail cases -----

    !========================================================================



    write(*,*)'ALL fauxnon TESTS PASSED!'
  end subroutine fauxnon_test
  !-----------------------------------------

end module fauxnon_class

