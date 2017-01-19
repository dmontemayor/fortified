!!!!!!!! silcrow ATTRIBUTE LIST !!!!!!!!
!§ mass double static 0
!§ E_b double static 0
!§ E_p double static 0
!§ omega_r double static 0
!§ omega_b double static 0
!§ omega_c double static 0
!left off here need to make silcrow object factory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>\brief
!! doublewell class
!!\details
!! A 1D asymmetric double harmonic well potential with inverted
!! harmonic barrier.
!!
!! The barrier potential (V_b) is centered at the origin q=0
!! with frequency omega_b and energy shift= Eb.
!! The reactant potenial (V_r) is centered at q=q_r<0
!! with frequency omega_r and no energy shift.
!! The product potenial (V_p) is centered at q=q_p>0
!! with frequency omega_p and with energy shift= Ep.
!!
!! The potentials are written
!! V_r(q)= 1/2 * mass * omega_r^2 * (q-q_r)^2
!! V_b(q)=-1/2 * mass * omega_b^2 * q^2 + Eb
!! V_p(q)= 1/2 * mass * omega_r^2 * (q-q_p)^2 +Ep
!!
!! The barrier potential seamlessly connects to the reactant potenial
!! at the point q_{r,b} and to the product potential at the point
!! q_{p,b} such that
!! V_r(q_{r,b}) = V_b(q_{r,b})
!! and
!! V_p(q_{p,b}) = V_b(q_{p,b})
!!
!! In addition,
!! The double well potential is differentiable at all points so
!! dV_r(q_{r,b})/dq = dV_b(q_{r,b})/dq
!! and
!! dV_p(q_{p,b})/dq = dV_b(q_{p,b})/dq
!!
!! which yields:
!! q_{r,b} = q_r * omega_r^2 / (omega_r^2 + omega_b^2)
!! and
!! q_{p,b} = q_p * omega_p^2 / (omega_p^2 + omega_b^2)
!! for the intersection points
!! with
!! q_r = sqrt(2/mass * Eb * (omega_r^2 + omega_b^2)) / (omega_r * omega_b)
!! and
!! q_p = sqrt(2/mass * (Eb-Ep) * (omega_p^2 + omega_b^2)) / (omega_p * omega_b)
!! for the well centers
!!
!! System can thus be described by the 5 parameters:
!! omega_r, omega_b, omega_p, Eb, and Ep
!<------------------------------------------------------------------------
module doublewell_class
  use type_kinds
  use atomicunits
  implicit none
  private

  public::doublewell, doublewell_test
  public::describe, check, make, kill, backup, update, reset, status

  type doublewell
     logical::initialized=.false.
     character(len=label)::name='doublewell'
     
     !***************      Enter doublewell attributes here     ***************!
     !integer(long),private::npt    !dof discretization
     real(double),private::mass    !mass of reaction
     real(double),private::E_b     !barrier potential energy height
     real(double),private::E_p     !product side potential energy difference
     real(double),private::omega_r !reactant side potential energy frequency
     real(double),private::omega_b !barrier potential energy frequency
     real(double),private::omega_p !product side potential energy frequency

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
  end type doublewell

  !> Creates the doublewell object.
  interface make
     module procedure doublewell_init
  end interface

  !> Destroys the doublewell object.
  interface kill
     module procedure doublewell_kill
  end interface

  !> Returns current state of the doublewell object.
  interface status
     module procedure doublewell_status
  end interface

  !> Returns a plain text description of the doublewell object.
  interface describe
     module procedure doublewell_describe
  end interface
  
  !> Returns the current state of the doublewell object.
  interface backup
     module procedure doublewell_backup
  end interface

  !> Recaluclates the doublewell object.
  interface update
     module procedure doublewell_update
  end interface

  !> Reinitializes the doublewell object.
  interface reset
     module procedure doublewell_reset
  end interface

  !> Checks the doublewell object.
  interface check
     module procedure doublewell_check
  end interface

contains
  !> \brief Get function for reactant side minimum
  real(double) function get_qr(this)
    type(doublewell),intent(in)::this
    get_qr=-sqrt(2.0_double*this%E_b*(this%omega_r**2+this%omega_b**2)/this%mass)&
         /(this%omega_r*this%omega_b)
  end function get_qr
  !====================  
  !> \brief Get function for product side minimum
  real(double) function get_qp(this)
    type(doublewell),intent(in)::this
    get_qp=sqrt(2.0_double*(this%E_b-this%E_p)*(this%omega_p**2+this%omega_b**2)/this%mass)&
         /(this%omega_p*this%omega_b)
  end function get_qp
  !====================  
  !> \brief Get function for reactant-barrier intersection
  real(double) function get_qrb(this)
    type(doublewell),intent(in)::this
    get_qrb=get_qr(this)*this%omega_r**2/(this%omega_r**2+this%omega_b**2)
  end function get_qrb
  !====================  
  !> \brief Get function for product-barrier intersection
  real(double) function get_qpb(this)
    type(doublewell),intent(in)::this
    get_qpb=get_qp(this)*this%omega_p**2/(this%omega_p**2+this%omega_b**2)
  end function get_qpb
  !====================  
  !======================================================================
  !> \brief Retruns a description of doublewell as a string.
  !> \param[in] THIS is the doublewell object.
  !======================================================================
  character(len=comment) function doublewell_describe(this)
    type(doublewell),intent(in)::this
    character(len=5)::FMT='(A)'

    write(doublewell_describe,FMT)'No description for doublewell has been provided.'
   
  end function doublewell_describe

  !======================================================================
  !> \brief Creates and initializes doublewell.
  !! \param THIS is the doublewell object.
  !! \param[in] FILE is an optional string containing the name of a
  !! previously backuped doublewell file.
  !=====================================================================
  subroutine doublewell_init(this,file)
    use filemanager
    use testing_class
    type(doublewell),intent(inout)::this
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
       
       !check if file is of type doublewell
       read(unit,*)header
       call assert(trim(header).EQ.this%name,msg='doublewell_init: bad input file header in file'//file)
       
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
       !this%npt=500
       this%mass=mp
       this%E_b=2000_double*invcm
       this%E_p=0_double*invcm
       this%omega_r=1200_double*invcm
       this%omega_b=1200_double*invcm
       this%omega_p=1200_double*invcm
       !***       example      ***
       !***set scalar parameter***
       !this%NNN=123
       !**************************
       
       !Use reset to make a default object
       call reset(this,state=1)
    end if


    !finished making now update object
    call update(this)

  end subroutine doublewell_init

  !======================================================================
  !> \brief Destroys the doublewell object.
  !> \param THIS is the doublewell object to be destroyed.
  !> \remarks kill is simply the reset method passed with a null flag 
  !====================================================================
  subroutine doublewell_kill(this)
    type(doublewell),intent(inout)::this
 
    call reset(this,0)

  end subroutine doublewell_kill

  !======================================================================
  !> \brief Computes the current state of doublewell object.
  !> \param THIS is the doublewell  object to be updated.
  !======================================================================
  subroutine doublewell_update(this)
    type(doublewell),intent(inout)::this

    !Recompute calculated variables that might have evolved

    !***************************      Example     *****************************
    !** attribute 'var' is always equall to the trace of the denisity matrix ** 
    ! this%var=0._double
    ! do istate=1,this%object%nstate
    !    this%var=this%var+this%den(istate,istate)
    ! end do
    !**************************************************************************

  end subroutine doublewell_update

  !======================================================================
  !> \brief Re-initiallizes the doublewell object.
  !> \param THIS is the doublewell  object to be re-initialized.
  !> \param STATE is an optional integer:
  !>        when 0, will create a null state by deallocating all dynamic
  !>        memory and returning the object to an un-initiallized state;
  !>        when not 0, will return the object to the default settings;
  !>        when not present, object will reset based on current scalar
  !>        parameters.
  !======================================================================
  subroutine doublewell_reset(this,state)
    type(doublewell),intent(inout)::this
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
          !this%mass=-1.0_double
          !this%E_b=huge(1_double)
          !this%E_p=huge(1_double)
          !this%omega_r=-1.0_double
          !this%omega_b=-1.0_double
          !this%omega_p=-1.0_double
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
    

  end subroutine doublewell_reset

  !======================================================================
  !> \brief Backups the current state of the doublewell object to file.
  !> \param[in] THIS is the doublewell  object to be updated.
  !> \param[in] FILE is a string containing the location of the backup file.
  !======================================================================
  subroutine doublewell_backup(this,file)
    use filemanager
    use string
    use testing_class
    type(doublewell),intent(in)::this
    character*(*),intent(in)::file
    integer(short)::unit
    logical::fileisopen
    integer(long)::i,j
    
    !check input file
    inquire(file=file,opened=fileisopen,number=unit)
    if(unit.LT.0)unit=newunit()
    if(.not.fileisopen)open(unit,file=file)
    
    !check doublewell object
    call assert(check(this).EQ.0,msg='doublewell object does not pass check.')

    !always write the data type on the first line
    write(unit,*)'doublewell'

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
  end subroutine doublewell_backup
  
  !======================================================================
  !> \brief Retrun the current state of doublewell as a string.
  !> \param[in] THIS is the doublewell object.
  !> \param[in] MSG is an optional string message to annotate the status.
  !======================================================================
  character(len=line) function doublewell_status(this,msg)
    type(doublewell),intent(in)::this
    character*(*),intent(in),optional::msg
    character(len=5)::FMT='(A)'
    
    !Edit the status prompt to suit your needs
    write(doublewell_status,FMT)'doublewell status is currently not available'
    
  end function doublewell_status

 !======================================================================
  !> \brief Checks the doublewell object.
  !> \param[in] THIS is the doublewell object to be checked.
  !> \return 0 if all checks pass or exit at first failed check and returm non zero.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function doublewell_check(this)
    use testing_class
    type(doublewell),intent(in)::this
    


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
    doublewell_check=0

    !check that object is initialized
    call assert(this%initialized&
         ,msg='doublewell_check: doublewell object not initialized.'&
         ,iostat=doublewell_check)
    if(doublewell_check.NE.0)return

    !check that object has correct name
    call assert(this%name.EQ.'doublewell'&
         ,msg='doublewell_check: doublewell name is not set.'&
         ,iostat=doublewell_check)
    if(doublewell_check.NE.0)return

    !Check all attributes are within acceptable values

    !**********   Example - check an object attribute 'that'  *********
    !call assert(check(this%that).EQ.0&
    !     ,msg='doublewell_check: that sub-object failed check'&
    !     ,iostat=doublewell_check)
    !if(doublewell_check.NE.0)return
    !***********************************************************************
    
    !***   Example - check an integer attribute 'ndim' is well behaved   ***
    !call assert(check(this%ndim).EQ.0&
    !     ,msg='doublewell_check: ndim failed check',iostat=doublewell_check)
    !if(doublewell_check.NE.0)return
    !***********************************************************************
 
    !*** Example - add a constrain that says 'ndim' can only be positive ***
    !call assert(this%ndim.GT.0&
    !     ,msg='doublewell_check: ndim is not positive',iostat=doublewell_check)
    !if(doublewell_check.NE.0)return
    !***********************************************************************

    !***  Example - check a real valued attribute 'var' is well behaved  ***
    !call assert(check(this%var).EQ.0&
    !     ,msg='doublewell_check: var failed check',iostat=doublewell_check)
    !if(doublewell_check.NE.0)return
    !***********************************************************************

    !***  Example - add a constrain that says 'var' can not be zero     ***
    !call assert(abs(this%var).GT.epsilon(this%var)&
    !     ,msg='doublewell_check: var is tiny',iostat=doublewell_check)
    !if(doublewell_check.NE.0)return
    !***********************************************************************

    !***  Example - check a real valued pointer attribute 'matrix'       ***
    !***            is well behaved                                      ***
    !call assert(check(this%matrix).EQ.0&
    !     ,msg='doublewell_check: matrix failed check',iostat=doublewell_check)
    !if(doublewell_check.NE.0)return
    !***********************************************************************

    !********* Example - check an NxM matrix has right dimensions **********
    !call assert(size(this%matrix).EQ.N*M&
    !     ,msg='doublewell_check: number of matrix elements not = N*M.'&
    !     ,iostat=doublewell_check)
    !if(doublewell_check.NE.0)return
    !***********************************************************************


!!$    !well behaved real attribute mass
!!$    
!!$    call assert(check(this%ndim).EQ.0&
!!$         ,msg='doublewell_check: ndim failed check',iostat=doublewell_check)
!!$    if(doublewell_check.NE.0)return






    
  end function doublewell_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the doublewell methods.
  !> \param[in] this is the doublewell object whose methods will be excercised.
  !> \return Nothing if all tests pass or 1 and a stop for the first failed test.
  !> \remark Will stop after first failed check.
  !======================================================================
  subroutine doublewell_test
    use testing_class
    use filemanager
    type(doublewell)::this
    character(len=label)::string
    integer(long)::unit
    
    !verify doublewell is compatible with current version
    include 'verification'

    !======Unit testing=====
    ! Testing of individual components. Typically done by
    ! the developer and not by testers as it requires
    ! detailed knowledge of the internal program design.

!!$    write(*,*)'test kill doublewell breaks static parameter mass.....'
!!$    call make(this) !make doublewell
!!$    call kill(this)
!!$    call assert(check(this%mass).NE.0,&
!!$         msg='kill doublewell does not break static parameter mass')
!!$
!!$    write(*,*)'test kill doublewell breaks static parameter E_b.....'
!!$    call make(this) !make doublewell
!!$    call kill(this)
!!$    call assert(check(this%E_b).NE.0,&
!!$         msg='kill doublewell does not break static parameter E_b')
!!$
!!$    write(*,*)'test kill doublewell breaks static parameter E_p.....'
!!$    call make(this) !make doublewell
!!$    call kill(this)
!!$    call assert(check(this%E_p).NE.0,&
!!$         msg='kill doublewell does not break static parameter E_p')
!!$
!!$    write(*,*)'test kill doublewell breaks static parameter omega_r.....'
!!$    call make(this) !make doublewell
!!$    call kill(this)
!!$    call assert(check(this%omega_r).NE.0,&
!!$         msg='kill doublewell does not break static parameter omega_r')
!!$
!!$    write(*,*)'test kill doublewell breaks static parameter omega_b.....'
!!$    call make(this) !make doublewell
!!$    call kill(this)
!!$    call assert(check(this%omega_b).NE.0,&
!!$         msg='kill doublewell does not break static parameter omega_b')
!!$
!!$    write(*,*)'test kill doublewell breaks static parameter omega_p.....'
!!$    call make(this) !make doublewell
!!$    call kill(this)
!!$    call assert(check(this%omega_p).NE.0,&
!!$         msg='kill doublewell does not break static parameter omega_p')

        
    !======Functional testing======
    ! This type of testing igores the internal parts
    ! and focuses on the output as per requirements.

    write(*,*)'test default potential is produced.....'
    ! mass=mp, Eb=2000*invcm, Ep=0, and w0=1200*invcm where omega_r=omega_b=omega_p=w0
    ! resulting in -qr=qp=0.8140 bohr and -qrb=qpb=qp/2=0.4070 bohr
    call make(this)
    call assert(get_qr(this),-0.8140_double,.0001_double&
         ,msg='default reactant side minimum is not -0.8140 bohr')
    call assert(get_qp(this),0.8140_double,.0001_double&
         ,msg='default reactant side minimum is not 0.8140 bohr')
    call assert(get_qrb(this),get_qr(this)/2_double,.0001_double&
         ,msg='default reactant-barrier intersection is not -0.4070 bohr')
    call assert(get_qpb(this),get_qp(this)/2_double,.0001_double&
         ,msg='default product-barrier intersection is not 0.4070 bohr')

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

    !!write(*,*)'test kill doublewell breaks static parameter NNN.....'
    !!call make(this) !make doublewell
    !!call kill(this)
    !!call assert(check(this%NNN).NE.0,&
    !!     msg='kill doublewell does not break static parameter NNN')
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
    !call make(this) !make doublewell
    !this%NNN=MMM    !manually set static parameter to non-default value
    !call system('rm -f doublewell.tmpfile*') !remove any previous backup file(s)
    !call backup(this,file='doublewell.tmpfile') !create backup file
    !call kill(this) !destroy doublewell
    !call make(this,file='doublewell.tmpfile') !make new doublewell from backup file
    !call system('rm -f doublewell.tmpfile*') !remove backup file(s)
    !!assert non default parameter is conserved
    !call assert(this%NNN.EQ.MMM,&
    !     msg='doublewell static parameter NNN is not stored properly')
    !call kill(this)    !destroy doublewell to clean up memory

    !write(*,*)'test make sets correct default value for static parameter NNN'
    !call make(this) !make doublewell
    !call assert(this%NNN.EQ.MMM,&
    !     msg='doublewell default static parameter NNN is not MMM')
    !call kill(this)

    !write(*,*)'test edge case for static parameter NNN breaks doublewell'
    !call make(this) !make doublewell
    !this%NNN=MMM
    !call assert(check(this).NE.0,&
    !     msg='edge case value MMM for static parameter NNN &
    !     &does not break doublewell')
    !call kill(this)

    !==================    dynamic arrays and pointers      ==================

    !write(*,*)'test make allocates memory for dyanamic pointer array PPP'
    !call make(this)
    !call assert(associated(this%PPP),msg='doublewell dynamic pointer array &
    !     &PPP remains associated after killed.')
    !call kill(this)

    !write(*,*)'test kill deallocates memory for dynamic pointer array PPP'
    !call make(this)
    !call kill(this)
    !call assert(.not.associated(this%PPP),msg='doublewell dynamic pointer array &
    !     &PPP remains associated after killed.')

    !write(*,*)'test dynamic pointer array PPP is saved properly in backup file'
    !call make(this) !make doublewell
    !this%PPP=YYY    !manually set dynamic pointer array to non-default value
    !call system('rm -f doublewell.tmpfile*') !remove any previous backup file(s)
    !call backup(this,file='doublewell.tmpfile') !create backup file
    !call kill(this) !destroy doublewell
    !call make(this,file='doublewell.tmpfile') !make new doublewell from backup file
    !call system('rm -f doublewell.tmpfile*') !remove backup file(s)
    !!assert non default parameter is conserved
    !call assert(all(this%PPP.EQ.YYY),&
    !     msg='doublewell dynamic pointer array PPP is not stored properly')
    !call kill(this)    !destroy doublewell to clean up memory

    !write(*,*)'test make allocates dynamic pointer array PPP of default size'
    !call make(this) !make doublewell
    !!assert dynamic pointer array has proper size for all dimensions
    !call assert(size(this%PPP,J).EQ.N),&
    !     msg='doublewell dynamic pointer array PPP is not of size N for &
    !     &dimension J')
    !call assert(size(this%PPP,I).EQ.N),&
    !     msg='doublewell dynamic pointer array PPP is not of size N for &
    !     &dimension I')
    !call kill(this)    !destroy doublewell to clean up memory
    
    !write(*,*)'test dynamic pointer array PPP can be resized by adjusting &
    !     & static parameter NNN then reseting with state=1'
    !call make(this) !make doublewell
    !this%NNN=MMM !adjust static parameter NNN to non-default value
    !call reset(this,state=1) !reset doublewell to reallocate dynamic memory
    !!assert dynamic pointer array size has changed properly
    !call assert(size(this%PPP,J).EQ.MMM),&
    !     msg='doublewell dynamic pointer array PPP did not change size upon &
    !     &reset with state=1')
    !call kill(this)    !destroy doublewell to clean up memory

    !write(*,*)'test edge case for dynamic pointer array PPP breaks doublewell'
    !call make(this) !make doublewell
    !this%PPP=YYY !set edge case value
    !assert edge case value breaks doublewell
    !call assert(check(this).NE.0,&
    !     msg='edge case value YYY for dynamic pointer array PPP &
    !     &does not break doublewell')
    !call kill(this)    !destroy doublewell to clean up memory

    !==================        calculated variables         ==================

    !write(*,*)'test make sets default value for calculated variable XXX'
    !call make(this) !make doublewell
    !call assert(this%XXX.EQ.YYY,&
    !     msg='doublewell default calculated variable XXX is not YYY')
    !call kill(this)

    !write(*,*)'test edge case for calculated variable XXX breaks doublewell'
    !call make(this) !make doublewell
    !this%XXX=YYY
    !call assert(check(this).NE.0,&
    !     msg='edge case value YYY for calculated variable XXX &
    !     &does not break doublewell')
    !call kill(this)
    !========================================================================





    write(*,*)'ALL doublewell TESTS PASSED!'
  end subroutine doublewell_test
  !-----------------------------------------

end module doublewell_class

