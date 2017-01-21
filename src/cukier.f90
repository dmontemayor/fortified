!>\brief
!! cukier class
!!\details
!! A proton coupled electron transfer model system as described by
!! Cukier, JPC 98, 2377 (1994).
!!
!! The system consists of one proton and one electron allowed to move in 1D between two points D and A representing a Donor and Acceptor sites respectively. Also present is a solvent represented as a boson bath with known spectral density.
!!
!! System  variables:
!! x      - the proton position
!! V(x)   - electronic coupling matrix (diabatic representation)
!!
!! System parameters:
!! V0     - electronic coupling energy scaling factor
!! gamma  - proton-electron decay distance
!! mass   - proton mass
!! Eb     - nuclear potential barrier height
!! omega0 - nuclear potential fundamental frequency in harmonic approx
!! D      - Dimensionless barrier height
!! a      - spring cutoff distance
!<------------------------------------------------------------------------
module cukier_class
  use type_kinds
  use atomicunits
  use wavelet_class

  implicit none
  private

  public::cukier, cukier_test
  public::describe, check, make, kill, backup, update, reset, status

  type cukier
     logical::initialized=.false.
     character(len=label)::name='cukier'

     !***************      Enter cukier attributes here     ***************!
     type(wavelet),private::proton !quantum mechanical proton
     integer(long)::npt      !proton dof discretization
     integer(long)::nstate !number of proton states
     real(double),dimension(:),pointer::adiabat !adiabatic state energies
     real(double)::xmin,xmax !proton dof boundaries

     !real(double)::rate
     real(double)::mass

     !type(wavelet),private::bath


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
  end type cukier

  !> Creates the cukier object.
  interface make
     module procedure cukier_init
  end interface

  !> Destroys the cukier object.
  interface kill
     module procedure cukier_kill
  end interface

  !> Returns current state of the cukier object.
  interface status
     module procedure cukier_status
  end interface

  !> Returns a plain text description of the cukier object.
  interface describe
     module procedure cukier_describe
  end interface

  !> Returns the current state of the cukier object.
  interface backup
     module procedure cukier_backup
  end interface

  !> Recaluclates the cukier object.
  interface update
     module procedure cukier_update
  end interface

  !> Reinitializes the cukier object.
  interface reset
     module procedure cukier_reset
  end interface

  !> Checks that the cukier object.
  interface check
     module procedure cukier_check
  end interface

contains
  !======================================================================
  !> \brief Retruns a description of cukier as a string.
  !> \param[in] THIS is the cukier object.
  !======================================================================
  character(len=comment) function cukier_describe(this)
    type(cukier),intent(in)::this
    character(len=5)::FMT='(A)'

    write(cukier_describe,FMT)'A proton coupled electron transfer model system as described in Cukier, JPC 98, 2377 (1994).'

  end function cukier_describe

  !======================================================================
  !> \brief Creates and initializes cukier.
  !! \param THIS is the cukier object.
  !! \param[in] FILE is an optional string containing the name of a
  !! previously backuped cukier file.
  !=====================================================================
  subroutine cukier_init(this,file)
    use filemanager
    use testing_class
    type(cukier),intent(inout)::this
    character*(*),intent(in),optional::file
    integer(long)::unit,i
    logical::fileisopen=.false.
    character(len=label)::header
    character(len=path)::infile

    !initialize all sub-objects
    !***      example     ***
    !call make(this%object)
    !************************
    call make(this%proton)

    !check input file
    if(present(file))then

       !open input file if not open and get unit number
       inquire(file=file,opened=fileisopen,number=unit)
       if(unit.LT.0)unit=newunit()
       if(.not.fileisopen)open(unit,file=file)

       !check if file is of type cukier
       read(unit,*)header
       call assert(trim(header).EQ.this%name&
            ,msg='cukier_init: bad input file header in file'//file)

       !READ static parameters
       !***          example            ***
       !***read scalar parameter from file***
       !read(unit,*)this%XXX
       !**************************************
       read(unit,*)this%npt    !proton dof discretization
       read(unit,*)this%nstate !number of proton states
       read(unit,*)this%xmin   !proton dof minval boundary
       read(unit,*)this%xmax   !proton dof maxval boundary

       !use reset to manage dynamic memory
       !, reset sub-objects, and set random parameters
       !based on current static parameters
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
       read(unit,*)infile
       call make(this%proton,file=trim(infile))

       !finished reading all attributes - now close backup file
       if(.not.fileisopen)close(unit)

       !declare initialization complete
       this%initialized=.true.

    else
       !Set static parameters to default settings
       !***       example      ***
       !***set scalar parameter***
       !this%XXX=YYY
       !**************************
       this%npt=255
       this%nstate=2
       this%xmin=-1.2_double*a0
       this%xmax=1.2_double*a0

       !Use reset to make a default object
       call reset(this,state=1)
    end if


    !finished making now update object
    call update(this)

  end subroutine cukier_init

  !======================================================================
  !> \brief Destroys the cukier object.
  !> \param THIS is the cukier object to be destroyed.
  !> use reset with state=0 flag to cleanly deallocate memory
  !====================================================================
  subroutine cukier_kill(this)
    type(cukier),intent(inout)::this

    call reset(this,state=0)

  end subroutine cukier_kill

  !======================================================================
  !> \brief Computes the current state of cukier object.
  !> \param THIS is the cukier  object to be updated.
  !======================================================================
  subroutine cukier_update(this)
    use filemanager
    use testing_class
    type(cukier),intent(inout)::this

    real(double),parameter::gamma=1.0!ion-electron interaction decay distance
    real(double),parameter::mass=mp  !mobile ion mass (H+)
    real(double),parameter::Eb=2000*invcm !double well barrier height
    real(double),parameter::w0=1200*invcm !well frequency in harmonic approx
    real(double),parameter::D=Eb/(hbar*w0)!Dimensionless barrier height
    real(double),parameter::a=sqrt(hbar/(mass*w0))!spring cutoff distance


    real(double)::V(0:this%npt-1)                   !proton potential
    real(double)::psi(0:this%nstate-1,0:this%npt-1) !proton wf

    real(double)::x,y,dx
    integer(long)::i


    !calculate proton potential
    dx=(this%xmax-this%xmin)/real(this%npt,double)
    do i=0,this%npt-1
       x=this%proton%grid(i,0)/a
       y=x*x/4
       V(i)=y*(1-y/(4*D))   !potential energy
    end do
    V=Eb-(hbar*w0)*V

    !initiate square wavefunction
    psi=1.0_double

    !initiate zero valued adiabatic energies
    this%adiabat=0.0_double

    !calculate proton states
    do i=0,this%nstate-1
       call solvestate(i,E=this%adiabat(i)&
            ,wf=psi(i,:),V=V,mass=mp,dx=dx)
    end do

    !!normalize wavefunctions
    !normalize electron wavefunction
    this%proton%wf(:,0)=psi(0,:)/sqrt(sum(psi(0,:)**2))

    !***************************      Example     *****************************
    !** attribute 'var' is always equall to the trace of the denisity matrix **
    ! this%var=0._double
    ! do istate=1,this%object%nstate
    !    this%var=this%var+this%den(istate,istate)
    ! end do
    !**************************************************************************

  end subroutine cukier_update

  !======================================================================
  !> \brief Re-initiallizes the cukier object.
  !> \param THIS is the cukier  object to be re-initialized.
  !> \param STATE is an optional integer:
  !>        when 0, will create a null state by deallocating all dynamic
  !>        memory and killing all sub-objects thus returning the object
  !>        to an un-initiallized state;
  !>        when not 0, will force reallocate dynamic memory and remake
  !>        all sub-objects based on current scalar parameters
  !>        (will create new objects thus break pointers);
  !>        when not present, object will reset based on current scalar
  !>        parameters (preserves pointers).
  !======================================================================
  subroutine cukier_reset(this,state)
    type(cukier),intent(inout)::this
    integer(long),intent(in),optional::STATE
    integer(long)::i
    real(double)::dx

    if(present(state))then
      !un-initialize cukier object until properly reset
      this%initialized=.false.
       if(state.EQ.0)then
          !nullify all pointers
          !******        Example - cleanup pointer attribute 'PPP'       ****
          !if(associated(this%PPP))nullify(this%PPP)
          !******************************************************************
          if(associated(this%adiabat))nullify(this%adiabat)

          !kill all objects
          !**** example **********
          !call kill(this%object)
          !***********************
          call kill(this%proton)

          !set static parameters to error values
          this%npt=-1
          this%nstate=-1
          this%xmin=huge(this%xmin)
          this%xmax=huge(this%xmax)

       else
         if(checkparam(this).EQ.0) then

           !reallocate dynamic memory
           !***  Example - cleanup pointer attribute 'PPP'     ***
           !***            then reallocate memory              ***
           !if(associated(this%PPP))nullify(this%PPP)
           !allocate(this%PPP(0:this%npt-1))
           !******************************************************
           if(associated(this%adiabat))nullify(this%adiabat)
           allocate(this%adiabat(0:this%nstate-1))

           !Set default dynamic memory values
           !***  Example - set values in pointer 'PPP' to zero ***
           !this%PPP(:)=0.0_double
           !******************************************************

           !overwrite sub-object default static parameters
           !***       example      ***
           !***set object static parameter***
           !this%object%XXX=123
           !**************************
           this%proton%npt=this%npt

           !reset all sub objects to correct any memory issues
           !***      example     ***
           !call reset(this%object,state=1)
           !************************
           call reset(this%proton,state=1)

           !overwrite sub-object default dynamic parameters
           !***       example      ***
           !***set object pointer array values***
           !this%object%PPP(:)=XXX
           !**************************
           !define proton grid
           !dx=0.018*angstrom
           dx=(this%xmax-this%xmin)/real(this%npt,double)
           do i=0,this%npt-1
             this%proton%grid(i,0)=(i-this%npt/2)*dx
           end do

           !declare initialization complete
           this%initialized=.true.

         end if
       end if
     end if

    !reset object based on current static parameters
    if(this%initialized)then

       !Sample Random parameters
       !***  Example - attribute 'var' samples a Gaussian random number
       ! this%var=gran()

       !Reset sub-objects
       !**** example **********
       !call reset(this%object)
       !***********************
       call reset(this%proton)

       !update now that object is fully reset and not in null state
       call update(this)
    end if

  end subroutine cukier_reset

  !======================================================================
  !> \brief Backups the current state of the cukier object to file.
  !> \param[in] THIS is the cukier  object to be updated.
  !> \param[in] FILE is a string containing the location of the backup file.
  !======================================================================
  subroutine cukier_backup(this,file)
    use filemanager
    use string
    type(cukier),intent(in)::this
    character*(*),intent(in)::file
    integer(short)::unit
    logical::fileisopen
    integer(long)::i,j

    !check input file
    inquire(file=file,opened=fileisopen,number=unit)
    if(unit.LT.0)unit=newunit()
    if(.not.fileisopen)open(unit,file=file)

    !always write the data type on the first line
    write(unit,*)'cukier'

    !******      Backup below all the derived type's attributes       ****
    !******         in the order the MAKE method reads them           ****

    !First, Static attributes
    !******          Example - Backup a scalar attribute            ******
    ! write(unit,*)this%var
    !*********************************************************************
    write(unit,*)this%npt
    write(unit,*)this%nstate
    write(unit,*)this%xmin
    write(unit,*)this%xmax


    !Second, Dynamic attributes
    !***       Example - Backup an NxM matrix attribute                ***
    ! write(unit,*)((this%matrix(i,j),j=1,M),i=1,N)
    !*********************************************************************


    !Last,objects
    !******              Example - Backup an object            ***********
    ! call backup(this%object,file//'.object')
    ! write(unit,*)quote(file//'.object')!write the object location
    !*********************************************************************
    call backup(this%proton,file//'.proton')
    write(unit,*)quote(file//'.proton')!write the object location


    !finished writing all attributes - now close backup file
    close(unit)
  end subroutine cukier_backup

  !======================================================================
  !> \brief Retrun the current state of cukier as a string.
  !> \param[in] THIS is the cukier object.
  !> \param[in] MSG is an optional string message to annotate the status.
  !======================================================================
  character(len=line) function cukier_status(this,msg)
    type(cukier),intent(in)::this
    character*(*),intent(in),optional::msg
    character(len=5)::FMT='(A)'

    !Edit the status prompt to suit your needs
    write(cukier_status,FMT)'cukier status is currently not available'

  end function cukier_status
  !======================================================================
  !> \brief Checks the cukier object parameters are good.
  !> \remark Will stop program after first failed check.
  !======================================================================
  integer(short)function checkparam(this)
    use testing_class
    type(cukier),intent(in)::this

    !initiate with no problems found
    checkparam=0

    !***   check npt is well behaved   ***
    call assert(check(this%npt).EQ.0&
         ,msg='checkparam: npt failed check',iostat=checkparam)
    if(checkparam.NE.0)return

    !***   check npt is positive   ***
    call assert(this%npt.GT.0&
         ,msg='checkparam: npt is not positive',iostat=checkparam)
    if(checkparam.NE.0)return

    !***   check nstate is well behaved   ***
    call assert(check(this%nstate).EQ.0&
         ,msg='checkparam: nstate failed check',iostat=checkparam)
    if(checkparam.NE.0)return

    !***   check nstate is positive   ***
    call assert(this%nstate.GT.0&
         ,msg='checkparam: nstate is not positive',iostat=checkparam)
    if(checkparam.NE.0)return

    !***   check xmin is well behaved   ***
    call assert(check(this%xmin).EQ.0&
          ,msg='checkparam: xmin failed check',iostat=checkparam)
    if(checkparam.NE.0)return

    !***   check xmax is well behaved   ***
    call assert(check(this%xmax).EQ.0&
         ,msg='checkparam: xmax failed check',iostat=checkparam)
    if(checkparam.NE.0)return

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
  !> \brief Checks the cukier object.
  !> \param[in] THIS is the cukier object to be checked.
  !> \return 0 if all checks pass or exit at first failed check and returm non zero.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function cukier_check(this)
    use testing_class
    type(cukier),intent(in)::this
    integer(long)::istate

    !initiate with no problems found
    cukier_check=0

    !check that object is initialized
    call assert(this%initialized&
         ,msg='cukier_check: cukier object not initialized.'&
         ,iostat=cukier_check)
    if(cukier_check.NE.0)return

    !check that object has correct name
    call assert(this%name.EQ.'cukier'&
         ,msg='cukier_check: cukier name is not set.'&
         ,iostat=cukier_check)
    if(cukier_check.NE.0)return

    !check all parameters are within acceptable values
    call assert(checkparam(this).EQ.0&
         ,msg='cukier_check: unacceptable parameters found!'&
         ,iostat=cukier_check)
    if(cukier_check.NE.0)return

    !Check all sub-objects
    !**********   Example - check an object attribute 'that'  *********
    !call assert(check(this%that).EQ.0&
    !     ,msg='cukier_check: that sub-object failed check'&
    !     ,iostat=cukier_check)
    !if(cukier_check.NE.0)return
    !***********************************************************************

    !***   check proton has one dof   ***
    call assert(this%proton%ndof.EQ.1&
         ,msg='cukier_check: proton does not have 1 dof.',iostat=cukier_check)
    if(cukier_check.NE.0)return

    !Check dynamic attributes
    !***  Example - check a real valued pointer attribute 'matrix'       ***
    !***            is well behaved                                      ***
    !call assert(check(this%matrix).EQ.0&
    !     ,msg='cukier_check: matrix failed check',iostat=cukier_check)
    !if(cukier_check.NE.0)return
    !***********************************************************************
    !********* Example - check an NxM matrix has right dimensions **********
    !call assert(size(this%matrix).EQ.N*M&
    !     ,msg='cukier_check: number of matrix elements not = N*M.'&
    !     ,iostat=cukier_check)
    !if(cukier_check.NE.0)return
    !***********************************************************************

    !***   check Adiabat is well behaved   ***
    call assert(check(this%adiabat).EQ.0&
         ,msg='cukier_check: adiabat failed check',iostat=cukier_check)
    if(cukier_check.NE.0)return

    !***   check Adiabats are ordered in ascending energy   ***
    istate=0
    do while(istate.LT.this%nstate-1.and.this%nstate.GT.1)
       call assert(all(this%adiabat(istate+1:this%nstate-1)&
            .GE.this%adiabat(istate))&
            ,msg='cukier_check: adiabats are not ordered in ascending energy',iostat=cukier_check)
       if(cukier_check.NE.0)return
       istate=istate+1
    end do

  end function cukier_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the cukier methods.
  !> \param[in] this is the cukier object whose methods will be excercised.
  !> \return Nothing if all tests pass or 1 and a stop for the first failed test.
  !> \remark Will stop after first failed check.
  !======================================================================
  subroutine cukier_test
    use testing_class
    use filemanager
    type(cukier)::this
    character(len=label)::string
    integer(long)::unit
    integer(short)::ierr

    !diagnostic vars
    integer(long)::i
    integer(long)::Rnpt
    real(double)::Rmin,Rmax,dR,R,x,y

    !verify cukier is compatible with current version
    include 'verification'


    write(*,*)'test proton is not classical...'
    call make(this)
    call assert(this%proton%npt.GT.1,msg='cukier proton is classical')
    call kill(this)

    write(*,*)'test kill cleans up proton memory'
    call make(this)
    call kill(this)
    call assert(check(this%proton).NE.0,msg='cukier proton remains viable after kill cukier.')

    write(*,*)'test negative npt breaks cukier object.'
    call make(this)
    this%npt=-1
    call assert(check(this).NE.0,msg='cukier object is not broken with negative npt')
    call kill(this)

    write(*,*)'test make sets correct default value for npt '
    call make(this)
    call assert(this%npt.EQ.255,msg='cukier default npt is not 255')
    call kill(this)

    write(*,*)'test npt attribute is stored properly in backup file'
    call make(this)
    this%npt=200!First, manually set cukier attributes to non-default values
    call reset(this,state=1)
    call system('rm -f cukier.tmpfile*')
    call backup(this,file='cukier.tmpfile')
    call kill(this)
    call make(this,file='cukier.tmpfile')
    !Then, assert non default attribute values are conserved
    call assert(this%npt.EQ.200,msg='cukier npt is not stored properly')
    call kill(this)
    call system('rm -f cukier.tmpfile*')

    write(*,*)'test negative nstate breaks cukier object.'
    call make(this)
    this%nstate=-1
    call assert(check(this).NE.0,msg='cukier object is not broken with negative ntate')
    call kill(this)

    write(*,*)'test make sets correct default value for nstate'
    call make(this)
    call assert(this%nstate.EQ.2,msg='cukier default nstate is not 2')
    call kill(this)

    write(*,*)'test nstate attribute is stored properly in backup file'
    call make(this)
    this%nstate=4!First, manually set cukier attributes to non-default values
    call reset(this,state=1)
    call system('rm -f cukier.tmpfile*')
    call backup(this,file='cukier.tmpfile')
    call kill(this)
    call make(this,file='cukier.tmpfile')
    !Then, assert non default attribute values are conserved
    call assert(this%nstate.EQ.4,msg='cukier nstate is not stored properly')
    call kill(this)
    call system('rm -f cukier.tmpfile*')

    write(*,*)'test make calculates normalized proton wf.'
    call make(this)
    call assert(real(sum(this%proton%wf(:,0)*conjg(this%proton%wf(:,0))))&
         ,1._double,tol=1E-8_double&
         ,msg='proton wave function is not normalized.')
    call kill(this)

    write(*,*)'test cuckier object breaks when 1st excited state is lower energy than ground state'
    call make(this)
    this%adiabat(0)=1.0
    this%adiabat(1)=0.0
    call assert(check(this).NE.0,msg='cukier object did not break when 1st excited state is lower energy than ground state.')

    write(*,*)'test make generates degenerate ground and 1st excited states for default settings.'
    call make(this)
    x=this%adiabat(0)
    y=this%adiabat(1)
    call assert((x-y).EQ.0.0_double.or.&
         abs(x-y)/y.LT.0.03_double,msg='ground and 1st excited states are not degenerate.')
    call kill(this)

    write(*,*)'test make sets correct default proton dof boundaries.'
    call make(this)
    call assert(this%xmin.EQ.-1.2_double*a0,msg='default xmin is not -1.2 a0.')
    call assert(this%xmax.EQ.1.2_double*a0,msg='default xmax is not 1.2 a0.')
    call kill(this)

    write(*,*)'test xmin attribute is stored properly in backup file'
    call make(this)
    this%xmin=-2.0_double!First, manually set cukier attributes to non-default values
    call reset(this,state=1)
    call system('rm -f cukier.tmpfile*')
    call backup(this,file='cukier.tmpfile')
    call kill(this)
    call make(this,file='cukier.tmpfile')
    !Then, assert non default attribute values are conserved
    call assert(this%xmin.EQ.-2.0_double,&
         msg='cukier xmin is not stored properly')
    call kill(this)
    call system('rm -f cukier.tmpfile*')

    write(*,*)'test xmax attribute is stored properly in backup file'
    call make(this)
    this%xmax=2.0_double!First, manually set cukier attributes to non-default values
    call reset(this,state=1)
    call system('rm -f cukier.tmpfile*')
    call backup(this,file='cukier.tmpfile')
    call kill(this)
    call make(this,file='cukier.tmpfile')
    !Then, assert non default attribute values are conserved
    call assert(this%xmax.EQ.2.0_double,&
         msg='cukier xmax is not stored properly')
    call kill(this)
    call system('rm -f cukier.tmpfile*')

    write(*,*)'test cukier model break when proton has more than 1 dof.'
    call make(this)
    this%proton%ndof=2
    call reset(this,state=1)
    call assert(check(this).NE.0,msg='cukier model does not break when&
         & proton has more than 1 dof.')
    call kill(this)

    write(*,*)'test make generates non-degenerate ground and 2st excited states for default settings.'
    call make(this)
    this%nstate=3
    call reset(this,state=1)
    x=this%adiabat(0)
    y=this%adiabat(2)
    call assert(abs(x-y).GT.0.1*x,msg='ground and 2st excited states are degenerate.')
    call kill(this)

    !write(*,*)'test model predicts correct isotope effect on PCET rate.'
    !call make(this)
    !this%mass=mp !H+
    !!this%D(0)=1.0_double !donor definition
    !!this%D(1)=1.0_double
    !!this%A(0)=1.0_double !acceptor definition
    !!this%A(1)=-1.0_double
    !call reset(this,state=1)
    !x=this%rate !kH
    !call kill(this)
    !call make(this)
    !this%mass=mp+mn !D+
    !!this%D(0)=1.0_double !donor definition
    !!this%D(1)=1.0_double
    !!this%A(0)=1.0_double !acceptor definition
    !!this%A(1)=-1.0_double
    !call reset(this,state=1)
    !y=this%rate !kD
    !call kill(this)
    !call assert(x/y,1.4_double,0.4_double,&
  !       msg='isotope effect kH/kD is not between 1.0-1.8')

    !================== consider the following tests ========================
    !-----  make tests -----
    !***          example          ****
    !write(*,*)'test make sets correct default values'
    !call make(this)
    !call assert(this%XXX.EQ.YYY,msg='cukier default XXX is not YYY')
    !call kill(this)
    !**********************************

    !----- kill tests -----
    !***          example          ****
    !write(*,*)'test kill cleans up dynamic memory and pointers'
    !call make(this)
    !call kill(this)
    !call assert(.not.associated(this%PPP),msg='cukier pointer PPP remains associated after killed.')
    !**********************************

    !----- backup tests -----
    !***          example          ****
    !write(*,*)'test attributes are stored properly stored in backup file'
    !call make(this)
    !this%var=XXX!First, manually set cukier attributes to non-default values
    !call system('rm -f cukier.tmpfile*')
    !call backup(this,file='cukier.tmpfile')
    !call kill(this)
    !call make(this,file='cukier.tmpfile')
    !call assert(this%var.EQ.XXX)!Then, assert non default attribute values are conserved
    !call system('rm -f cukier.tmpfile*')
    !call kill(this)
    !**********************************

    !----- status tests -----

    !----- update tests -----

    !----- reset tests -----

    !----- fail cases -----

    !========================================================================



    write(*,*)'ALL cukier TESTS PASSED!'
  end subroutine cukier_test
  !-----------------------------------------

end module cukier_class
