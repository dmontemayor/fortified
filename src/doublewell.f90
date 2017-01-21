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
  use wavelet_class
  implicit none
  private

  public::doublewell, doublewell_test
  public::describe, check, make, kill, backup, update, reset, status

  type doublewell
     logical::initialized=.false.
     character(len=label)::name='doublewell'

     !***************      Enter doublewell attributes here     ***************!
     real(double)::E_b     !barrier potential energy height
     real(double)::E_p     !product side potential energy difference
     real(double)::omega_r !reactant side potential energy frequency
     real(double)::omega_b !barrier potential energy frequency
     real(double)::omega_p !product side potential energy frequency
     real(double)::mass    !mass of reaction
     real(double),dimension(:),pointer::V       !reaction potential
     real(double)::xmin,xmax     !reaction coordinate boundaries
     integer(long)::npt          !reaction coordinate discretization
     integer(long)::nstate       !number of quantum states
     real(double),dimension(:),pointer::adiabat !quantum state energies
     real(double),dimension(:,:),pointer::psi   !quantum state wavefunctions
     type(wavelet)::rxn          !reaction wf in delta function representation
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
    integer(long)::i,j

    !initialize all sub-objects
    !***      example     ***
    !call make(this%object)
    !************************
    call make(this%rxn)

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
       read(unit,*)this%E_b     !barrier potential energy height
       read(unit,*)this%E_p     !product side potential energy difference
       read(unit,*)this%omega_r !reactant side potential energy frequency
       read(unit,*)this%omega_b !barrier potential energy frequency
       read(unit,*)this%omega_p !product side potential energy frequency
       read(unit,*)this%mass    !mass of reaction
       read(unit,*)this%xmin    !reaction coordinate boundaries
       read(unit,*)this%xmax    !reaction coordinate boundaries
       read(unit,*)this%npt     !reaction coordinate discretization
       read(unit,*)this%nstate  !number of quantum states

       !use reset to manage dynamic memory
       !, reset sub-objects, and set calculated variable seeds
       call reset(this,state=1)

       !READ dynamic array values
       !***      example     ***
       !read(unit,*)(this%PPP(i),i=0,this%NNN-1)
       !************************
       !read(unit,*)(this%V(i),i=0,this%npt-1)          !reaction potential
       !read(unit,*)(this%adiabat(i),i=0,this%nstate-1) !rxn Eigenstate energies
       !read(unit,*)((this%psi(i,j),j=0,this%npt-1),i=0,this%nstate-1) !rxn wf

       !READ sub-objects
       !***      example     ***
       !read(unit,*)infile
       !call make(this%object,file=trim(infile))
       !************************
       !read(unit,*)infile
       !call make(this%rxn,file=trim(infile)) !reaction eignenstate wavefunctions

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
       this%E_b=2000_double*invcm     !barrier potential energy height
       this%E_p=0_double              !product side potential energy difference
       this%omega_r=1200_double*invcm !reactant side potential energy frequency
       this%omega_b=1200_double*invcm !barrier potential energy frequency
       this%omega_p=1200_double*invcm !product side potential energy frequency
       this%mass=mp                   !mass of reaction
       this%xmin=-1.2_double*a0       !reaction coordinate maxiumum boundary
       this%xmax=1.2_double*a0        !reaction coordinate maxiumum boundary
       this%npt=255                   !reaction coordinate discretization
       this%nstate=1                  !number of quantum states

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

    real(double)::V(0:this%npt-1)                   !reaction potential
    real(double)::psi(0:this%nstate-1,0:this%npt-1) !reaction wf
    real(double)::qrb,qpb,qr,qp
    real(double)::omega,x0,shift,konst
    real(double)::x,y,dx
    integer(long)::i

    !get well extrema and intersection points
    qr=get_qr(this)
    qp=get_qp(this)
    qrb=get_qrb(this)
    qpb=get_qpb(this)

    !calculate reaction potential
    dx=(this%xmax-this%xmin)/real(this%npt,double)
    do i=0,this%npt-1
      x=this%rxn%grid(i,0)
      if(x.LT.qrb)then
        !reactant potential
        x0=qr
        omega=this%omega_r
        shift=0_double
        konst=0.5_double
      else
        if(x.LT.qpb)then
          !barrier potential
          x0=0_double
          omega=this%omega_b
          shift=this%E_b
          konst=-0.5_double
        else
          !product potential
          x0=qp
          omega=this%omega_p
          shift=this%E_p
          konst=0.5_double
        end if
      end if
      !piecewise potential energy
      V(i)=konst*this%mass*(omega*(x-x0))**2+shift
    end do

    !initiate square wavefunction
    psi=1.0_double

    !initiate zero valued adiabatic energies
    this%adiabat=0.0_double

    !calculate reaction quantum states
    do i=0,this%nstate-1
       call solvestate(i,E=this%adiabat(i)&
            ,wf=psi(i,:),V=V,mass=this%mass,dx=dx)
            !normalize wavefuntion
            psi(i,:)=psi(i,:)/sqrt(sum(psi(i,:)**2))
    end do

    !return normalized wavefunction and reaction potential
    this%rxn%wf(:,0)=psi(0,:)
    this%V=V
    this%psi=0.0_double
    if(all(psi.EQ.psi))this%psi=psi

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
    integer(long)::i
    real(double)::dx

    if(present(state))then
      !un-initialize doublewell object until properly reset
      this%initialized=.false.
       if(state.EQ.0)then
          !nullify all pointers
          !******        Example - cleanup pointer attribute 'PPP'       ****
          !if(associated(this%PPP))nullify(this%PPP)
          !******************************************************************
          if(associated(this%V))nullify(this%V)
          if(associated(this%adiabat))nullify(this%adiabat)
          if(associated(this%psi))nullify(this%psi)

          !kill all sub-objects
          !**** example **********
          !call kill(this%object)
          !***********************
          call kill(this%rxn)

          !set all scalar parameters to error values
          !**** example **********
          !this%nstate=-1
          !***********************
          this%E_b=huge(1_double)     !barrier potential energy height
          this%E_p=huge(1_double)     !product side potential energy difference
          this%omega_r=-1_double      !reactant side potential energy frequency
          this%omega_b=-1_double      !barrier potential energy frequency
          this%omega_p=-1_double      !product side potential energy frequency
          this%mass=-1_double         !mass of reaction
          this%xmin=huge(1_double)    !reaction coordinate min boundary
          this%xmax=huge(1_double)    !reaction coordinate max boundary
          this%npt=-1                 !reaction coordinate discretization
          this%nstate=-1              !number of quantum states
       else
         if(checkparam(this).EQ.0) then
           !reallocate dynamic memory
           !***  Example - cleanup pointer attribute 'PPP'     ***
           !***            then reallocate memory              ***
           !if(associated(this%PPP))nullify(this%PPP)
           !allocate(this%PPP(0:this%npt-1))
           !******************************************************
           if(associated(this%V))nullify(this%V)
           allocate(this%V(0:this%npt-1))
           if(associated(this%adiabat))nullify(this%adiabat)
           allocate(this%adiabat(0:this%nstate-1))
           if(associated(this%psi))nullify(this%psi)
           allocate(this%psi(0:this%nstate-1,0:this%npt-1))

           !Set default dynamic memory values
           !***  Example - set values in pointer 'PPP' to zero ***
           !this%PPP(:)=0.0_double
           !******************************************************


           !overwrite sub-object default static parameters
           !***       example      ***
           !***set object static parameter***
           !this%object%NNN=123
           !**************************
           this%rxn%npt=this%npt

           !reset all sub objects to correct any memory issues
           !***      example     ***
           !call reset(this%object,state=1)
           !************************
           call reset(this%rxn,state=1)

           !overwrite sub-object default dynamic array
           !***       example      ***
           !***set object pointer array values***
           !this%object%PPP(:)=XXX
           !**************************
           !define reaction coordinate grid
           !dx=0.018*angstrom
           dx=(this%xmax-this%xmin)/real(this%npt,double)
           do i=0,this%npt-1
             this%rxn%grid(i,0)=this%xmin+i*dx
           end do

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
       call reset(this%rxn)

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
    write(unit,*)this%E_b     !barrier potential energy height
    write(unit,*)this%E_p     !product side potential energy difference
    write(unit,*)this%omega_r !reactant side potential energy frequency
    write(unit,*)this%omega_b !barrier potential energy frequency
    write(unit,*)this%omega_p !product side potential energy frequency
    write(unit,*)this%mass    !mass of reaction
    write(unit,*)this%xmin    !reaction coordinate max boundary
    write(unit,*)this%xmax    !reaction coordinate max boundary
    write(unit,*)this%npt     !reaction coordinate discretization
    write(unit,*)this%nstate  !number of quantum states


    !Second, Dynamic arrays
    !***       Example - Backup an NxM matrix                          ***
    ! write(unit,*)((this%matrix(i,j),j=1,M),i=1,N)
    !*********************************************************************
    !write(unit,*)(this%V(i),i=0,this%npt-1)          !reaction potential
    !write(unit,*)(this%adiabat(i),i=0,this%nstate-1) !rxn Eigenstate energies
    !write(unit,*)((this%psi(i,j),j=0,this%npt-1),i=0,this%nstate-1) !rxn wfs


    !Last,objects
    !******              Example - Backup an object            ***********
    ! call backup(this%object,file//'.object')
    ! write(unit,*)quote(file//'.object')!write the object location
    !*********************************************************************
    !call backup(this%rxn,file//'.rxn')  !backup the reaction object
    !write(unit,*)quote(file//'.rxn')    !write the object location


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
  !> \brief Checks the doublewell object parameters are good.
  !> \remark Will stop program after first failed check.
  !======================================================================
  integer(short)function checkparam(this)
    use testing_class
    type(doublewell),intent(in)::this

    !initiate with no problems found
    checkparam=0

    !barrier potential energy height
    !is well behaved
    call assert(check(this%E_b).EQ.0&
         ,msg='checkparam: E_b failed check',iostat=checkparam)
    if(checkparam.NE.0)return

    !product side potential energy difference
    !is well behaved
    call assert(check(this%E_p).EQ.0&
         ,msg='checkparam: E_p failed check',iostat=checkparam)
    if(checkparam.NE.0)return

    !reactant side potential energy frequency
    !is well behaved
    call assert(check(this%omega_r).EQ.0&
         ,msg='checkparam: omega_r failed check',iostat=checkparam)
    if(checkparam.NE.0)return
    !is positive
    call assert(abs(this%omega_r).GT.0_double&
         ,msg='checkparam: omega_r is not positive',iostat=checkparam)
    if(checkparam.NE.0)return

    !barrier potential energy frequency
    !is well behaved
    call assert(check(this%omega_b).EQ.0&
         ,msg='checkparam: omega_b failed check',iostat=checkparam)
    if(checkparam.NE.0)return
    !is positive
    call assert(abs(this%omega_b).GT.0_double&
         ,msg='checkparam: omega_b is not positive',iostat=checkparam)
    if(checkparam.NE.0)return

    !product side potential energy frequency
    !is well behaved
    call assert(check(this%omega_p).EQ.0&
         ,msg='checkparam: omega_p failed check',iostat=checkparam)
    if(checkparam.NE.0)return
    !is positive
    call assert(abs(this%omega_p).GT.0_double&
         ,msg='checkparam: omega_p is not positive',iostat=checkparam)
    if(checkparam.NE.0)return

    !mass of reaction
    !is well behaved
    call assert(check(this%mass).EQ.0&
         ,msg='checkparam: mass failed check',iostat=checkparam)
    if(checkparam.NE.0)return
    !is positive
    call assert(abs(this%mass).GT.0_double&
         ,msg='checkparam: mass is not positive',iostat=checkparam)
    if(checkparam.NE.0)return

    !reaction coordinate min boundary
    !is well behaved
    call assert(check(this%xmin).EQ.0&
         ,msg='checkparam: xmin failed check',iostat=checkparam)
    if(checkparam.NE.0)return

    !reaction coordinate max boundary
    !is well behaved
    call assert(check(this%xmax).EQ.0&
         ,msg='checkparam: xmax failed check',iostat=checkparam)
    if(checkparam.NE.0)return

    !reaction coordinate discretization
    !is well behaved
    call assert(check(this%npt).EQ.0&
    ,msg='checkparam: npt failed check',iostat=checkparam)
    if(checkparam.NE.0)return
    !is positive
    call assert(this%npt.GT.0&
    ,msg='checkparam: npt is not positive',iostat=checkparam)
    if(checkparam.NE.0)return

    !number of quantum states
    !is well behaved
    call assert(check(this%nstate).EQ.0&
    ,msg='checkparam: nstate failed check',iostat=checkparam)
    if(checkparam.NE.0)return
    !is positive
    call assert(this%nstate.GT.0&
    ,msg='checkparam: nstate is not positive',iostat=checkparam)
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
  !> \brief Checks the doublewell object.
  !> \param[in] THIS is the doublewell object to be checked.
  !> \return 0 if all checks pass or exit at first failed check and returm non zero.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function doublewell_check(this)
    use testing_class
    type(doublewell),intent(in)::this
    integer(long)::istate

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

    !check all parameters are within acceptable values
    call assert(checkparam(this).EQ.0&
         ,msg='doublewell_check: unacceptable parameters found!'&
         ,iostat=doublewell_check)
    if(doublewell_check.NE.0)return

    !Check all sub-objects
    !**********   Example - check an object attribute 'that'  *********
    !call assert(check(this%that).EQ.0&
    !     ,msg='doublewell_check: that sub-object failed check'&
    !     ,iostat=doublewell_check)
    !if(doublewell_check.NE.0)return
    !***********************************************************************
    !***   check rxn has one dof   ***
    call assert(this%rxn%ndof.EQ.1&
         ,msg='doublewell_check: rxn does not have 1 dof.',iostat=doublewell_check)
    if(doublewell_check.NE.0)return

    !Check dynamic attributes
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
    !***   check Potential is well behaved   ***
    call assert(check(this%V).EQ.0&
         ,msg='doublewell_check: V failed check',iostat=doublewell_check)
    if(doublewell_check.NE.0)return

    !***   check Potential is size npt  ***
    call assert(size(this%V).EQ.this%npt&
         ,msg='doublewell_check: V is not size npt',iostat=doublewell_check)
    if(doublewell_check.NE.0)return

    !***   check Adiabat is well behaved   ***
    call assert(check(this%adiabat).EQ.0&
         ,msg='doublewell_check: adiabat failed check',iostat=doublewell_check)
    if(doublewell_check.NE.0)return

    !***   check Adiabat dim 1 is size nstate  ***
    call assert(size(this%adiabat,1).EQ.this%nstate&
         ,msg='doublewell_check: adiabat dim 1 is not size nstate'&
        ,iostat=doublewell_check)
    if(doublewell_check.NE.0)return

    !***   check Adiabats are ordered in ascending energy   ***
    istate=0
    do while(istate.LT.this%nstate-1.and.this%nstate.GT.1)
       call assert(all(this%adiabat(istate+1:this%nstate-1)&
            .GE.this%adiabat(istate))&
            ,msg='doublewell_check: adiabats are not in ascending order'&
            ,iostat=doublewell_check)
       if(doublewell_check.NE.0)return
       istate=istate+1
    end do

    !***   check psi is well behaved   ***
    call assert(check(this%psi).EQ.0&
         ,msg='doublewell_check: psi failed check',iostat=doublewell_check)
    if(doublewell_check.NE.0)return

    !***   check psi dim 1 is size nstate  ***
    call assert(size(this%psi,1).EQ.this%nstate&
         ,msg='doublewell_check: psi dim 1 is not size nstate'&
        ,iostat=doublewell_check)
    if(doublewell_check.NE.0)return

    !***   check psi dim 2 is size npt ***
    call assert(size(this%psi,2).EQ.this%npt&
         ,msg='doublewell_check: psi dim 2 is not size npt'&
        ,iostat=doublewell_check)
    if(doublewell_check.NE.0)return



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
    real(double)::norm,tol

    !verify doublewell is compatible with current version
    include 'verification'

    !======Unit testing=====
    ! Testing of individual components. Typically done by
    ! the developer and not by testers as it requires
    ! detailed knowledge of the internal program design.


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
    call assert(this%V(this%npt/2),2000_double*invcm,1_double*invcm&
         ,msg='default potential is not 2000 invcm at origin')
    call kill(this)

    write(*,*)'test wavefunction is generated for default potential'
    call make(this)
    call assert(real(this%rxn%wf(0,0),double).LT.&
        real(this%rxn%wf(1,0),double)&
        ,msg='ground state wavfunction not generated')
    call assert(real(sum(this%rxn%wf(:,0)*conjg(this%rxn%wf(:,0))),double)&
        ,1.0_double,(this%xmax-this%xmin)/real(this%npt,double)&
        ,msg='ground state not normalized')
    call kill(this)
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
