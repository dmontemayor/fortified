!>\brief
!! metiu class
!!\details
!! A charge transfer model system as described by Shin and Metiu, JCP 102, 9285 (1995).
!!
!! The system consists of three ions and one electron. Ions A and B are not
!! allowed to move however ion C and the electron can move along the line
!! connecting the two fixed ions. 
!!
!! System variables include:
!! the electron position x
!! the moving ion position R
!<------------------------------------------------------------------------
module metiu_class
  use type_kinds
  use atomicunits
  use wavelet_class

  implicit none
  private

  public::metiu, metiu_test
  public::describe, check, make, kill, backup, update, reset, status

  type metiu
     logical::initialized=.false.
     character(len=label)::name='metiu'
     
     !***************      Enter metiu attributes here     ***************!
     type(wavelet),private::ionA,ionB,ionC,electron !all matter

     !model parameters
     integer(long)::npt    !electron dof discretization
     integer(long)::nstate !number of electronic states
     real(double)::R       !ionC position                                      
     real(double)::V_b     !barrier on ground adiabatic state                  
     real(double)::R_c     !ion-electron interaction potential cut-off distance

     real(double),dimension(:),pointer::Eigenvalue !adiabatic state energies

     
!!$    !scan potential
!!$    integer(long),parameter::npt=256,nstate=2
!!$    real(double),parameter::Rmin=-4.99*Angstrom,Rmax=4.99*Angstrom
!!$    real(double),parameter::dR=(Rmax-Rmin)/real(npt)
!!$    real(double),parameter::R1=-5.0_double*Angstrom  !ion A position
!!$    real(double),parameter::R2= 5.0_double*Angstrom  !ion B position
!!$    real(double)::R  !ion C position
!!$    real(double)::x  !electron position
!!$    real(double)::y  !electron-ion distance
!!$    real(double)::V_N !Coulombic interaction of ion C with fixed ions
!!$    real(double)::V_e(0:npt-1)!electronic potential parameterized by R
!!$    real(double)::E(0:nstate-1,0:npt-1)  !adiabatic energies parameterized by R
!!$    real(double)::psi(0:nstate-1,0:npt-1) !adiabatic wave functions for given R


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
  end type metiu

  !> Creates the metiu object.
  interface make
     module procedure metiu_init
  end interface

  !> Destroys the metiu object.
  interface kill
     module procedure metiu_kill
  end interface

  !> Returns current state of the metiu object.
  interface status
     module procedure metiu_status
  end interface

  !> Returns a plain text description of the metiu object.
  interface describe
     module procedure metiu_describe
  end interface
  
  !> Returns the current state of the metiu object.
  interface backup
     module procedure metiu_backup
  end interface

  !> Recaluclates the metiu object.
  interface update
     module procedure metiu_update
  end interface

  !> Reinitializes the metiu object.
  interface reset
     module procedure metiu_reset
  end interface

  !> Checks that the metiu object.
  interface check
     module procedure metiu_check
  end interface

contains
  !======================================================================
  !> \brief Retruns a description of metiu as a string.
  !> \param[in] THIS is the metiu object.
  !======================================================================
  character(len=comment) function metiu_describe(this)
    type(metiu),intent(in)::this
    character(len=5)::FMT='(A)'

    write(metiu_describe,FMT)'A charge transfer model system as described by Shin and Metiu, JCP 102, 9285 (1995)'
   
  end function metiu_describe

  !======================================================================
  !> \brief Creates and initializes metiu.
  !! \param THIS is the metiu object.
  !! \param[in] FILE is an optional string containing the name of a
  !! previously backuped metiu file.
  !=====================================================================
  subroutine metiu_init(this,file)
    use filemanager
    use testing_class
    type(metiu),intent(inout)::this
    character*(*),intent(in),optional::file
    integer(long)::unit,i
    logical::fileisopen=.false.
    character(len=label)::header
    character(len=path)::infile

    !initialize all sub-objects
    !***      example     ***
    !call make(this%object)
    !************************
    call make(this%ionA)
    call make(this%ionB)
    call make(this%ionC)
    call make(this%electron)
    
    !check input file
    if(present(file))then 
       
       !open input file if not open and get unit number
       inquire(file=file,opened=fileisopen,number=unit)
       if(unit.LT.0)unit=newunit()
       if(.not.fileisopen)open(unit,file=file)
       
       !check if file is of type metiu
       read(unit,*)header
       call assert(trim(header).EQ.this%name,msg='metiu_init: bad input file header in file'//file)
              
       !READ static parameters
       !***          example            ***
       !***read scalar parameter from file***
       !read(unit,*)this%XXX
       !**************************************
       read(unit,*)this%npt !electron dof discretization
       read(unit,*)this%nstate !number of electronic states
       read(unit,*)this%R   !ionC position
       read(unit,*)this%R_c !ion-electron interaction potential cut-off distance

       
       !use reset to manage dynamic memory, reset sub-objects, and set random parameters
       !based on current static parameters
       call reset(this,state=1)

       !READ dynamic array values
       !***      example     ***
       !read(unit,*)(this%PPP(i),i=0,this%XXX-1)
       !************************
       !read(unit,*)(this%Eigenvalue(i),i=0,this%nstate-1)
    
       !READ sub-objects
       !***      example     ***
       !read(unit,*)infile
       !call make(this%object,file=trim(infile))
       !************************
       read(unit,*)infile
       call make(this%ionA,file=trim(infile))
       read(unit,*)infile
       call make(this%ionB,file=trim(infile))
       read(unit,*)infile
       call make(this%ionC,file=trim(infile))
       read(unit,*)infile
       call make(this%electron,file=trim(infile))

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
       this%npt=128
       this%nstate=2
       this%R=0.0_double
       this%R_c=1.5_double*angstrom
       
       !Use reset to make a default object
       call reset(this,state=1)
    end if

    

  end subroutine metiu_init

  !======================================================================
  !> \brief Destroys the metiu object.
  !> \param THIS is the metiu object to be destroyed.
  !> use reset with state=0 flag to cleanly deallocate memory
  !====================================================================
  subroutine metiu_kill(this)
    type(metiu),intent(inout)::this
 
    call reset(this,state=0)

  end subroutine metiu_kill

  !======================================================================
  !> \brief Computes the current state of metiu object.
  !> \param THIS is the metiu  object to be updated.
  !======================================================================
  subroutine metiu_update(this)
    use filemanager
    use testing_class
    type(metiu),intent(inout)::this

    !scan potential
    !integer(long),parameter::npt=256,nstate=2
    !real(double),parameter::Rmin=-4.99*Angstrom,Rmax=4.99*Angstrom
    !real(double),parameter::dR=(Rmax-Rmin)/real(this%npt)
    !real(double),parameter::R1=-5.0_double*Angstrom  !ion A position
    !real(double),parameter::R2= 5.0_double*Angstrom  !ion B position
    real(double)::R1 !ion A position
    real(double)::R2 !ion B position
    real(double)::Rmin !ion A position
    real(double)::Rmax !ion B position
    real(double)::dR !dof increment
    real(double)::R  !ion C position
    real(double)::x  !electron position
    real(double)::y  !electron-ion distance
    real(double)::V_N !Coulombic interaction of ion C with fixed ions
    real(double)::V_e(0:this%npt-1)!electronic potential parameterized by R
    real(double)::E(0:this%nstate-1,0:this%npt-1)  !adiabatic energies parameterized by R
    real(double)::psi(0:this%nstate-1,0:this%npt-1) !adiabatic wave functions for given R
    integer(long)::unit,i,j,unit2

    R1=this%ionA%grid(0,0)
    R2=this%ionB%grid(0,0)
    Rmin=this%ionA%grid(0,0)
    Rmax=this%ionB%grid(0,0)
    dR=(Rmax-Rmin)/real(this%npt)



    unit=newunit()
    open(unit,file='metiu_pot.dat')
    open(unit2,file='metiu_eng.dat')
    do i=0,this%npt-1
       R=Rmin+i*dR
       !ion-ion interaction
       V_N=kC/abs(R-R1)+kC/abs(R-R2)

       !calculate electron-ion interaction
       V_e=0.0_double
       do j=0,this%npt-1
          x=Rmin+j*dR
          !interaction with ion A
          y=abs(x-R1)
          if(y.GT.epsilon(y))&
               V_e(j)=V_e(j)-erf(y/this%R_c)/y
          !interaction with ion B
          y=abs(x-R2)
          if(y.GT.epsilon(y))&
               V_e(j)=V_e(j)-erf(y/this%R_c)/y
          !interaciton with ion C
          y=abs(x-R)
          if(y.GT.epsilon(y))&
               V_e(j)=V_e(j)-erf(y/this%R_c)/y
       end do
       V_e=V_e+V_N !add ion-ion interaction
       call solvestate(0,E=E(0,i),wf=psi(0,:),V=V_e,mass=me,dx=dR)
       call solvestate(1,E=E(1,i),wf=psi(1,:),V=V_e,mass=me,dx=dR)

       do j=0,this%npt-1
          x=Rmin+j*dR
          write(unit,*)R/Angstrom,x/Angstrom,V_e(j)/eV,psi(0,j),psi(1,j)
       end do
       write(unit,*)
       write(unit2,*)R/Angstrom,E(0,i)/eV,E(1,i)/eV,V_N/eV
       !write(*,*)i,E(0,i)/eV,E(1,i)/eV
    end do
    close(unit)
    close(unit2)

!left off here

    !calculate V_b
    write(*,*)E(0,this%npt/2)/eV,minval(E(0,:))/eV,maxval(E(0,this%npt/4:this%npt*3/4))/eV
    this%V_b=maxval(E(0,this%npt/4:this%npt*3/4))-minval(E(0,:))

    !!calculate Delta
    !this%Delta=E(1,this%npt/2)-E(0,this%npt/2)

    !Recompute any attribute values that might have evolved



    !****************************      Example     ******************************
    !*** attribute 'var' is always equall to the trace of the denisity matrix *** 
    ! this%var=0._double
    ! do istate=1,this%object%nstate
    !    this%var=this%var+this%den(istate,istate)
    ! end do
    !****************************************************************************

  end subroutine metiu_update

  !======================================================================
  !> \brief Re-initiallizes the metiu object.
  !> \param THIS is the metiu  object to be re-initialized.
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
  subroutine metiu_reset(this,state)
    type(metiu),intent(inout)::this
    integer(long),intent(in),optional::STATE
    
    if(present(state))then
       if(state.EQ.0)then
          !nullify all pointers
          !******        Example - cleanup pointer attribute 'PPP'       ****
          !if(associated(this%PPP))nullify(this%PPP)
          !******************************************************************
          if(associated(this%Eigenvalue))nullify(this%Eigenvalue)
          
          !kill all objects
          !**** example **********
          !call kill(this%object)
          !***********************
          call kill(this%ionA)
          call kill(this%ionB)
          call kill(this%ionC)
          call kill(this%electron)
          
          !un-initialized metiu object
          this%initialized=.false.
       else

          !allocate dynamic memory
          !***  Example - cleanup pointer attribute 'PPP'     ***
          !***            then reallocate memory              ***
          !if(associated(this%PPP))nullify(this%PPP)
          !allocate(this%PPP(0:this%npt-1))
          !******************************************************
          if(associated(this%Eigenvalue))nullify(this%Eigenvalue)
          allocate(this%Eigenvalue(0:this%nstate-1))

          !Set default dynamic memory values
          !***  Example - set values in pointer 'PPP' to zero ***
          !this%PPP(:)=0.0_double
          !******************************************************
                    
          !overwrite sub-object default static parameters
          !***       example      ***
          !***set object static parameter***
          !this%object%XXX=123
          !**************************
          this%electron%npt=this%npt
          
          !reset all sub objects to correct any memory issues
          !***      example     ***
          !call reset(this%object,state=1)
          !************************
          call reset(this%ionA,state=1)
          call reset(this%ionB,state=1)
          call reset(this%ionC,state=1)
          call reset(this%electron,state=1)
          
          !overwrite sub-object default dynamic parameters
          !***       example      ***
          !***set object pointer array values***
          !this%object%PPP(:)=XXX
          !**************************
          this%ionA%grid=-5.0_double*angstrom
          this%ionB%grid=5.0_double*angstrom
          this%ionC%grid=this%R

          !declare initialization complete
          this%initialized=.true.
          
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
       call reset(this%ionA)
       call reset(this%ionB)
       call reset(this%ionC)
       call reset(this%electron)
       
    end if
    
  end subroutine metiu_reset

  !======================================================================
  !> \brief Backups the current state of the metiu object to file.
  !> \param[in] THIS is the metiu  object to be updated.
  !> \param[in] FILE is a string containing the location of the backup file.
  !======================================================================
  subroutine metiu_backup(this,file)
    use filemanager
    use string
    type(metiu),intent(in)::this
    character*(*),intent(in)::file
    integer(short)::unit
    logical::fileisopen
    integer(long)::i,j
    
    !check input file
    inquire(file=file,opened=fileisopen,number=unit)
    if(unit.LT.0)unit=newunit()
    if(.not.fileisopen)open(unit,file=file)
    
    !always write the data type on the first line
    write(unit,*)'metiu'

    !******      Backup below all the derived type's attributes       ****
    !******         in the order the MAKE method reads them           ****

    !First, Static attributes
    !******          Example - Backup a scalar attribute            ******
    ! write(unit,*)this%var
    !*********************************************************************
    write(unit,*)this%npt
    write(unit,*)this%nstate
    write(unit,*)this%R   
    write(unit,*)this%R_c 

    
    !Second, Dynamic attributes
    !***       Example - Backup an NxM matrix attribute                ***
    ! write(unit,*)((this%matrix(i,j),j=1,M),i=1,N)
    !*********************************************************************


    !Last,objects
    !******              Example - Backup an object            ***********
    ! call backup(this%object,file//'.object')
    ! write(unit,*)quote(file//'.object')!write the object location
    !*********************************************************************
    call backup(this%ionA,file//'.ionA')
    write(unit,*)quote(file//'.ionA')!write the object location
    call backup(this%ionB,file//'.ionB')
    write(unit,*)quote(file//'.ionB')!write the object location
    call backup(this%ionC,file//'.ionC')
    write(unit,*)quote(file//'.ionC')!write the object location
    call backup(this%electron,file//'.electron')
    write(unit,*)quote(file//'.electron')!write the object location

     
    !finished writing all attributes - now close backup file
    close(unit)  
  end subroutine metiu_backup
  
  !======================================================================
  !> \brief Retrun the current state of metiu as a string.
  !> \param[in] THIS is the metiu object.
  !> \param[in] MSG is an optional string message to annotate the status.
  !======================================================================
  character(len=line) function metiu_status(this,msg)
    type(metiu),intent(in)::this
    character*(*),intent(in),optional::msg
    character(len=5)::FMT='(A)'
    
    !Edit the status prompt to suit your needs
    write(metiu_status,FMT)'metiu status is currently not available'
    
  end function metiu_status

 !======================================================================
  !> \brief Checks the metiu object.
  !> \param[in] THIS is the metiu object to be checked.
  !> \return 0 if all checks pass or exit at first failed check and returm non zero.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function metiu_check(this)
    use testing_class
    type(metiu),intent(in)::this

    !initiate with no problems found 
    metiu_check=0

    !check that object is initialized
    call assert(this%initialized&
         ,msg='metiu_check: metiu object not initialized.'&
         ,iostat=metiu_check)
    if(metiu_check.NE.0)return

    !check that object has correct name
    call assert(this%name.EQ.'metiu'&
         ,msg='metiu_check: metiu name is not set.'&
         ,iostat=metiu_check)
    if(metiu_check.NE.0)return

    !Check all attributes are within acceptable values

    !***   check npt is well behaved   ***
    call assert(check(this%npt).EQ.0&
         ,msg='metiu_check: npt failed check',iostat=metiu_check)
    if(metiu_check.NE.0)return

    !***   check nstate is well behaved   ***
    call assert(check(this%nstate).EQ.0&
         ,msg='metiu_check: nstate failed check',iostat=metiu_check)
    if(metiu_check.NE.0)return

    !***   check nstate is positive   ***
    call assert(this%nstate.GT.0&
         ,msg='metiu_check: nstate is not positive',iostat=metiu_check)
    if(metiu_check.NE.0)return

    !***   check R is well behaved   ***
    call assert(check(this%R).EQ.0&
         ,msg='metiu_check: R failed check',iostat=metiu_check)
    if(metiu_check.NE.0)return

    !***   check R is between -5 and 5 angstrom   ***
    call assert(this%R.GT.-5*angstrom .and. this%R.LT.5*angstrom&
         ,msg='metiu_check: R not in range (-5:5) angstrom',iostat=metiu_check)
    if(metiu_check.NE.0)return

    !***   check R_c is well behaved   ***
    call assert(check(this%R_c).EQ.0&
         ,msg='metiu_check: R_c failed check',iostat=metiu_check)
    if(metiu_check.NE.0)return

    !***   check R_C is positive   ***
    call assert(this%R_c.GT.epsilon(this%R_c)&
         ,msg='metiu_check: R_c is not positive',iostat=metiu_check)
    if(metiu_check.NE.0)return
    
    !***   check Eigenvalue is well behaved   ***
    call assert(check(this%Eigenvalue).EQ.0&
         ,msg='metiu_check: Eigenvalue failed check',iostat=metiu_check)
    if(metiu_check.NE.0)return

    !**********   Example - check an object attribute 'primitive'  *********
    !call assert(check(this%primitive).EQ.0&
    !     ,msg='metiu_check: primitive sub-object failed check'&
    !     ,iostat=metiu_check)
    !if(metiu_check.NE.0)return
    !***********************************************************************
    
    !***   Example - check an integer attribute 'ndim' is well behaved   ***
    !call assert(check(this%ndim).EQ.0&
    !     ,msg='metiu_check: ndim failed check',iostat=metiu_check)
    !if(metiu_check.NE.0)return
    !***********************************************************************
 
    !*** Example - add a constrain that says 'ndim' can only be positive ***
    !call assert(this%ndim.GT.0&
    !     ,msg='metiu_check: ndim is not positive',iostat=metiu_check)
    !if(metiu_check.NE.0)return
    !***********************************************************************

    !***  Example - check a real valued attribute 'var' is well behaved  ***
    !call assert(check(this%var).EQ.0&
    !     ,msg='metiu_check: var failed check',iostat=metiu_check)
    !if(metiu_check.NE.0)return
    !***********************************************************************

    !***  Example - add a constrain that says 'var' can not be zero     ***
    !call assert(abs(this%var).GT.epsilon(this%var)&
    !     ,msg='metiu_check: var is tiny',iostat=metiu_check)
    !if(metiu_check.NE.0)return
    !***********************************************************************

    !***  Example - check a real valued pointer attribute 'matrix'       ***
    !***            is well behaved                                      ***
    !call assert(check(this%matrix).EQ.0&
    !     ,msg='metiu_check: matrix failed check',iostat=metiu_check)
    !if(metiu_check.NE.0)return
    !***********************************************************************

    !********* Example - check an NxM matrix has right dimensions **********
    !call assert(size(this%matrix).EQ.N*M&
    !     ,msg='metiu_check: number of matrix elements not = N*M.'&
    !     ,iostat=metiu_check)
    !if(metiu_check.NE.0)return
    !***********************************************************************

  end function metiu_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the metiu methods.
  !> \param[in] this is the metiu object whose methods will be excercised.
  !> \return Nothing if all tests pass or 1 and a stop for the first failed test.
  !> \remark Will stop after first failed check.
  !======================================================================
  subroutine metiu_test
    use testing_class
    use filemanager
    type(metiu)::this
    character(len=label)::string
    integer(long)::unit
    
    !verify metiu is compatible with current version
    include 'verification'

    write(*,*)'test ion A is classical...'
    call make(this)
    call assert(this%ionA%npt.EQ.1,msg='metiu ion A is not classical')
    call kill(this)
    write(*,*)'test kill cleans up ion A memory'
    call make(this)
    call kill(this)
    call assert(check(this%ionA).NE.0,msg='metiu ionA remains viable after kill metiu.')

    write(*,*)'test ion B is classical...'
    call make(this)
    call assert(this%ionB%npt.EQ.1,msg='metiu ion B is not classical')
    call kill(this)
    write(*,*)'test kill cleans up ion B memory'
    call make(this)
    call kill(this)
    call assert(check(this%ionB).NE.0,msg='metiu ionB remains viable after kill metiu.')

    write(*,*)'test ion C is classical...'
    call make(this)
    call assert(this%ionC%npt.EQ.1,msg='metiu ion C is not classical')
    call kill(this)
    write(*,*)'test kill cleans up ion C memory'
    call make(this)
    call kill(this)
    call assert(check(this%ionC).NE.0,msg='metiu ionC remains viable after kill metiu.')

    write(*,*)'test electron is not classical...'
    call make(this)
    call assert(this%electron%npt.GT.1,msg='metiu electron is classical')
    call kill(this)
    write(*,*)'test kill cleans up electron memory'
    call make(this)
    call kill(this)
    call assert(check(this%electron).NE.0,msg='metiu electron remains viable after kill metiu.')

    write(*,*)'test make sets correct default value for npt '
    call make(this)
    call assert(this%npt.EQ.128,msg='metiu default npt is not 128')
    call kill(this)

    write(*,*)'test npt attribute is stored properly in backup file'
    call make(this)
    this%npt=200!First, manually set metiu attributes to non-default values
    call reset(this,state=1)
    call system('rm -f metiu.tmpfile')
    call backup(this,file='metiu.tmpfile')
    call kill(this)
    call make(this,file='metiu.tmpfile')
    !Then, assert non default attribute values are conserved
    call assert(this%npt.EQ.200,msg='metiu npt is not stored properly')
    call kill(this)    
    call system('rm -f metiu.tmpfile')

    write(*,*)'test make sets correct default position for ionA '
    call make(this)
    call assert(this%ionA%grid(0,0),-5*angstrom&
         ,msg='metiu default ionA position is not -5 angstrom')
    call kill(this)

    write(*,*)'test make sets correct default density for ionA '
    call make(this)
    call assert(real(this%ionA%wf(0)*conjg(this%ionA%wf(0))),1.0_double&
         ,msg='metiu default ionA density is not 1')
    call kill(this)

    write(*,*)'test ionA attribute is stored properly in backup file'
    call make(this)
    this%ionA%npt=2!First, manually set metiu attributes to non-default values
    call reset(this,state=1)
    call system('rm -f metiu.tmpfile')
    call backup(this,file='metiu.tmpfile')
    call kill(this)
    call make(this,file='metiu.tmpfile')
    !Then, assert non default attribute values are conserved
    call assert(this%ionA%npt.EQ.2,msg='metiu ionA is not stored properly')
    call kill(this)    
    call system('rm -f metiu.tmpfile')

    write(*,*)'test make sets correct default position for ionB '
    call make(this)
    call assert(this%ionB%grid(0,0),5*angstrom&
         ,msg='metiu default ionB position is not 5 angstrom')
    call kill(this)

    write(*,*)'test make sets correct default density for ionB '
    call make(this)
    call assert(real(this%ionB%wf(0)*conjg(this%ionB%wf(0))),1.0_double&
         ,msg='metiu default ionB density is not 1')
    call kill(this)

    write(*,*)'test ionB attribute is stored properly in backup file'
    call make(this)
    this%ionB%npt=2!First, manually set metiu attributes to non-default values
    call reset(this,state=1)
    call system('rm -f metiu.tmpfile')
    call backup(this,file='metiu.tmpfile')
    call kill(this)
    call make(this,file='metiu.tmpfile')
    !Then, assert non default attribute values are conserved
    call assert(this%ionB%npt.EQ.2,msg='metiu ionB is not stored properly')
    call kill(this)    
    call system('rm -f metiu.tmpfile')

    write(*,*)'test make sets correct default position for ionC '
    call make(this)
    call assert(this%ionC%grid(0,0),0*angstrom&
         ,msg='metiu default ionC position is not 0 angstrom')
    call kill(this)

    write(*,*)'test make sets correct default density for ionC '
    call make(this)
    call assert(real(this%ionC%wf(0)*conjg(this%ionC%wf(0))),1.0_double&
         ,msg='metiu default ionC density is not 1')
    call kill(this)

    write(*,*)'test ionC attribute is stored properly in backup file'
    call make(this)
    this%ionC%npt=2!First, manually set metiu attributes to non-default values
    call reset(this,state=1)
    call system('rm -f metiu.tmpfile')
    call backup(this,file='metiu.tmpfile')
    call kill(this)
    call make(this,file='metiu.tmpfile')
    !Then, assert non default attribute values are conserved
    call assert(this%ionC%npt.EQ.2,msg='metiu ionC is not stored properly')
    call kill(this)    
    call system('rm -f metiu.tmpfile')

    write(*,*)'test electron attribute is stored properly in backup file'
    call make(this)
    this%npt=200!First, manually set metiu attributes to non-default values
    call reset(this,state=1)
    call system('rm -f metiu.tmpfile')
    call backup(this,file='metiu.tmpfile')
    call kill(this)
    call make(this,file='metiu.tmpfile')
    !Then, assert non default attribute values are conserved
    call assert(this%electron%npt.EQ.200,msg='metiu electron is not stored properly')
    call kill(this)
    call system('rm -f metiu.tmpfile')

    write(*,*)'test make sets correct default value for R '
    call make(this)
    call assert(this%R,0.0_double,msg='metiu default ionC poistion is not 0')
    call kill(this)

    write(*,*)'test R attribute is stored properly in backup file'
    call make(this)
    this%R=1.0!First, manually set metiu attributes to non-default values
    call reset(this,state=1)
    call system('rm -f metiu.tmpfile')
    call backup(this,file='metiu.tmpfile')
    call kill(this)
    call make(this,file='metiu.tmpfile')
    !Then, assert non default attribute values are conserved
    call assert(this%R,1.0_double,msg='metiu R is not stored properly')
    call kill(this)    
    call system('rm -f metiu.tmpfile')


    write(*,*)'test make sets correct default value for R_c '
    call make(this)
    call assert(this%R_c,1.5*angstrom&
         ,msg='metiu default cut-off distance is not 1.5 angstrom')
    call kill(this)

    write(*,*)'test R_c attribute is stored properly in backup file'
    call make(this)
    this%R_c=1.0_double!First, manually set metiu attributes to non-default values
    call reset(this,state=1)
    call system('rm -f metiu.tmpfile')
    call backup(this,file='metiu.tmpfile')
    call kill(this)
    call make(this,file='metiu.tmpfile')
    !Then, assert non default attribute values are conserved
    call assert(this%R_c,1.0_double,msg='metiu R_c is not stored properly')
    call kill(this)    
    call system('rm -f metiu.tmpfile')

    write(*,*)'test make sets correct default value for nstate '
    call make(this)
    call assert(this%nstate.EQ.2,msg='metiu default nstate is not 2')
    call kill(this)

    write(*,*)'test nstate attribute is stored properly in backup file'
    call make(this)
    this%nstate=4!First, manually set metiu attributes to non-default values
    call reset(this,state=1)
    call system('rm -f metiu.tmpfile')
    call backup(this,file='metiu.tmpfile')
    call kill(this)
    call make(this,file='metiu.tmpfile')
    !Then, assert non default attribute values are conserved
    call assert(this%nstate.EQ.4,msg='metiu nstate is not stored properly')
    call kill(this)    
    call system('rm -f metiu.tmpfile')

    write(*,*)'test make sets correct default groundstate energy.'
    call make(this)
    call assert(this%Eigenvalue(0),0.0_double&
         ,msg='metiu default ground state energy is not 0.')
    call kill(this)




    !Reproduce V_b and Delta in figure 2 in JCP 102, 9285 (1995)
    write(*,*)'test V_b=0.87eV for R_c=1.50 Angstrom'
    call make(this)
    this%R_c=1.50*angstrom
    !this%R_c=2.50*angstrom
    call update(this)
    write(*,*)this%V_b,0.87*eV,abs(this%V_b-0.87*eV)/(0.87*eV)
    call assert(abs(this%V_b-0.87*eV)/(0.87*eV).LT.0.05,msg='metiu figure 2a mismatch V_b does&
         not equal 0.87eV within 5% error, for R_c=1.5 Angstrom')
    call kill(this)




    

    !================== consider the following tests ========================
    !-----  make tests -----
    !***          example          ****
    !write(*,*)'test make sets correct default values'
    !call make(this)
    !call assert(this%XXX.EQ.YYY,msg='metiu default XXX is not YYY')
    !call kill(this)
    !**********************************

    !----- kill tests -----
    !***          example          ****
    !write(*,*)'test kill cleans up dynamic memory and pointers'
    !call make(this)
    !call kill(this)
    !call assert(.not.associated(this%PPP),msg='metiu pointer PPP remains associated after killed.')
    !**********************************

    !----- backup tests -----
    !***          example          ****
    !write(*,*)'test attributes are stored properly stored in backup file'
    !call make(this)
    !this%var=XXX!First, manually set metiu attributes to non-default values
    !call system('rm -f metiu.tmpfile')
    !call backup(this,file='metiu.tmpfile')
    !call kill(this)
    !call make(this,file='metiu.tmpfile')
    !call assert(this%var.EQ.XXX)!Then, assert non default attribute values are conserved
    !call kill(this)    
    !**********************************

    !----- status tests -----

    !----- update tests -----

    !----- reset tests -----
    
    !----- fail cases -----

    !========================================================================



    write(*,*)'ALL metiu TESTS PASSED!'
  end subroutine metiu_test
  !-----------------------------------------

end module metiu_class

