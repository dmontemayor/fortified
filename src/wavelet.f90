!>\brief
!! wavelet class
!!\details
!! A wavefunction in vector representation of npt points with ndof degrees
!! of freedom.
!! The wavefunction is a list of npt points by 2+ndof columns. The first
!! two columns represent the real and imaginary weights of the vector
!! denoted by the following ndof columns.
!<------------------------------------------------------------------------
module wavelet_class
 use type_kinds
  implicit none
  private

  public::wavelet, wavelet_test
  public::describe, check, make, kill, backup, update, reset, status

  type wavelet
     logical::initialized=.false.
     character(len=label)::name='wavelet'
     
     !***************      Enter wavelet attributes here     ***************!
     integer(long)::ndof  !degrees of freedom
     integer(long)::npt  !number of points

     !wavefunction in delta function rep
     real(double),dimension(:,:),pointer::wf 
     !real(double),dimension(:,:),pointer::grid !delta function coordinates

     !***********************************************************************!
     !=======================================================================!
     !***   Here are some example attributes your derived type may have   ***!
     ! type(primitive)::primitive                      !an other type        !
     ! integer(short)::ierr                            !a short integer      !
     ! integer(long)::ndof                             !a long integer       !
     ! real(double)::var                               !a real variable      !
     ! complex(double)::zed                            !a complex variable   !
     ! real(double),dimension(:,:),pointer::matrix     !a real matrix        !
     ! complex(double),dimension(:,:),pointer::Zmatrix !a complex matrix     !
     !***********************************************************************!
  end type wavelet

  !> Creates the wavelet object.
  interface make
     module procedure wavelet_init
  end interface

  !> Destroys the wavelet object.
  interface kill
     module procedure wavelet_kill
  end interface

  !> Returns current state of the wavelet object.
  interface status
     module procedure wavelet_status
  end interface

  !> Returns a plain text description of the wavelet object.
  interface describe
     module procedure wavelet_describe
  end interface
  
  !> Returns the current state of the wavelet object.
  interface backup
     module procedure wavelet_backup
  end interface

  !> Recaluclates the wavelet object.
  interface update
     module procedure wavelet_update
  end interface

  !> Reinitializes the wavelet object.
  interface reset
     module procedure wavelet_reset
  end interface

  !> Checks that the wavelet object.
  interface check
     module procedure wavelet_check
  end interface

contains
  !======================================================================
  !> \brief Retruns a description of wavelet as a string.
  !> \param[in] THIS is the wavelet object.
  !======================================================================
  character(len=comment) function wavelet_describe(this)
    type(wavelet),intent(in)::this
    character(len=5)::FMT='(A)'

    write(wavelet_describe,FMT)'The wavelet is the basis of any quantum&
         & mechanical object. It exists in NDOF spatial dimensions&
         & distributed over NPT points in space. The wavefunction itself&
         & WF is complex valued.'
  end function wavelet_describe

  !======================================================================
  !> \brief Creates and initializes wavelet.
  !! \param THIS is the wavelet object.
  !! \param[in] FILE is an optional string containing the name of a
  !! previously backuped wavelet file.
  !=====================================================================
  subroutine wavelet_init(this,file)
    use filemanager
    use testing_class
    type(wavelet),intent(inout)::this
    character*(*),intent(in),optional::file
    integer(long)::unit
    logical::fileisopen=.false.
    character(len=label)::header
    integer(long)::i,j
    
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
       
       !check if file is of type wavelet
       read(unit,*)header
       call assert(trim(header).EQ.this%name,msg='wavelet_init: bad input file header in file'//file)
       
       !read static parameters
       !***             example            ***
       !***read scalar parameters from file***
       !read(unit,*)this%XXX
       !**************************************
       read(unit,*)this%ndof
       read(unit,*)this%npt
       
       !use reset to manage dynamic memory, reset sub-objects, and set random parameters
       call reset(this,state=1)
       
       !READ dynamic array values
       !***      example     ***
       !read(unit,*)(this%PPP(i),i=0,this%XXX-1)
       !************************
       read(unit,*)((this%wf(i,j),j=0,this%ndof+1),i=0,this%npt-1)
       !read(unit,*)((this%grid(i,j),j=0,this%ndof-1),i=0,this%npt-1)
       
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
       this%ndof=1
       this%npt=1
       
       !Use reset to make a default object
       call reset(this,state=1)
    end if

  end subroutine wavelet_init

  !======================================================================
  !> \brief Destroys the wavelet object.
  !> \param THIS is the wavelet object to be destroyed.
  !====================================================================
  subroutine wavelet_kill(this)
    type(wavelet),intent(inout)::this
 
    call reset(this,state=0)

  end subroutine wavelet_kill

  !======================================================================
  !> \brief Computes the current state of wavelet object.
  !> \param THIS is the wavelet  object to be updated.
  !======================================================================
  subroutine wavelet_update(this)
    type(wavelet),intent(inout)::this

    !Recompute any attribute values that might have evolved

    !****************************      Example     ******************************
    !*** attribute 'var' is always equall to the trace of the denisity matrix *** 
    ! this%var=0._double
    ! do istate=1,this%primitive%nstate
    !    this%var=this%var+this%den(istate,istate)
    ! end do
    !****************************************************************************

  end subroutine wavelet_update


  !======================================================================
  !> \brief Re-initiallizes the wavelet object.
  !> \param THIS is the wavelet  object to be re-initialized.
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
  subroutine wavelet_reset(this,state)
    type(wavelet),intent(inout)::this
    integer(long),intent(in),optional::STATE
    
    if(present(state))then
       if(state.EQ.0)then
          !nullify all pointers
          !******        Example - cleanup pointer attribute 'PPP'       ****
          !if(associated(this%PPP))nullify(this%PPP)
          !******************************************************************
          if(associated(this%wf))nullify(this%wf)
          !if(associated(this%grid))nullify(this%grid)
          
          !kill all objects
          !**** example **********
          !call kill(this%object)
          !***********************
          
          !set static parameters to error values
          this%ndof=-1
          this%npt=-1

          !un-initialized metiu object
          this%initialized=.false.
       else

          !allocate dynamic memory
          !***  Example - cleanup pointer attribute 'PPP'     ***
          !***            then reallocate memory              ***
          !if(associated(this%PPP))nullify(this%PPP)
          !allocate(this%PPP(0:this%npt-1))
          !******************************************************
          if(associated(this%wf))nullify(this%wf)
          allocate(this%wf(0:this%npt-1,0:this%ndof+1))

          !if(associated(this%grid))nullify(this%grid)
          !allocate(this%grid(0:this%npt-1,0:this%ndof-1))
          
          !Set default dynamic memory values
          !***  Example - set values in pointer 'PPP' to zero ***
          !this%PPP(:)=0.0_double
          !******************************************************
          !this%wf=(1.0_double,0.0_double)
          !this%grid=0.0_double
                    
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
       
       !Reset sub-objects
       !**** example **********
       !call reset(this%object)
       !***********************
       
    end if
    

  end subroutine wavelet_reset

  !======================================================================
  !> \brief Backups the current state of the wavelet object to file.
  !> \param[in] THIS is the wavelet  object to be updated.
  !> \param[in] FILE is a string containing the location of the backup file.
  !======================================================================
  subroutine wavelet_backup(this,file)
    use filemanager
    use testing_class
    use string
    type(wavelet),intent(in)::this
    character*(*),intent(in)::file
    integer(short)::unit
    logical::fileisopen
    integer(long)::i,j
    
    !check input file
    inquire(file=file,opened=fileisopen,number=unit)
    if(unit.LT.0)unit=newunit()
    if(.not.fileisopen)open(unit,file=file)
    
    !check wavelet object
    call assert(check(this).EQ.0,msg='wavelet object does not pass check.')

    !always write the data type on the first line
    write(unit,*)'wavelet'

    !******      Backup below all the derived type's attributes       ****
    !******         in the order the MAKE method reads them           ****

    !scalars     
    !******          Example - Backup a scalar attribute            ******
    ! write(unit,*)this%var
    !*********************************************************************
    write(unit,*)this%ndof
    write(unit,*)this%npt


    !dynamics arrays
    !***       Example - Backup an NxM matrix attribute                ***
    ! write(unit,*)((this%matrix(i,j),j=0,size(this%matrix,2)-1),i=0,size(this%matrix,1)-1)
    !*********************************************************************
    write(unit,*)((this%wf(i,j),j=0,size(this%wf,2)-1),i=0,size(this%wf,1)-1)
    !write(unit,*)((this%grid(i,j),j=0,size(this%grid,2)-1),i=0,size(this%grid,1)-1)


    !objects
    !******              Example - Backup an object            ***********
    ! call backup(this%primitive,file//'.primitive')
    ! write(unit,*)quote(file//'.primitive')!write the object location
    !*********************************************************************


    
    !finished writing all attributes - now close backup file
    close(unit)  
  end subroutine wavelet_backup
  
  !======================================================================
  !> \brief Retrun the current state of wavelet as a string.
  !> \param[in] THIS is the wavelet object.
  !> \param[in] MSG is an optional string message to annotate the status.
  !======================================================================
  character(len=line) function wavelet_status(this,msg)
    type(wavelet),intent(in)::this
    character*(*),intent(in),optional::msg
    character(len=5)::FMT='(A)'
    
    !Edit the status prompt to suit your needs
    write(wavelet_status,FMT)'wavelet status is currently not available'
    
  end function wavelet_status

 !======================================================================
  !> \brief Checks the wavelet object.
  !> \param[in] THIS is the wavelet object to be checked.
  !> \return 0 if all checks pass or exit at first failed check and returm non zero.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function wavelet_check(this)
    use testing_class
    type(wavelet),intent(in)::this

    !initiate with no problems found 
    wavelet_check=0

    !check that object is initialized
    call assert(this%initialized&
         ,msg='wavelet_check: wavelet object not initialized.'&
         ,iostat=wavelet_check)
    if(wavelet_check.NE.0)return

    !check that object has correct name
    call assert(this%name.EQ.'wavelet'&
         ,msg='wavelet_check: wavelet name is not set.'&
         ,iostat=wavelet_check)
    if(wavelet_check.NE.0)return
    
    !Check all attributes are within acceptable values

    !----------  ndof
    call assert(check(this%ndof).EQ.0&
         ,msg='wavelet_check: ndof failed check',iostat=wavelet_check)
    if(wavelet_check.NE.0)return

    call assert(this%ndof.LE.3&
         ,msg='wavelet_check: ndof>3 not allowed',iostat=wavelet_check)
    if(wavelet_check.NE.0)return

    call assert(this%ndof.GT.0&
         ,msg='wavelet_check: ndof is not positive.',iostat=wavelet_check)
    if(wavelet_check.NE.0)return


    !----------  npt
    call assert(check(this%npt).EQ.0&
         ,msg='wavelet_check: npt failed check',iostat=wavelet_check)
    if(wavelet_check.NE.0)return

    call assert(this%npt.GT.0&
         ,msg='wavelet_check: npt is not positive.',iostat=wavelet_check)
    if(wavelet_check.NE.0)return


    !----------- wf
    ! Description:
    ! wf is a quantum mechanical wavefunction. 
    ! wf should be integrable, complex, and in the delta representation.
    ! wf is a set of NPT complex weighted delta functions.
    ! delta function corrdianates are specified with the GRID attribute. 
    !------------

!!$    ! complex elements of wf should all be well behaved.
!!$    call assert(check(this%wf).EQ.0&
!!$         ,msg='wavelet_check: wf failed check',iostat=wavelet_check)
!!$    if(wavelet_check.NE.0)return
!!$
!!$    ! wf should form positive definate density
!!$    call assert(real(sum(this%wf*conjg(this%wf)),double).GT.epsilon(1.0_double)&
!!$         ,msg='wavelet_check: wf does not have positive density',iostat=wavelet_check)
!!$    if(wavelet_check.NE.0)return
!!$    
!!$    ! wf should be represented by NPT delta functions
!!$    call assert(size(this%wf).EQ.this%npt&
!!$         ,msg='wavelet_check: wf is not the same size as npt',iostat=wavelet_check)
!!$    if(wavelet_check.NE.0)return


!!$    !----------- grid
!!$    ! Description:
!!$    ! grid specifies the corrdinates of the NPT delta functions representing WF 
!!$    ! grid should be real and NPT x NDOF.
!!$    !------------
!!$
!!$    ! real elements of grid should all be well behaved.
!!$    call assert(check(this%grid).EQ.0&
!!$         ,msg='wavelet_check: grid failed check',iostat=wavelet_check)
!!$    if(wavelet_check.NE.0)return
!!$
!!$    ! grid should specify NPT x NDOF delta function coordinates
!!$    call assert(size(this%grid,1).EQ.this%npt&
!!$         ,msg='wavelet_check: grid dimension 1 is not the same size as npt',iostat=wavelet_check)
!!$    if(wavelet_check.NE.0)return
!!$    call assert(size(this%grid,2).EQ.this%ndof&
!!$         ,msg='wavelet_check: grid dimension 2 is not the same size as ndof',iostat=wavelet_check)
!!$    if(wavelet_check.NE.0)return



    !**********   Example - check an object attribute 'primitive'  *********
    !call assert(check(this%primitive).EQ.0&
    !     ,msg='wavelet_check: primitive sub-object failed check'&
    !     ,iostat=wavelet_check)
    !if(wavelet_check.NE.0)return
    !***********************************************************************
    
    !***   Example - check an integer attribute 'ndof' is well behaved   ***
    !call assert(check(this%ndof).EQ.0&
    !     ,msg='wavelet_check: ndof failed check',iostat=wavelet_check)
    !if(wavelet_check.NE.0)return
    !***********************************************************************
 
    !*** Example - add a constrain that says 'ndof' can only be positive ***
    !call assert(this%ndof.GT.0&
    !     ,msg='wavelet_check: ndof is not positive',iostat=wavelet_check)
    !if(wavelet_check.NE.0)return
    !***********************************************************************

    !***  Example - check a real valued attribute 'var' is well behaved  ***
    !call assert(check(this%var).EQ.0&
    !     ,msg='wavelet_check: var failed check',iostat=wavelet_check)
    !if(wavelet_check.NE.0)return
    !***********************************************************************

    !***  Example - add a constrain that says 'var' can not be zero     ***
    !call assert(abs(this%var).GT.epsilon(this%var)&
    !     ,msg='wavelet_check: var is tiny',iostat=wavelet_check)
    !if(wavelet_check.NE.0)return
    !***********************************************************************

    !***  Example - check a real valued pointer attribute 'matrix'       ***
    !***            is well behaved                                      ***
    !call assert(check(this%matrix).EQ.0&
    !     ,msg='wavelet_check: matrix failed check',iostat=wavelet_check)
    !if(wavelet_check.NE.0)return
    !***********************************************************************

    !********* Example - check an NxM matrix has right dimensions **********
    !call assert(size(this%matrix).EQ.N*M&
    !     ,msg='wavelet_check: number of matrix elements not = N*M.'&
    !     ,iostat=wavelet_check)
    !if(wavelet_check.NE.0)return
    !***********************************************************************

  end function wavelet_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the wavelet methods.
  !> \param[in] this is the wavelet object whose methods will be excercised.
  !> \return Nothing if all tests pass or 1 and a stop for the first failed test.
  !> \remark Will stop after first failed check.
  !======================================================================
  subroutine wavelet_test
    use testing_class
    use filemanager
    type(wavelet)::this
    character(len=label)::string
    integer(long)::unit
    
    !verify wavelet is compatible with current version
    include 'verification'

    write(*,*)'test NDOF attribute is stored properly in backup file'
    call make(this)
    this%ndof=2!First, manually set wavelet attributes to non-default values
    call reset(this,state=1)
    call system('rm -f wavelet.tmpfile*')
    call backup(this,file='wavelet.tmpfile')
    call kill(this)
    call make(this,file='wavelet.tmpfile')
    !Then, assert non default attribute values are conserved
    call assert(this%ndof.EQ.2,msg='wavelet NDOF not stored properly')
    call kill(this)    
    call system('rm -f wavelet.tmpfile*')

    write(*,*)'test defaut degrees of freedom is 1'
    call make(this)
    call assert(this%ndof.EQ.1,msg='wavelet default NDOF is not 1.')
    call kill(this)

    write(*,*)'test wavelet breaks with &
         & number of degrees of freedom greater than 3.'
    call make(this)
    this%ndof=4
    call reset(this,state=1)
    call assert(check(this).NE.0,msg='wavelet does not break with&
         & NDOF>3.')
    call kill(this)

    write(*,*)'test wavelet breaks with non-positive &
         & number of degrees of freedom.'
    call make(this)
    this%ndof=0
    call reset(this,state=1)
    call assert(check(this).NE.0,msg='wavelet does not break with&
         & non-positive NDOF.')
    call kill(this)

    write(*,*)'test NPT attribute is stored properly in backup file'
    call make(this)
    this%npt=2!First, manually set wavelet attributes to non-default values
    call reset(this,state=1)
    call system('rm -f wavelet.tmpfile*')
    call backup(this,file='wavelet.tmpfile')
    call kill(this)
    call make(this,file='wavelet.tmpfile')
    !Then, assert non default attribute values are conserved
    call assert(this%npt.EQ.2,msg='wavelet NPT not stored properly')
    call kill(this)    
    call system('rm -f wavelet.tmpfile*')

    write(*,*)'test defaut number of points is 1'
    call make(this)
    call assert(this%npt.EQ.1,msg='wavelet default NPT is not 1.')
    call kill(this)

    write(*,*)'test wavelet breaks with huge number of points.'
    call make(this)
    this%npt=huge(this%npt)
    call reset(this,state=1)
    call assert(check(this).NE.0,msg='wavelet does not break with&
         & huge NPT.')
    call kill(this)

    write(*,*)'test wavelet breaks with non-positive &
         & number of points.'
    call make(this)
    this%npt=0
    call reset(this,state=1)
    call assert(check(this).NE.0,msg='wavelet does not break with&
         & non-positive NPT.')
    call kill(this)

    write(*,*)'test make allocates wf pointer'
    call make(this)
    call assert(associated(this%wf),msg='wavelet pointer WF is not&
         & associated after make.')
    call kill(this)

    write(*,*)'test kill cleans up wf pointer'
    call make(this)
    call kill(this)
    call assert(.not.associated(this%wf),msg='wavelet pointer WF remains associated after killed.')

    write(*,*)'test WF attribute is stored properly in backup file'
    call make(this)
    this%wf=0.5!First, manually set wavelet attributes to non-default values
    call update(this)
    call system('rm -f wavelet.tmpfile*')
    write(*,*)this%wf
    write(*,*)size(this%wf,1),size(this%wf,2)
    call backup(this,file='wavelet.tmpfile')
    call kill(this)
    call make(this,file='wavelet.tmpfile')
    !Then, assert non default attribute values are conserved
    write(*,*)this%wf
    write(*,*)size(this%wf,1),size(this%wf,2)
    call assert(all(this%wf.EQ.0.5),msg='wavelet WF not stored properly')
    call kill(this)    
    call system('rm -f wavelet.tmpfile*')



!!$
!!$    write(*,*)'test wavelet breaks with non-positive basis set size.'
!!$    call make(this)
!!$    this%npt=0
!!$    call reset(this,state=1)
!!$    call assert(check(this).NE.0,msg='wavelet does not break with non-positive basis set size.')
!!$    call kill(this)
!!$
!!$    write(*,*)'test default basis set size is 1.'
!!$    call make(this)
!!$    call assert(this%npt.EQ.1,msg='wavelet default basis set size is not 1.')
!!$    call kill(this)

!!$    write(*,*)'test wavelet breaks if density is not positive definite'
!!$    call make(this)
!!$    this%wf=0._double
!!$    call assert(check(this).NE.0,msg='wavelet does not break with non-positive definite density.')
!!$    call kill(this)
!!$

!!$    write(*,*)'test kill cleans up grid pointer'
!!$    call make(this)
!!$    call kill(this)
!!$    call assert(.not.associated(this%grid),msg='wavelet pointer GRID remains associated after killed.')


!!$    write(*,*)'test wavelet breaks if WF size is not basis set size NPT.'
!!$    call make(this)
!!$    this%npt=size(this%wf)+1
!!$    call assert(check(this).NE.0,msg='wavelet does not break when NPT is not the same size as WF')
!!$    call kill(this)

!!$    write(*,*)'test wavelet breaks if GRID size is not NPT x NDOF.'
!!$    call make(this)
!!$    this%npt=size(this%grid,1)+1
!!$    this%ndof=size(this%grid,2)+1
!!$    call assert(check(this).NE.0,msg='wavelet does not break when GRID is not the same size as NPT x NDOF')
!!$    call kill(this)



!!$    write(*,*)'test GRID attribute is stored properly in backup file'
!!$    call make(this)
!!$    this%grid=1.5!First, manually set wavelet attributes to non-default values
!!$    call system('rm -f wavelet.tmpfile')
!!$    call backup(this,file='wavelet.tmpfile')
!!$    call kill(this)
!!$    call make(this,file='wavelet.tmpfile')
!!$    call system('rm -f wavelet.tmpfile')
!!$    !Then, assert non default attribute values are conserved
!!$    call assert(all(this%grid.EQ.1.5),msg='wavelet GRID not stored properly')
!!$    call kill(this)    

    !================== consider the following tests ========================
    !-----  make tests -----
    !***          example          ****
    !write(*,*)'test make sets correct default values'
    !call make(this)
    !call assert(this%XXX.EQ.YYY,msg='wavelet default XXX is not YYY')
    !call kill(this)
    !**********************************

    !----- kill tests -----
    !***          example          ****
    !write(*,*)'test kill cleans up dynamic memory and pointers'
    !call make(this)
    !call kill(this)
    !call assert(.not.associated(this%PPP),msg='wavelet pointer PPP remains associated after killed.')
    !**********************************

    !----- backup tests -----
    !***          example          ****
    !write(*,*)'test attributes are stored properly stored in backup file'
    !call make(this)
    !this%var=XXX!First, manually set wavelet attributes to non-default values
    !call system('rm -f wavelet.tmpfile')
    !call backup(this,file='wavelet.tmpfile')
    !call kill(this)
    !call make(this,file='wavelet.tmpfile')
    !call system('rm -f wavelet.tmpfile')
    !!Then, assert non default attribute values are conserved
    !call assert(this%var.EQ.XXX)
    !call kill(this)    
    !**********************************

    !----- status tests -----

    !----- update tests -----

    !----- reset tests -----
    
    !----- fail cases -----

    !========================================================================
    write(*,*)'ALL wavelet TESTS PASSED!'
  end subroutine wavelet_test
  !-----------------------------------------

end module wavelet_class

