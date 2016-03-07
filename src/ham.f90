!>\brief
!! ham class
!!\details
!! Hamiltonian representing the greatest common denominator among propagators.
!! All propagators must operate on the system defined exclusively by this Hamiltonian.
!<------------------------------------------------------------------------
module ham_class
  use type_kinds
  use quasiparticle_class
  implicit none
  private

  public::ham, ham_test
  public::make, kill, display, backup, update, reset, check


  type ham
     logical::initialized=.false.

     !space
     integer(long)::ndim           !number of spatial dimensions
     real(long),pointer::Lcell(:)  !length of orthorhombic unit cell in each dimension
     logical,pointer::PBC(:)       !periodic boundary conditions in each dimension

     !matter and energy
     integer(long)::nstate                !number of quantum states
     integer(long)::npart                 !number of particles
     type(quasiparticle),pointer::Z(:,:)  !npart particles in ndim spatial dimensions

  end type ham

  !> Creates the ham object.
  interface make
     module procedure ham_init
  end interface

  !> Destroys the ham object.
  interface kill
     module procedure ham_kill
  end interface

  !> Displays the current state of the ham object.
  interface display
     module procedure ham_display
  end interface

  !> Backups the current state of the ham object.
  interface backup
     module procedure ham_backup
  end interface

  !> Recaluclates the ham object.
  interface update
     module procedure ham_update
  end interface

  !> Reinitializes the ham object.
  interface reset
     module procedure ham_reset
  end interface

  !> Checks that the ham object.
  interface check
     module procedure ham_check
  end interface

contains

  !======================================================================
  !> \brief Creates and initializes the ham object.
  !> \param this is the ham object to be initialized.
  !> \param[in] file is an optional string containing the name of a previously backupd ham file.
!!$  !> \remark If no input file is provided the user must manually initialize THIS using stout.
  !=====================================================================
  subroutine ham_init(this,file,ndim,nstate,npart)
    use filemanager
    use testing_class
    type(ham),intent(inout)::this
    character*(*),intent(in),optional::file
    integer(long),intent(in),optional::ndim,nstate,npart
    integer(long)::unit
    logical::fileisopen=.false.
    character(len=label)::header
    integer(long)::i,j

    !set scalar parameters
    if(present(file))then 
       
       !check file status    
       inquire(file=file,opened=fileisopen,number=unit)
       if(unit.LT.0)unit=newunit()
       if(.not.fileisopen)open(unit,file=file)
    
       !check if file is of type ham
       read(unit,*)header
       call assert(trim(header).EQ.'ham',msg='ham_init: bad input file header in file'//file)
       
       !read scalar parameters
       read(unit,*)this%ndim
       read(unit,*)this%nstate
       read(unit,*)this%npart
    else
       !set default scalar parameters
       this%ndim=1
       if(present(ndim))this%ndim=ndim
       this%nstate=2
       if(present(nstate))this%nstate=nstate
       this%npart=1
       if(present(npart))this%npart=npart
    end if

    !allocate dynamic arrays
    if(associated(this%Lcell))nullify(this%Lcell)
    allocate(this%Lcell(0:this%ndim-1))
    if(associated(this%PBC))nullify(this%PBC)
    allocate(this%PBC(0:this%ndim-1))
    if(associated(this%Z))nullify(this%Z)
    allocate(this%Z(0:this%npart-1,0:this%ndim-1))

    !set dynamic arrays
    if(present(file))then 
       !read dynamic arrays
       read(unit,*)(this%Lcell(i),i=0,this%ndim-1)
       read(unit,*)(this%PBC(i),i=0,this%ndim-1)
       do i=0,this%npart-1
          do j=0,this%ndim-1
             call make(this%Z(i,j),file=file)
          end do
       end do
    else
       !set default array values
       this%Lcell=0._double
       this%PBC=.false.
       do i=0,this%npart-1
          do j=0,this%ndim-1
             call make(this%Z(i,j))
          end do
       end do
    end if
    
    !finished reading all attributes - now close backup file
    if(present(file))close(unit)


    !declare initialization complete
    this%initialized=.true.

  end subroutine ham_init

  !======================================================================
  !> \brief Destroys the ham object.
  !> \param this is the ham object to be destroyed.
  !====================================================================
  subroutine ham_kill(this)
    type(ham),intent(inout)::this
    integer(long)::i,j 

    !nullify all pointers to dynamic arrays
    if(associated(this%Lcell))nullify(this%Lcell)
    if(associated(this%PBC))nullify(this%PBC)
    do i=0,this%npart-1
       do j=0,this%ndim-1
          call kill(this%Z(i,j))
       end do
    end do
    if(associated(this%Z))nullify(this%Z)

    !un-initialized ham object
    this%initialized=.false.

  end subroutine ham_kill

  !======================================================================
  !> \brief Computes the current state of ham object.
  !> \param this is the ham  object to be updated.
  !======================================================================
  subroutine ham_update(this)
    type(ham),intent(inout)::this

!!$    call Note('Begin ham_update.')
!!$    
!!$    !Primitives usually dont get updated
!!$    !    call update(this%primitive)    
!!$
!!$
!!$    !******   Recompute any attribute values that might have evolved   ******!
!!$
!!$
!!$
!!$
!!$
!!$    !************************************************************************!
!!$    !========================================================================!
!!$    !*****    Example - attribute 'var' is always equall to the trace   *****!
!!$    !                   of the primitive's denisity                          !
!!$    !                                                                        !
!!$    ! this%var=0._double                                                     !
!!$    ! do istate=1,this%primitive%nstate                                             !
!!$    !    this%var=this%var+this%primitive%den(istate,istate)                        !
!!$    ! end do                                                                 !
!!$    !                                                                        !
!!$    !************************************************************************!
!!$
!!$    !Usually leave out final check before we exit
!!$    !if(check(this).EQ.1)call stop('ham_update: failed final check!')

  end subroutine ham_update

  !======================================================================
  !> \brief Re-initiallizes the ham object.
  !> \param this is the ham  object to be re-initialized.
  !======================================================================
  subroutine ham_reset(this)
    type(ham),intent(inout)::this
!!$    integer(long)::istate,jstate
!!$
!!$    call Note('Begin ham_reset.')
!!$
!!$    !reset the primitive
!!$    call reset(this%primitive)
!!$
!!$    !****  Reset any attributes to satisfy re-initiation conditions   ****! 
!!$
!!$
!!$
!!$
!!$
!!$    !************************************************************************!
!!$    !========================================================================!
!!$    !********    Example - attribute 'var' is always initially a     ********!
!!$    !                      Gaussian random number                            !
!!$    !                                                                        !
!!$    ! this%var=gran()                                                        !
!!$    !                                                                        !
!!$    !************************************************************************!
!!$
!!$    !update now that we have changed some values and do a final check
!!$    call update(this)
!!$    if(check(this).EQ.1)call stop('ham_reset: failed final check!')

  end subroutine ham_reset

  !======================================================================
  !> \brief Backups the current state of the ham object to file.
  !> \param[in] this is the ham  object to be updated.
  !> \param[in] file is a string containing the location of the backup file.
  !======================================================================
  subroutine ham_backup(this,file)
    use filemanager
    type(ham),intent(in)::this
    character*(*),intent(in)::file
    integer(short)::unit
    logical::fileisopen
    integer(long)::i,j
    
    !check file status
    inquire(file=file,opened=fileisopen,number=unit)
    if(unit.LT.0)unit=newunit()
    if(.not.fileisopen)open(unit,file=file)
    
    !always write the data type on the first line
    write(unit,*)'ham'

    !backup scalar parameters
    write(unit,*)this%ndim
    write(unit,*)this%nstate
    write(unit,*)this%npart

    !backup dynamic arrays
    write(unit,*)(this%Lcell(i),i=0,this%ndim-1)
    write(unit,*)(this%PBC(i),i=0,this%ndim-1)

    !quasiparticle list
    do i=0,this%npart-1
       do j=0,this%ndim-1
          call backup(this%Z(i,j),file=file)
       end do
    end do
    !finished saving all attributes - now close backup file
    close(unit)

  end subroutine ham_backup

  !======================================================================
  !> \brief Retrun the ham object as a single line record entry.
  !> \param[in] this is the ham object.
  !> \param[in] msg is an optional string message to annotate the displayed object.
  !======================================================================
  character(len=line) function ham_display(this,msg)
    type(ham),intent(in)::this
    character*(*),intent(in),optional::msg
    character(len=5)::FMT='(A)'

    write(ham_display,FMT)'ham'

   
  end function ham_display

  !======================================================================
  !> \brief Checks the ham object.
  !> \param[in] this is the ham object to be checked.
  !> \return Nothing if all checks pass or 1 and a warn for the first failed check.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function ham_check(this)
    use testing_class
    type(ham),intent(in)::this
    integer:: i,j

    !initiate with no problems found 
    ham_check=0
    !call Note('Checking ham.')

    !check that object is initialized
    call assert(this%initialized,msg='ham_check: ham object not initialized.',iostat=ham_check)
    if(ham_check.NE.0)return

    !check that ndim is not NAN
    call assert(this%ndim.EQ.this%ndim,msg='ham_check: ndim is NAN',iostat=ham_check)
    if(ham_check.NE.0)return

    !check that ndim is not huge
    call assert(this%ndim.LT.huge(1_long),msg='ham_check: ndim is huge',iostat=ham_check)
    if(ham_check.NE.0)return

    !check that ndim is positive
    call assert(this%ndim.GT.0,msg='ham_check: ndim is not greater than zero',iostat=ham_check)
    if(ham_check.NE.0)return

    !check Lcell points to something
    call assert(associated(this%Lcell),msg='ham_check: Lcell is not associated',iostat=ham_check)
    if(ham_check.NE.0)return

    !check Lcell is of dimension ndim
    call assert(size(this%Lcell).EQ.this%ndim,msg='ham_check: Lcell is not of dimension this%ndim',iostat=ham_check)
    if(ham_check.NE.0)return

    !check Lcell is not NAN
    call assert(all(this%Lcell.EQ.this%Lcell),msg='ham_check: Lcell has NAN values',iostat=ham_check)
    if(ham_check.NE.0)return

    !check Lcell is not huge
    call assert(all(this%Lcell.LT.huge(1._double)),msg='ham_check: Lcell has huge values',iostat=ham_check)
    if(ham_check.NE.0)return

    !check Lcell is not negative (zero vauled dimension means that cell length is infinite in that dimension
    call assert(all(this%Lcell.GE.0),msg='ham_check: Lcell has negative values',iostat=ham_check)
    if(ham_check.NE.0)return

    !check periodic boundary conditions are associated
    call assert(associated(this%PBC),msg='ham_check: PBC is not associated',iostat=ham_check)
    if(ham_check.NE.0)return

    !check PBC is of dimension ndim
    call assert(size(this%PBC).EQ.this%ndim,msg='ham_check: PBC is not of dimension this%ndim',iostat=ham_check)
    if(ham_check.NE.0)return

    !check PBC is false for any dimension with cell length zero
    call assert(all(.not.(this%Lcell.EQ.0.and.this%PBC)),&
         msg='ham_check: PBC is set for infinite length dimension',iostat=ham_check) 
    if(ham_check.NE.0)return

    !check number of quantum states is not NAN
    call assert(this%nstate.EQ.this%nstate,msg='ham_check: nstate is NAN',iostat=ham_check)
    if(ham_check.NE.0)return

    !check number of quantum states is not huge
    call assert(this%nstate.LT.huge(1_long),msg='ham_check: nstate is huge',iostat=ham_check)
    if(ham_check.NE.0)return

    !check number of quantum states is greater than 1
    call assert(this%nstate.GT.1,msg='ham_check: nstate is not greater than 1',iostat=ham_check)
    if(ham_check.NE.0)return

    !check number of quasi-particles is not NAN
    call assert(this%npart.EQ.this%npart,msg='ham_check: npart is NAN',iostat=ham_check)
    if(ham_check.NE.0)return

    !check number of quasi-particles is not huge
    call assert(this%npart.LT.huge(1_long),msg='ham_check: npart is huge',iostat=ham_check)
    if(ham_check.NE.0)return

    !check number of quasi-particles is positive
    call assert(this%npart.GT.0,msg='ham_check: npart is not positive integer',iostat=ham_check)
    if(ham_check.NE.0)return

    !check particle array points to something
    call assert(associated(this%Z),msg='ham_check: quasi-particle array Z is not associated',iostat=ham_check)
    if(ham_check.NE.0)return

    !check particle arrray has dimensions npart by ndim
    call assert(size(this%Z,1).EQ.this%npart,msg='ham_check: Z array dimension 1 is not this%npart',iostat=ham_check)
    if(ham_check.NE.0)return
    call assert(size(this%Z,2).EQ.this%ndim,msg='ham_check: Z array dimension 2 is not this%ndim',iostat=ham_check)
    if(ham_check.NE.0)return

    !check particle array elements
    do i=0,this%npart-1
       do j=0,this%ndim-1
          call assert(check(this%Z(i,j)).EQ.0&
               ,msg='ham_check: quasi-particle element in Z array has failed check',iostat=ham_check)
          if(ham_check.NE.0)return
       end do
    end do

  end function ham_check

  !======================================================================
  !> \brief Tests the ham methods.
  !> \param[in] this is the ham object whose methods will be excercised.
  !> \return Nothing if all tests pass or 1 and a stop for the first failed test.
  !> \remark Will stop after first failed check.
  !======================================================================
  subroutine ham_test
    use testing_class
    type(ham)::this
    character(len=label)::string
    character(len=line)::record
    integer(long)::i,j    


    write(*,*)'test ham can be checked prior to being made.'
    call assert(check(this).EQ.1,msg='check ham object does not return 1 prior to make.')

    write(*,*)'test ham can be created.'
    call make(this)
    call assert(check(this).EQ.0,msg='ham object was not created properly.')
    write(*,*)'test ham kill method sets initiallization flag to false.'
    call kill(this)
    call assert(.not.this%initialized,msg='ham object remains initialized after killed.')
    write(*,*)'test kill method cleans up dynamic memory and pointers'
    call assert(.not.associated(this%Lcell),msg='ham pointer Lcell remains associated after killed.')
    call assert(.not.associated(this%PBC),msg='ham pointer PBC remains associated after killed.')
    call assert(.not.associated(this%Z),msg='ham pointer Z remains associated after killed.')

    write(*,*)'check default values are set correctly'
    call make(this)
    call assert(this%ndim.EQ.1,msg='default ndim is not 1')
    call assert(this%Lcell(0).EQ.0,msg='default Lcell is not 0')
    call assert(.not.this%PBC(0),msg='default PBC is not false')
    call assert(this%nstate.EQ.2,msg='default nstate is not 2')
    call assert(this%npart.EQ.1,msg='default npart is not 1')
    call kill(this)

    write(*,*)'test ham can be displayed '
    call make(this)
    record=display(this)
    call assert(trim(record).EQ.'ham',msg='ham object does not display properly.')
    call kill(this)

    write(*,*)'test ham can be backupd'
    call make(this)
    call system('rm -f testham.tmpfile')
    call backup(this,file='testham.tmpfile')
    call assert('testham.tmpfile',msg='ham backup file was not created.')
    call system('rm -f testham.tmpfile')
    call kill(this)

    write(*,*)'test ham can be created with backup file'
    call make(this)
    call backup(this,file='testham.tmpfile')
    call kill(this)
    call make(this,file='testham.tmpfile')
    call assert(check(this).EQ.0,msg='ham object was not created properly from backupfile.')
    call system('rm -f testham.tmpfile')
    
    write(*,*)'test backup file begins with ham object name in first line'
    call make(this)
    call system('rm -f testham.tmpfile')
    call backup(this,file='testham.tmpfile')
    open(123,file='testham.tmpfile')
    read(123,*)string
    call assert(trim(string).EQ.'ham',msg='save file does not have ham label on first line.')
    call system('rm -f testham.tmpfile')
    call kill(this)
    
    write(*,*)'test ham attributes are properly backupd in backupfile'
    call make(this,ndim=3,nstate=3,npart=2)
    !set non default attribute values
    this%Lcell=1._double
    this%PBC=.True.
    do i=0,this%npart-1
       do j=0,this%ndim-1
         call make(this%Z(i,j),npt=2)
          this%Z(i,j)%mass=2.0_double
          this%Z(i,j)%charge=2.0_double
          this%Z(i,j)%den=.5_double
          this%Z(i,j)%pos=1.0_double
          this%Z(i,j)%mom=1.0_double
          this%Z(i,j)%force=1.0_double
       end do
    end do
    call system('rm -f testham.tmpfile')
    call backup(this,file='testham.tmpfile')
    call kill(this)
    call make(this,file='testham.tmpfile')
    !assert non defulat attribute values are conserved
    call assert(this%ndim.EQ.3,msg='ham ndim is not 3')
    call assert(all(this%Lcell.EQ.1.0_double),msg='ham Lcell is not 1')
    call assert(all(this%PBC),msg='ham PBC is not true')
    call assert(this%nstate.EQ.3,msg='ham nstate is not 3')
    call assert(this%npart.EQ.2,msg='ham npart is not 2')
    do i=0,this%npart-1
       do j=0,this%ndim-1
          call assert(this%Z(i,j)%mass,2.0_double,msg='ham quasiparticle mass is not 2')
          call assert(this%Z(i,j)%charge,2.0_double,msg='ham quasiparticle charge is not 2')
          call assert(this%Z(i,j)%den(0),0.5_double,msg='ham quasiparticle density is not 0.5')
          call assert(this%Z(i,j)%pos(0),1.0_double,msg='ham quasiparticle position is not 1')
          call assert(this%Z(i,j)%mom(0),1.0_double,msg='ham quasiparticle momentum is not 1')
          call assert(this%Z(i,j)%force(0),1.0_double,msg='ham quasiparticle force is not 1')
       end do
    end do
    call kill(this)    

    write(*,*)'test ham can be updated'
    call make(this)
    call update(this)
    call assert(check(this).EQ.0,msg='ham object was not updated properly.')
    call kill(this)
    
    write(*,*)'test ham can be resetted'
    call make(this)
    call reset(this)
    call assert(check(this).EQ.0,msg='ham object was not resetted properly.')
    call kill(this)
    

    write(*,*)'ALL ham TESTS PASSED!'
  end subroutine ham_test
  !-----------------------------------------

end module ham_class

