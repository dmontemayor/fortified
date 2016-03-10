!>\brief
!! quasiparticle class
!!\details
!! Use this quasiparticle to help you start building a new object and methods.
!! Do not remove the \a make, \a kill, \a status, \a backup, \a update, \a reset, \a check, and \a test methods.
!! These methods must be present for the class to integrate properly.
!! You can leave these methods blank if you do not want those particular funcitionalities, although, this is not advised.
!! Follow the examples commented in the source code to define the various attributes of your class.
!! DO NOT OVERWRITE THIS FILE. Make a copy and title it \a myclass.f90 where 'myclass' is your choice
!! one-word title for your new class.
!! Finally, Replace instances of the string 'quasiparticle' in this quasiparticle with the same
!! one-word title 'myclass' you chose.
!<------------------------------------------------------------------------
module quasiparticle_class
  use type_kinds
  implicit none
  private

  public::quasiparticle, quasiparticle_test
  public::make, kill, status, backup, update, reset, check

  type quasiparticle                  !quasi-particle with single degree of freedom
     logical::initialized=.false.
     character(len=label)::name='quasiparticle'
     real(double)::mass               !quasi-particle total mass
     real(double)::charge             !quasi-particle total charge
     integer(long)::npt               !density is on a grid of npt points to accomodate delocalization
     real(double),pointer::den(:)     !quasi-particle density for each point on grid 
     real(double),pointer::pos(:)     !quasi-particle position on grid
     real(double),pointer::mom(:)     !quasi-particle momentum on grid
     real(double),pointer::force(:)   !quasi-particle force on grid
  end type quasiparticle

  !> Creates the quasiparticle object.
  interface make
     module procedure quasiparticle_init
  end interface

  !> Destroys the quasiparticle object.
  interface kill
     module procedure quasiparticle_kill
  end interface

  !> Returns the current state of the quasiparticle object.
  interface status
     module procedure quasiparticle_status
  end interface

  !> Returns a description of the quasiparticle object.
  interface describe
     module procedure quasiparticle_describe
  end interface

  !> Backups the current state of the quasiparticle object.
  interface backup
     module procedure quasiparticle_backup
  end interface

  !> Recaluclates the quasiparticle object.
  interface update
     module procedure quasiparticle_update
  end interface

  !> Reinitializes the quasiparticle object.
  interface reset
     module procedure quasiparticle_reset
  end interface

  !> Checks that the quasiparticle object.
  interface check
     module procedure quasiparticle_check
  end interface

contains
  !======================================================================
  !> \brief Retruns a description of template as a string.
  !> \param[in] this is the template object.
  !======================================================================
  character(len=comment) function quasiparticle_describe(this)
    type(quasiparticle),intent(in)::this
    character(len=5)::FMT='(A)'

    write(quasiparticle_describe,FMT)'No description for quasiparticle has been provided.'
   
  end function quasiparticle_describe

  !======================================================================
  !> \brief Creates and initializes the quasiparticle object.
  !> \param this is the quasiparticle object to be initialized.
  !> \param[in] file is an optional string containing the name of a previously backupd quasiparticle file.
!!$  !> \remark If no input file is provided the user must manually initialize THIS using stout.
  !=====================================================================
  subroutine quasiparticle_init(this,file,npt)
    use filemanager
    use testing_class
    type(quasiparticle),intent(inout)::this
    character*(*),intent(in),optional::file
    integer(long),intent(in),optional::npt
    integer(long)::unit
    logical::fileisopen=.false.
    character(len=label)::header
    integer(long)::i

    !set scalar parameters
    if(present(file))then 
       
       !check file status    
       inquire(file=file,opened=fileisopen,number=unit)
       if(unit.LT.0)unit=newunit()
       if(.not.fileisopen)open(unit,file=file)
    
       !check if file is of type quasiparticle
       read(unit,*)header
       call assert(trim(header).EQ.'quasiparticle',msg='quasiparticle_init: bad input file header in file'//file)
       
       !read scalar parameters
       read(unit,*)this%mass
       read(unit,*)this%charge
       read(unit,*)this%npt
    else
       !set default scalar parameters
       this%mass=1.0_double
       this%charge=0.0_double
       this%npt=1
       if(present(npt))this%npt=npt
    end if

    !allocate dynamic arrays
    if(associated(this%den))nullify(this%den)
    allocate(this%den(0:this%npt-1))
    if(associated(this%pos))nullify(this%pos)
    allocate(this%pos(0:this%npt-1))
    if(associated(this%mom))nullify(this%mom)
    allocate(this%mom(0:this%npt-1))
    if(associated(this%force))nullify(this%force)
    allocate(this%force(0:this%npt-1))

    !set dynamic arrays
    if(present(file))then 
       !read dynamic arrays
       read(unit,*)(this%den(i),i=0,this%npt-1)
       read(unit,*)(this%pos(i),i=0,this%npt-1)
       read(unit,*)(this%mom(i),i=0,this%npt-1)
       read(unit,*)(this%force(i),i=0,this%npt-1)
    else
       !set default array values
       this%den=1._double/real(this%npt,double)
       this%pos=0._double
       this%mom=0._double
       this%force=0._double
    end if
    
    !finished reading all attributes - now close backup file
    if(.not.fileisopen)close(unit)

    !declare initialization complete
    this%initialized=.true.

  end subroutine quasiparticle_init

  !======================================================================
  !> \brief Destroys the quasiparticle object.
  !> \param this is the quasiparticle object to be destroyed.
  !====================================================================
  subroutine quasiparticle_kill(this)
    type(quasiparticle),intent(inout)::this
 
    !nullify pointers to dynamic arrays
    if(associated(this%den))nullify(this%den)
    if(associated(this%pos))nullify(this%pos)
    if(associated(this%mom))nullify(this%mom)
    if(associated(this%force))nullify(this%force)

    !un-initialized quasiparticle object
    this%initialized=.false.

  end subroutine quasiparticle_kill

  !======================================================================
  !> \brief Computes the current state of quasiparticle object.
  !> \param this is the quasiparticle  object to be updated.
  !======================================================================
  subroutine quasiparticle_update(this)
    type(quasiparticle),intent(inout)::this

!!$    call Note('Begin quasiparticle_update.')
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
!!$    !if(check(this).EQ.1)call stop('quasiparticle_update: failed final check!')

  end subroutine quasiparticle_update

  !======================================================================
  !> \brief Re-initiallizes the quasiparticle object.
  !> \param this is the quasiparticle  object to be re-initialized.
  !======================================================================
  subroutine quasiparticle_reset(this)
    type(quasiparticle),intent(inout)::this
!!$    integer(long)::istate,jstate
!!$
!!$    call Note('Begin quasiparticle_reset.')
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
!!$    if(check(this).EQ.1)call stop('quasiparticle_reset: failed final check!')

  end subroutine quasiparticle_reset

  !======================================================================
  !> \brief Backups the current state of the quasiparticle object to file.
  !> \param[in] this is the quasiparticle  object to be updated.
  !> \param[in] file is a string containing the location of the backup file.
  !> \param[in] append is a boolean input determining if file is appended to or overwritten
  !======================================================================
  subroutine quasiparticle_backup(this,file)
    use filemanager

    type(quasiparticle),intent(in)::this
    character*(*),intent(in)::file
    integer(short)::unit
    logical::fileisopen
    integer(long)::i

    !check file    
    inquire(file=file,opened=fileisopen,number=unit)
    if(unit.LT.0)unit=newunit()
    if(.not.fileisopen)open(unit,file=file)
    
    !always write the data type on the first line
    write(unit,*)'quasiparticle'

    !backup scalar parameters
    write(unit,*)this%mass
    write(unit,*)this%charge
    write(unit,*)this%npt

    !backup dynamic arrays
    write(unit,*)(this%den(i),i=0,this%npt-1)
    write(unit,*)(this%pos(i),i=0,this%npt-1)
    write(unit,*)(this%mom(i),i=0,this%npt-1)
    write(unit,*)(this%force(i),i=0,this%npt-1)

    !finished saving all attributes - now close backup file if previously closed
    if(.not.fileisopen)close(unit)

  end subroutine quasiparticle_backup

  !======================================================================
  !> \brief Retrun the quasiparticle object as a string.
  !> \param[in] this is the quasiparticle object.
  !> \param[in] msg is an optional string message to annotate the status.
  !> \remark will only return mass charge and zeroth point position and momentum
  !======================================================================
  character(len=line) function quasiparticle_status(this,msg)
    type(quasiparticle),intent(in)::this
    character*(*),intent(in),optional::msg
    character(len=10)::FMT='(4(E15.8))'

    write(quasiparticle_status,FMT)this%mass,this%charge,this%pos(0),this%mom(0)
   
  end function quasiparticle_status

  !======================================================================
  !> \brief Checks the quasiparticle object.
  !> \param[in] this is the quasiparticle object to be checked.
  !> \return 0 if all checks pass or exit at first failed check and returm non zero.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function quasiparticle_check(this)
    use testing_class
    type(quasiparticle),intent(in)::this

    !initiate with no problems found 
    quasiparticle_check=0

    !check that object is initialized
    call assert(this%initialized,msg='quasiparticle_check: quasiparticle object not initialized.',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check mass is not NAN
    call assert(this%mass.EQ.this%mass,msg='quasiparticle_check: mass is NAN',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check that mass is positive
    call assert(this%mass.GT.0,msg='quasiparticle_check: mass is not positive',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check that mass is not huge
    call assert(this%mass.LT.huge(1._double),msg='quasiparticle_check: mass is huge',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check charge is not NAN
    call assert(this%charge.EQ.this%charge,msg='quasiparticle_check: charge is NAN',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check that abs(charge) is not huge
    call assert(abs(this%charge).LT.huge(1._double),msg='quasiparticle_check: abs(charge) is huge',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check npt is not NAN
    call assert(this%npt.EQ.this%npt,msg='quasiparticle_check: npt is NAN',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check that npt is positive
    call assert(this%npt.GT.0,msg='quasiparticle_check: npt is not positive',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check that npt is not huge
    call assert(this%npt.LT.huge(1_long),msg='quasiparticle_check: npt is huge',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check that density points to something
    call assert(associated(this%den),msg='quasiparticle_check: density is not associated',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check density is dimension npt
    call assert(size(this%den).EQ.this%npt,msg='quasiparticle_check: density is not of size npt',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check that density does not contain NANs
    call assert(all(this%den.EQ.this%den),msg='quasiparticle_check: density contains NANs',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check that density does not contain huge values
    call assert(all(abs(this%den).LT.huge(1.0_double))&
         ,msg='quasiparticle_check: density has huge values',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check that density is conserved
    call assert(sum(this%den),1.0_double&
         ,msg='quasiparticle_check: particle density is not conserved',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check that pos points to something
    call assert(associated(this%pos),msg='quasiparticle_check: pos is not associated',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check pos is dimension npt
    call assert(size(this%pos).EQ.this%npt,msg='quasiparticle_check: pos is not of size npt',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check that pos does not contain NANs
    call assert(all(this%pos.EQ.this%pos),msg='quasiparticle_check: pos contains NANs',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check that pos does not contain huge values
    call assert(all(abs(this%pos).LT.huge(1.0_double))&
         ,msg='quasiparticle_check: pos has huge values',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check that mom points to something
    call assert(associated(this%mom),msg='quasiparticle_check: mom is not associated',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check mom is dimension npt
    call assert(size(this%mom).EQ.this%npt,msg='quasiparticle_check: mom is not of size npt',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check that mom does not contain NANs
    call assert(all(this%mom.EQ.this%mom),msg='quasiparticle_check: mom contains NANs',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check that mom does not contain huge values
    call assert(all(abs(this%mom).LT.huge(1.0_double))&
         ,msg='quasiparticle_check: mom has huge values',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check that force points to something
    call assert(associated(this%force),msg='quasiparticle_check: force is not associated',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check force is dimension npt
    call assert(size(this%force).EQ.this%npt,msg='quasiparticle_check: force is not of size npt',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check that force does not contain NANs
    call assert(all(this%force.EQ.this%force),msg='quasiparticle_check: force contains NANs',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

    !check that force does not contain huge values
    call assert(all(abs(this%force).LT.huge(1.0_double))&
         ,msg='quasiparticle_check: force has huge values',iostat=quasiparticle_check)
    if(quasiparticle_check.NE.0)return

  end function quasiparticle_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the quasiparticle methods.
  !> \param[in] this is the quasiparticle object whose methods will be excercised.
  !> \return Nothing if all tests pass or 1 and a stop for the first failed test.
  !> \remark Will stop after first failed check.
  !======================================================================
  subroutine quasiparticle_test
    use testing_class
    use filemanager
    type(quasiparticle)::this
    character(len=label)::string
    integer(long)::unit
    
    !verify quasiparticle is compatible with current version
    include 'verification'

    !----- additional make tests -----
    write(*,*)'test quasiparticle make sets correct default values'
    call make(this)
    call assert(this%mass.EQ.1,msg='default mass is not 1')
    call assert(this%charge.EQ.0,msg='default charge is not 0')
    call assert(this%npt.EQ.1,msg='default npt is not 1')
    call assert(this%den(0).EQ.1,msg='default den is not 1')
    call assert(this%pos(0).EQ.0,msg='default pos is not 0')
    call assert(this%mom(0).EQ.0,msg='default mom is not 0')
    call assert(this%force(0).EQ.0,msg='default force is not 0')
    call kill(this)

    write(*,*)'test quasiparticle can be made with predefined npt'
    call make(this,npt=10)
    call assert(this%npt.EQ.10,msg='quasiparticle npt does not equal 10')
    call assert(check(this).EQ.0,msg='quasiparticle failed check upon make with npt=10')
    call kill(this)


    !----- additional kill tests -----
    write(*,*)'test quasiparticle kill method cleans up dynamic memory and pointers'
    call assert(.not.associated(this%den),msg='quasiparticle pointer den remains associated after killed.')
    call assert(.not.associated(this%pos),msg='quasiparticle pointer pos remains associated after killed.')
    call assert(.not.associated(this%mom),msg='quasiparticle pointer mom remains associated after killed.')
    call assert(.not.associated(this%force),msg='quasiparticle pointer force remains associated after killed.')


    !----- additional backup tests -----
    write(*,*)'test quasiparticle attributes are properly backupd in backupfile'
    call make(this,npt=2)
    !set non default attribute values
    this%mass=2._double
    this%charge=1._double
    this%den=0.5_double
    this%pos=1._double
    this%mom=1._double
    this%force=1._double
    call system('rm -f testquasiparticle.tmpfile')
    call backup(this,file='testquasiparticle.tmpfile')
    call kill(this)
    call make(this,file='testquasiparticle.tmpfile')
    !assert non defulat attribute values are conserved
    call assert(this%mass.EQ.2,msg='quasiparticle mass is not 1')
    call assert(this%charge.EQ.1,msg='quasiparticle charge is not 1')
    call assert(this%npt.EQ.2,msg='quasiparticle npt is not 2')
    call assert(all(this%den.EQ.0.5),msg='quasiparticle den is not 0.5')
    call assert(all(this%pos.EQ.1),msg='quasiparticle pos is not 1')
    call assert(all(this%mom.EQ.1),msg='quasiparticle mom is not 1')
    call assert(all(this%force.EQ.1),msg='quasiparticle force is not 1')
    call kill(this)    

    
    !----- additional update tests -----

    
    !----- additional reset tests -----
    
    !----- additional status tests -----
    write(*,*)'test quasiparticle status format'
    call make(this)
    call assert(trim(status(this)).EQ.' 0.10000000E+01 0.00000000E+00 0.00000000E+00 0.00000000E+00'&
         ,msg='quasiparticle status not properly formatted.')
    call kill(this)

    write(*,*)'ALL quasiparticle TESTS PASSED!'
  end subroutine quasiparticle_test
  !-----------------------------------------

end module quasiparticle_class

