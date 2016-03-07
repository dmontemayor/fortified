!>\brief
!! Class ffn
!!\details
!! Feed foward network FFN is a type of layer distinguished by its connecting
!! weights to external input nodes. 
!<------------------------------------------------------------------------
module FFN_class
  use type_kinds
  use layer_class
  implicit none
  private

  public::FFN, FFN_test
  public::make, kill, display, backup, update, reset, check, link

  type FFN
     logical::initialized=.false.
     character(len=label)::name='ffn'
     type(layer)::layer

     integer(long)::nsource,nmark
     real(double),pointer::zerothsource
     type(realptr),dimension(:),pointer::source
     real(double),dimension(:,:),pointer::W
     real(double),dimension(:),pointer::error
     type(realptr),dimension(:),pointer::markerror
     type(realptr),dimension(:,:),pointer::WT

  end type FFN

  !> Creates the FFN object.
  interface make
     module procedure FFN_init
  end interface

  !> Destroys the FFN object.
  interface kill
     module procedure FFN_kill
  end interface

  !> Displays the current state of the FFN object.
  interface display
     module procedure FFN_display
  end interface

  !> Backups the current state of the FFN object.
  interface backup
     module procedure FFN_backup
  end interface

  !> Recaluclates the FFN object.
  interface update
     module procedure FFN_update
  end interface

  !> Reinitializes the FFN object.
  interface reset
     module procedure FFN_reset
  end interface

  !> Checks that the FFN object.
  interface check
     module procedure FFN_check
  end interface

  !> Links that the FFN object.
  interface link
     module procedure FFN_link
     module procedure FFN_FFN_link
  end interface

contains

  !======================================================================
  !> \brief Links the FFN object to a prior FFN.
  !> \param this is the FFN object.
  !> \param sourceFFN is the prior FFN object.
  !=====================================================================
  subroutine FFN_FFN_link(this,sourceffn)
    use rand_class
    type(FFN),intent(inout)::this
    type(ffn),intent(inout)::sourceffn
    type(realptr),dimension(:),pointer::tmpsource
    real(double),dimension(:,:),allocatable::tmpW
    type(realptr),dimension(:),pointer::tmpmarkerror
    type(realptr),dimension(:,:),pointer::tmpWT
    integer(long)::i,j,N,M,L

!    !check sourceffn
!    if(check(sourceffn).NE.0)return

    !get current size of sources and sourceffn nodes
    N=size(this%source)
    M=size(sourceffn%layer%node)
    L=size(this%layer%node)

    !allocate temp pointer array with size equal to the sum of the current source and new sourceffn
    if(associated(tmpsource))nullify(tmpsource)
    allocate(tmpsource(0:N+M-1))

    !copy source pointer targets
    do i=0,N-1
       tmpsource(i)%ptr=>this%source(i)%ptr
    end do
    do i=1,M
       tmpsource(i+N-1)%ptr=>sourceffn%layer%node(i)
    end do

    !nullify source pointer array and allocate a new one with expanded size
    nullify(this%source)
    allocate(this%source(0:N+M-1))

    !copy source pointer targets
    do i=0,N+M-1
       this%source(i)%ptr => tmpsource(i)%ptr
    end do

    !update source counter
    this%nsource=N+M

    !nullify tmpsource
    nullify(tmpsource)

    !create temp weight matrix
    if(allocated(tmpW))deallocate(tmpW)
    allocate(tmpW(L,0:N-1))

    !temp backup weight matrix
    tmpW=this%W

    !resize weight matrix
    nullify(this%W)
    allocate(this%W(L,0:N+M-1))
    !this%W=0.0_double
    do i=1,L
       do j=0,N+M-1
          this%W(i,j)=gran()
       end do
    end do

    !return original weights
    this%W(:,0:N-1)=tmpW

    !nullify temp weights
    deallocate(tmpW)

    !source ffn stuff below

    !get current number of sourceffn marks and new marks (this ffn nodes)
    N=0
    if(associated(sourceffn%markerror))N=size(sourceffn%markerror)
    M=size(this%layer%node)
    L=size(sourceffn%layer%node)

    !allocate temp pointer array for sourceffn markerrors that point to this error array 
    if(associated(tmpmarkerror))nullify(tmpmarkerror)
    allocate(tmpmarkerror(N+M))

    !copy markerror pointer targets
    i=0
    do while(i.LT.N)
       i=i+1
       tmpmarkerror(i)%ptr=>sourceffn%markerror(i)%ptr
    end do
    do i=N+1,M
       tmpmarkerror(i)%ptr=>this%error(i)
    end do

    !nullify markerror pointer array and allocate a new one with expanded size
    nullify(sourceffn%markerror)
    allocate(sourceffn%markerror(N+M))

    !copy markerror pointer targets
    do i=1,N+M
       sourceffn%markerror(i)%ptr => tmpmarkerror(i)%ptr
    end do

    !update mark counter
    sourceffn%nmark=N+M

    !nullify tmpsource
    nullify(tmpmarkerror)

    !create temp markweight matrix
    if(associated(tmpWT))nullify(tmpWT)

    !temp backup markweight matrix targets if any
    if(N.GT.0)then
       allocate(tmpWT(N,L))
       do i=1,N
          do j=1,L
             tmpWT(i,j)%ptr=>sourceffn%WT(i,j)%ptr
          end do
       end do
    end if

    !resize markweight matrix
    nullify(sourceffn%WT)
    allocate(sourceffn%WT(N+M,L))

    !return original markweights
    do j=1,L
       i=0
       do while (i.LT.N)
          i=i+1
          sourceffn%WT(i,j)%ptr=>tmpWT(i,j)%ptr
       end do
       do i=1,M !new marks
          sourceffn%WT(N+i,j)%ptr=>this%W(i,(this%nsource-1)-L+j)
       end do
    end do

    !nullify temp markweights
    if(associated(tmpWT))nullify(tmpWT)

  end subroutine FFN_FFN_link
  !======================================================================
  !> \brief Links the FFN object to a prior sourcelayer.
  !> \param this is the FFN object.
  !> \param sourcelayer is a layer object.
  !=====================================================================
  subroutine FFN_link(this,sourcelayer)
    use rand_class
    type(FFN),intent(inout)::this
    type(layer),intent(in)::sourcelayer
    type(realptr),dimension(:),pointer::tmpsource
    real(double),dimension(:,:),allocatable::tmpW
    integer(long)::i,j,N,M,L

    !check sourcelayer
    if(check(sourcelayer).NE.0)return

    !get sizes of current sources and new sourcelayer
    N=size(this%source)
    M=size(sourcelayer%node)

    !allocate temp pointer array with size equal to the sum of the current source and new sourcelayer
    if(associated(tmpsource))nullify(tmpsource)
    allocate(tmpsource(0:N+M-1))

    !copy pointer targets
    do i=0,N-1
       tmpsource(i)%ptr=>this%source(i)%ptr
    end do
    do i=1,M
       tmpsource(i+N-1)%ptr=>sourcelayer%node(i)
    end do

    !nullify source pointer array and allocate a new one with expanded size
    nullify(this%source)
    allocate(this%source(0:N+M-1))

    !copy pointer targets
    do i=0,N+M-1
       this%source(i)%ptr => tmpsource(i)%ptr
    end do

    !update source counter
    this%nsource=N+M

    !nullify tmpsource
    nullify(tmpsource)

    !create temp weight matrix
    L=size(this%layer%node)
    if(allocated(tmpW))deallocate(tmpW)
    allocate(tmpW(L,0:N-1))

    !temp backup weight matrix
    tmpW=this%W

    !resize weight matrix
    nullify(this%W)
    allocate(this%W(L,0:N+M-1))
    !this%W=0.0_double
    do i=1,L
       do j=0,N+M-1
          this%W(i,j)=gran()
       end do
    end do

    !return original weights
    this%W(:,0:N-1)=tmpW

    !nullify temp weights
    deallocate(tmpW)

  end subroutine FFN_link

  !======================================================================
  !> \brief Creates and initializes the FFN object.
  !> \param this is the FFN object to be initialized.
  !> \param[in] file is an optional string containing the name of a previously backupd FFN file.
!!$  !> \remark If no input file is provided the user must manually initialize THIS using stout.
  !=====================================================================
  subroutine FFN_init(this,N,activation,file)
    use rand_class
    type(FFN),intent(inout)::this
    integer(long),optional,intent(in)::N
    character(len=*),optional,intent(in)::activation
    character*(*),intent(in),optional::file

    integer::i,M

    !declare layer size
    M=1
    if(present(N))M=N
    
    if(present(activation))then
       call make(this%layer,N=M,activation=activation)
    else
       call make(this%layer,N=M)
    end if

    !create zeroth source target
    if(associated(this%zerothsource))nullify(this%zerothsource)
    allocate(this%zerothsource)
    this%zerothsource=1.0_double

    !start source list with zerothsource
    if(associated(this%source))nullify(this%source)
    allocate(this%source(0:0))
    this%source(0)%ptr => this%zerothsource

    !set source counter
    this%nsource=1

    !create weight matrix
    if(associated(this%W))nullify(this%W)
    allocate(this%W(size(this%layer%node),0:0))
    !this%W=0.0_double
    do i=1,size(this%layer%node)
       this%W(i,0)=gran()
    end do

    !create error vector
    if(associated(this%error))nullify(this%error)
    allocate(this%error(0:size(this%layer%node)))

    !clear markerror vector
    if(associated(this%markerror))nullify(this%markerror)

    !set mark counter
    this%nmark=0

    !clear markweight matrix
    if(associated(this%WT))nullify(this%WT)

    !declare initialization complete
    this%initialized=.true.

  end subroutine FFN_init

  !======================================================================
  !> \brief Destroys the FFN object.
  !> \param this is the FFN object to be destroyed.
  !====================================================================
  subroutine FFN_kill(this)
    type(FFN),intent(inout)::this
 
    !kill the layer primitive
    call kill(this%layer)

    if(associated(this%source))nullify(this%source)
    if(associated(this%zerothsource))nullify(this%zerothsource)
    if(associated(this%W))nullify(this%W)
    if(associated(this%error))nullify(this%error)
    if(associated(this%markerror))nullify(this%markerror)
    if(associated(this%WT))nullify(this%WT)


    !un-initialize ffn object
    this%initialized=.false.

  end subroutine FFN_kill

  !======================================================================
  !> \brief Computes the current state of FFN object.
  !> \param this is the FFN  object to be updated.
  !======================================================================
  subroutine FFN_update(this,derivative)
    type(FFN),intent(inout)::this
    real(double),allocatable::source(:)
    logical,intent(in),optional::derivative

    integer(long)::i
    logical::df

    df=.false.
    if(present(derivative))df=derivative

    !copy source vaules into real array so we can do matmul
    if(allocated(source))deallocate(source)
    allocate(source(0:this%nsource-1))
    do i=0,this%nsource-1
       source(i)=this%source(i)%ptr
    end do

    !update ffn layer input according to weight
    this%layer%input(:)=matmul(this%W(:,:),source(:))

    !clean up copy of source state
    if(allocated(source))deallocate(source)

    !update ffn layer state according to input
    call update(this%layer,df)

  end subroutine FFN_update

  !======================================================================
  !> \brief Re-initiallizes the FFN object.
  !> \param this is the FFN  object to be re-initialized.
  !======================================================================
  subroutine FFN_reset(this)
    type(FFN),intent(inout)::this
!!$    integer(long)::istate,jstate
!!$
!!$    call Note('Begin FFN_reset.')
!!$
!!$    !reset the layer
!!$    call reset(this%layer)
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
!!$    if(check(this).EQ.1)call stop('FFN_reset: failed final check!')

  end subroutine FFN_reset

  !======================================================================
  !> \brief Backups the current state of the FFN object to file.
  !> \param[in] this is the FFN  object to be updated.
  !> \param[in] file is a string containing the location of the backup file.
  !======================================================================
  subroutine FFN_backup(this,file)
    type(FFN),intent(in)::this
    character*(*),intent(in)::file

!!$    integer(short)::unit
!!$    logical::usedunit      
!!$
!!$    call note('Begin FFN_backup.')
!!$    call Note('input file= '//file)
!!$    if(check(this).NE.0)then
!!$       call warn('FFN_backup: failed check.','not saving object.')
!!$    else
!!$
!!$       !assign a unique unit label
!!$       unit=newunit()
!!$
!!$       !open backup file
!!$       open(unit,file=file)
!!$
!!$       !always write the data type on the first line
!!$       write(unit,*)'FFN'
!!$
!!$       !backup the layer
!!$       call backup(this%layer,file//'.layer')
!!$
!!$       !write the location of the layer
!!$       write(unit,*)quote(file//'.layer')
!!$
!!$       !******      Backup below all the derived type's attributes       ******!
!!$       !******         in the order the MAKE command reads them         ******!
!!$
!!$
!!$
!!$
!!$
!!$       !*********************************************************************!
!!$       !=====================================================================!
!!$       !******      Example - Backup an attribute called 'var  '    ***********!
!!$       ! write(unit,*)this%var                                               !
!!$       !*********************************************************************!
!!$       !=====================================================================!
!!$       !***  Example - Backup an NxM matrix attribute called 'matrix'  ********!
!!$       ! write(unit,*)((this%matrix(i,j),j=1,M),i=1,N)                       !
!!$       !*********************************************************************!
!!$
!!$
!!$       !finished saving all attributes - now close backup file
!!$       close(unit)
!!$    end if
  end subroutine FFN_backup


  !======================================================================
  !> \brief Displays the FFN object.
  !> \param[in] this is the FFN object.
  !> \param[in] msg is an optional string message to preface the displayed object.
  !======================================================================
  subroutine FFN_display(this,msg)
    type(FFN),intent(in)::this
    character*(*),intent(in),optional::msg

!!$    call Note('Begin FFN_display.')
!!$
!!$    if(check(this).NE.0)then
!!$       call warn('FFN_display: failed check','displaying nothing.')
!!$       return
!!$    end if
!!$
!!$    write(Dunit,*)'____________________________________________________'
!!$    write(Dunit,*)'-------------------   FFN   -------------------'
!!$    if(present(msg))write(Dunit,*)msg
!!$    write(Dunit,*)'____________________________________________________'
!!$    
!!$    
!!$    !****    Display the derived type's attributes here if you want   ****!
!!$    
!!$    
!!$    
!!$    !*********************************************************************!
!!$    !=====================================================================!
!!$    !*******        Example display scalar attribute 'var'      **********!
!!$    !                                                                     !
!!$    ! write(Dunit,*)'VAR=',this%var                                       !
!!$    !                                                                     !
!!$    !*********************************************************************!
!!$    !=====================================================================!
!!$    !*******    Example display NxM matrix attribute 'matrix'   **********!
!!$    !                                                                     !
!!$    ! write(Dunit,*)'MATRIX=',this%matrix   !a simple example             !
!!$    !                                                                     !
!!$    ! write(Dunit,*)'MATRIX='               !a better example:            !
!!$    ! do i=1,N                          !write each row on a new line     !
!!$    !    write(Dunit,*)i,(this%matrix(i,j),j=1,M)                         !
!!$    ! end do                                                              !
!!$    !*********************************************************************!
!!$
!!$    call display(this%layer,msg='FFN layer')
!!$    write(Dunit,*)'===================================================='

  end subroutine FFN_display

  !======================================================================
  !> \brief Checks the FFN object.
  !> \param[in] this is the FFN object to be checked.
  !> \return Nothing if all checks pass or 1 and a warn for the first failed check.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function FFN_check(this)
    use testing_class
    type(FFN),intent(in)::this

    integer(long)::i,j

    !initiate with no problems found 
    FFN_check=0

    !check that object is initialized
    call assert(this%initialized,msg='FFN_check: FFN object not initialized.',iostat=FFN_check)
    if(FFN_check.NE.0)return

    !check the layer
    call assert(check(this%layer).EQ.0,msg='FFN_check: failed layer check!',iostat=FFN_check)
    if(FFN_check.NE.0)return

    !check zerothsource
    call assert(associated(this%zerothsource),msg='FFN_check: zerothsource is not associated',iostat=FFN_check)
    if(FFN_check.NE.0)return
    call assert(this%zerothsource.EQ.1.0_double,msg='FFN_check: zerothsource is not 1.0',iostat=FFN_check)
    if(FFN_check.NE.0)return

    !check source
    call assert(associated(this%source),msg='FFN_check: source is not associated',iostat=FFN_check)
    if(FFN_check.NE.0)return
    call assert(size(this%source).GE.1,msg='FFN_check: source size is less than 1',iostat=FFN_check)
    if(FFN_check.NE.0)return
    do i=0,size(this%source)-1
       call assert(this%source(i)%ptr.EQ.this%source(i)%ptr&
            ,msg='FFN_check: source elements point to NaN values',iostat=FFN_check)
       if(FFN_check.NE.0)return
    end do

    !check source counter
    call assert(size(this%source).EQ.this%nsource,msg='FFN_check: source counter does not equal size of source pointer array',&
         iostat=FFN_check)
    if(FFN_check.NE.0)return

    !check weight matrix
    call assert(associated(this%W),msg='FFN_check: weight matrix W is not associated',iostat=FFN_check)
    if(FFN_check.NE.0)return
    call assert(size(this%W,1).EQ.this%layer%N&
         ,msg='FFN_check: W dim 1 size does not equal number of layer nodes',iostat=FFN_check)
    if(FFN_check.NE.0)return
    call assert(size(this%W,2).EQ.this%nsource&
         ,msg='FFN_check: W dim 1 size does not equal number of sources',iostat=FFN_check)
    if(FFN_check.NE.0)return
    call assert(all(this%W.EQ.this%W)&
         ,msg='FFN_check: Weight elements have NaN values',iostat=FFN_check)
    if(FFN_check.NE.0)return

    !check error vector
    call assert(associated(this%error),msg='FFN_check: error vector is not associated',iostat=FFN_check)
    if(FFN_check.NE.0)return
    call assert(size(this%error).EQ.this%layer%N+1&
         ,msg='FFN_check: error vector size does not equal number of layer nodes +1',iostat=FFN_check)
    if(FFN_check.NE.0)return
    call assert(all(this%error.EQ.this%error)&
         ,msg='FFN_check: error vector elements have NaN values',iostat=FFN_check)
    if(FFN_check.NE.0)return

    !check mark counter
    call assert(this%nmark.GE.0,msg='mark counter is less negative.',iostat=FFN_check)
    if(FFN_check.NE.0)return

    !if marks are present
    if(this%nmark.GT.0)then
       
       !check markerror
       call assert(associated(this%markerror),msg='mark error is not associated while mark counter is positive.',iostat=FFN_check)
       if(FFN_check.NE.0)return
       call assert(size(this%markerror).EQ.this%nmark,msg='FFN_check: mark counter does not equal size of markerror pointer array',&
            iostat=FFN_check)
       if(FFN_check.NE.0)return
       do i=1,this%nmark
          call assert(this%markerror(i)%ptr.EQ.this%markerror(i)%ptr&
               ,msg='FFN_check: markerror elements point to NaN values',iostat=FFN_check)
          if(FFN_check.NE.0)return
       end do

       !check markweight matrix
       call assert(size(this%WT,1).EQ.this%nmark&
            ,msg='FFN_check: WT dim 1 size does not equal number of marks',iostat=FFN_check)
       if(FFN_check.NE.0)return
       call assert(size(this%WT,2).EQ.this%layer%N&
            ,msg='FFN_check: WT dim 2 size does not equal number of layer nodes',iostat=FFN_check)
       if(FFN_check.NE.0)return
       do i=1,this%nmark
          do j=1,this%layer%N
             call assert(this%WT(i,j)%ptr.EQ.this%WT(i,j)%ptr&
                  ,msg='FFN_check: mark weight elements have NaN values',iostat=FFN_check)
          end do
       end do
       if(FFN_check.NE.0)return
    end if

  end function FFN_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the FFN methods.
  !> \param[in] this is the FFN object whose methods will be excercised.
  !> \return Nothing if all tests pass or 1 and a stop for the first failed test.
  !> \remark Will stop after first failed check.
  !======================================================================
  subroutine FFN_test
    use testing_class
    type(FFN)::this,sourceffn,sourceffn2,markffn,markffn2
    type(layer)::sourcelayer,sourcelayer2

    integer i,j
    character(len=label)::string

    !verify ffn is compatible with current version
    include 'verification'

    write(*,*)'test ffn can be created with 4 nodes.'
    call make(this,N=4)
    call assert(check(this).EQ.0,msg='ffn object was not created properly.')
    call kill(this)

    write(*,*)'test ffn kill method sets initiallization flag to false.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.this%initialized,msg='ffn object remains initialized after killed.')

    write(*,*)'test ffn creates layer primitive properly with 4 nodes.'
    call make(this,N=4)
    call assert(check(this%layer).EQ.0,msg='ffn layer primitive was not created properly')
    call kill(this)

    write(*,*)'test ffn kill method sets layer primitive initiallization flag to false.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.this%layer%initialized,msg='ffn layer primitive object remains initialized after killed.')

    write(*,*)'test ffn make method allocates source array.'
    call make(this,N=4)
    call assert(associated(this%source),msg='ffn make mehtod doe snot allocate source array.')
    call kill(this)

    write(*,*)'test ffn make method allocates zerothsource.'
    call make(this,N=4)
    call assert(associated(this%zerothsource),msg='ffn zerothsource array is not associated after make.')
    call kill(this)

    write(*,*)'test ffn kill method nullifies zerothsource.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%zerothsource),msg='ffn zerothsource array remains associated after killed.')

    write(*,*)'test ffn make method allocates source array.'
    call make(this,N=4)
    call assert(associated(this%source),msg='ffn make method does not allocate source array.')
    call kill(this)

    write(*,*)'test ffn kill method nullifies source array.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%source),msg='ffn source array remains associated after killed.')

    write(*,*)'test source nodes can be added from a layer.'
    call make(sourcelayer,N=4)
    call make(this,N=2)
    call link(this,sourcelayer)
    call assert(size(this%source).EQ.4+1,msg='ffn source nodes from layer does not equal 4+1.')
    call kill(this)
    call kill(sourcelayer)

    write(*,*)'test source nodes can be added from a ffn.'
    call make(sourceffn,N=4)
    call make(this,N=2)
    call link(this,sourceffn)
    call assert(size(this%source).EQ.4+1,msg='ffn source nodes from ffn does not equal 4+1.')
    call kill(this)
    call kill(sourcelayer)

    write(*,*)'test changing sourcelayer nodes will change source nodes.'
    call make(sourcelayer,N=1)
    sourcelayer%node=0._double
    call make(this,N=2)
    call link(this,sourcelayer)
    sourcelayer%node=0.5_double
    call assert(this%source(1)%ptr.EQ.0.5,msg='ffn source pointer does not point to source node value.')
    call kill(this)
    call kill(sourcelayer)

    write(*,*)'test ffn make method allocates weight matrix.'
    call make(this,N=4)
    call assert(associated(this%W),msg='ffn weight matrix is not allocated after make.')
    call kill(this)

    write(*,*)'test ffn kill method nullifies weight matrix.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%W),msg='ffn weight matrix remains associated after killed.')

    write(*,*)'test ffn object passes check after link to layer method'
    call make(sourcelayer,N=1)
    call make(this,N=2)
    call link(this,sourcelayer)
    call assert(check(this).EQ.0,msg='ffn object failed check after link method.')
    call kill(this)
    call kill(sourcelayer)

    write(*,*)'test ffn object passes check after link to ffn method'
    call make(sourceffn,N=1)
    call make(this,N=2)
    call link(this,sourceffn)
    call assert(check(this).EQ.0,msg='ffn object failed check after link to ffn method.')
    call kill(this)
    call kill(sourcelayer)

    write(*,*)'test source ffn object passes check after link to ffn method'
    call make(sourceffn,N=1)
    call make(this,N=2)
    call link(this,sourceffn)
    call assert(check(sourceffn).EQ.0,msg='source ffn object failed check after link to ffn method.')
    call kill(this)
    call kill(sourcelayer)

    write(*,*)'test original weights are preserved after link to layer method'
    call make(sourcelayer,N=1)
    call make(this,N=2)
    call link(this,sourcelayer)
    this%W=0.5_double
    call make(sourcelayer2,N=1)
    call link(this,sourcelayer2)
    call assert(this%W(1,1).EQ.0.5_double,msg='ffn original weights are not preserved after link to layer method.')
    call kill(this)
    call kill(sourcelayer)
    call kill(sourcelayer2)

    write(*,*)'test original weights are preserved after link to ffn method'
    call make(sourceffn,N=1)
    call make(this,N=2)
    call link(this,sourceffn)
    this%W=0.5_double
    call make(sourceffn2,N=1)
    call link(this,sourceffn2)
    call assert(this%W(1,1).EQ.0.5_double,msg='ffn original weights are not preserved after link to ffn method.')
    call kill(this)
    call kill(sourceffn)
    call kill(sourceffn2)

!!$    write(*,*)'test new weights are 0.0 after link to layer method'
!!$    call make(sourcelayer,N=1)
!!$    call make(this,N=2)
!!$    call link(this,sourcelayer)
!!$    this%W=0.5_double
!!$    call make(sourcelayer2,N=1)
!!$    call link(this,sourcelayer2)
!!$    call assert(this%W(1,2).EQ.0.0_double,msg='ffn new weights are 0.0 after link to layer method.')
!!$    call kill(this)
!!$    call kill(sourcelayer)
!!$    call kill(sourcelayer2)
!!$
!!$    write(*,*)'test new weights are 0.0 after link to ffn method'
!!$    call make(sourceffn,N=1)
!!$    call make(this,N=2)
!!$    call link(this,sourceffn)
!!$    this%W=0.5_double
!!$    call make(sourceffn2,N=1)
!!$    call link(this,sourceffn2)
!!$    call assert(this%W(1,2).EQ.0.0_double,msg='ffn new weights are 0.0 after link to ffn method.')
!!$    call kill(this)
!!$    call kill(sourceffn)
!!$    call kill(sourceffn2)

    write(*,*)'test ffn make method allocates error vector.'
    call make(this,N=4)
    call assert(associated(this%error),msg='ffn error vector is not associated after make.')
    call kill(this)

    write(*,*)'test ffn kill method nullifies error vector.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%error),msg='ffn error vector remains associated after killed.')

    write(*,*)'test ffn kill method nullifies markerror vector.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%markerror),msg='ffn markerror vector remains associated after killed.')

    write(*,*)'test sourceffn markerror vector is allocated after link to ffn method.'
    call make(this,N=5)
    call make(markffn,N=2)
    call link(markffn,this)
    call assert(associated(this%markerror),msg='source ffn markerror vector not allocated after link.')
    call kill(markffn)
    call kill(this)

    write(*,*)'test ffn link method reallocates source ffn markerror vector correctly.'
    call make(this,N=5)
    call make(markffn,N=2)
    call link(markffn,this)
    call assert(size(this%markerror).EQ.2,msg='size of markerror vector not equal to markffn nodes after link.')
    call kill(markffn)
    call kill(this)

    write(*,*)'test ffn link method correctly sets markerror pointer to markffn error.'
    call make(this,N=2)
    call make(markffn,N=3)
    markffn%error=0.5
    call link(markffn,this)
    call assert(this%markerror(3)%ptr.EQ.0.5_double,msg='ffn markerror pointer value does not equal markffn error after link.')
    call kill(markffn)
    call kill(this)

    write(*,*)'test ffn kill method nullifies WT matrix.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%WT),msg='ffn WT matrix remains associated after killed.')

    write(*,*)'test ffn link method allocates WT matrix.'
    call make(this,N=5)
    call make(markffn,N=2)
    call link(markffn,this)
    call assert(associated(this%WT),msg='ffn WT matrix is not allocated after link.')
    call kill(markffn)
    call kill(this)

    write(*,*)'test ffn link method reallocates WT matrix correctly.'
    call make(this,N=5)
    call make(markffn,N=2)
    call link(markffn,this)
    call assert(size(this%WT,1).EQ.2,msg='size of WT matrix dim 1 not equal to markffn nodes after link.')
    call assert(size(this%WT,2).EQ.5,msg='size of WT matrix dim 2 not equal to ffn nodes after link.')
    call kill(markffn)
    call kill(this)

    write(*,*)'test ffn link method correctly sets points to WT.'
    call make(this,N=3)
    call make(markffn,N=4)
    call link(markffn,this)
    markffn%W=0.5
    call assert(this%WT(4,3)%ptr.EQ.0.5_double,msg='WT matrix elements do not equal ffn Weights after link.')
    call kill(markffn)
    call kill(this)

    write(*,*)'test ffn link method correctly sets points to WT given two sources.'
    call make(this,N=3)
    call make(sourceffn,N=3)
    call make(markffn,N=4)
    markffn%W(:,:)=0.1_double
    call link(markffn,sourceffn)
    markffn%W(:,1:3)=0.5_double
    call link(markffn,this)
    markffn%W(:,4:6)=0.2_double
!!$    do i=1,4
!!$       write(*,*)(sourceffn%WT(i,j)%ptr,j=1,2)
!!$    end do
!!$    write(*,*)
!!$    do i=1,4
!!$       write(*,*)(this%WT(i,j)%ptr,j=1,3)
!!$    end do
!!$    write(*,*)
!!$    do i=1,4
!!$       write(*,*)(markffn%W(i,j),j=0,5)
!!$    end do
    call assert(sourceffn%WT(4,3)%ptr.EQ.0.5_double,msg='sourceffn WT elements do not equal ffn Weights after link.')
    call assert(this%WT(4,3)%ptr.EQ.0.2_double,msg='ffn WT elements do not equal ffn Weights after link.')
    call kill(markffn)
    call kill(sourceffn)
    call kill(this)

    write(*,*)'test ffn link method correctly sets points to WT given two marks.'
    call make(this,N=3)
    call make(markffn,N=2)
    call make(markffn2,N=4)
    markffn%W(:,:)=0.1_double
    markffn2%W(:,:)=0.2_double
    call link(markffn,this)
    markffn%W(:,1:3)=0.3_double
    call link(markffn2,this)
    markffn2%W(:,1:3)=0.4_double
!!$    do i=1,6
!!$       write(*,*)(this%WT(i,j)%ptr,j=1,3)
!!$    end do
!!$    write(*,*)
!!$    do i=1,2
!!$       write(*,*)(markffn%W(i,j),j=0,3)
!!$    end do
!!$    write(*,*)
!!$    do i=1,4
!!$       write(*,*)(markffn2%W(i,j),j=0,3)
!!$    end do
    call assert(this%WT(2,1)%ptr.EQ.0.3_double,msg='ffn WT elements do not equal markffn Weights after link.')
    call assert(this%WT(6,1)%ptr.EQ.0.4_double,msg='ffn WT elements do not equal markffn2 Weights after link.')
    call kill(markffn)
    call kill(markffn2)
    call kill(this)

    write(*,*)'test update method sums previous layer input according weight'
    call make(sourcelayer,N=2)
    sourcelayer%node=3.0_double
    call make(this,N=2)
    call link(this,sourcelayer)
    this%W(:,:)=0._double
    this%W(:,1:2)=1._double
    call update(this)
    call assert(this%layer%input(1),6.0_double,msg='ffn update method does not update layer input correctly.')
    call kill(this)
    call kill(sourcelayer)

    write(*,*)'test ffn can be created with tanh activation and properly sets layer activation.'
    call make(this,N=2,activation='tanh')
    call assert(trim(this%layer%activation).EQ.'tanh',msg='Setting ffn activation does not set layer activation.')
    call kill(this)

    write(*,*)'test ffn can be updated with derivative flag and dnode is set properly'
    call make(this,N=2)
    this%layer%dnode=0.3_double
    call update(this,derivative=.true.)
    call assert(this%layer%dnode(1).EQ.1.0_double,msg='update with derivative flag does not update layer dnode properly.')
    call kill(this)

  end subroutine FFN_test
  !-----------------------------------------

end module FFN_class

