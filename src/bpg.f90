!>\brief
!! Class bpg
!!\details
!! Bipartite graph with nodes of one partition (top) belonging to its
!! defining layer and the nodes of the other partition (bottom) beloging
!! to a source layer.
!! When a bpg is linked (stacked on top) of this bpg, the top layer of this
!! bpg will become the bottom layer of the stacked bpg. The stacked bpg top
!! layer nodes are referred this bpg's marks whose derivatives are used
!! for backpropagation.
!<------------------------------------------------------------------------
module BPG_class
  use type_kinds
  use layer_class
  implicit none
  private

  public::BPG, BPG_test
  public::make, kill, display, store, update, reset, check, link

  type BPG
     logical::initialized=.false.
     type(layer)::layer

     integer(long)::nsource,nmark
     real(double),pointer::zerothsource
     type(realptr),dimension(:),pointer::source
     real(double),dimension(:,:),pointer::W
     real(double),dimension(:),pointer::error
     type(realptr),dimension(:),pointer::markerror
     type(realptr),dimension(:,:),pointer::WT

  end type BPG

  !> Creates the BPG object.
  interface make
     module procedure BPG_init
  end interface

  !> Destroys the BPG object.
  interface kill
     module procedure BPG_kill
  end interface

  !> Displays the current state of the BPG object.
  interface display
     module procedure BPG_display
  end interface

  !> Stores the current state of the BPG object.
  interface store
     module procedure BPG_store
  end interface

  !> Recaluclates the BPG object.
  interface update
     module procedure BPG_update
  end interface

  !> Reinitializes the BPG object.
  interface reset
     module procedure BPG_reset
  end interface

  !> Checks that the BPG object.
  interface check
     module procedure BPG_check
  end interface

  !> Links that the BPG object.
  interface link
     module procedure BPG_link
     module procedure BPG_BPG_link
  end interface

contains

  !======================================================================
  !> \brief Links the BPG object to a prior BPG.
  !> \param this is the BPG object.
  !> \param sourceBPG is the prior BPG object.
  !=====================================================================
  subroutine BPG_BPG_link(this,sourcebpg)
    use rand_class
    type(BPG),intent(inout)::this
    type(bpg),intent(inout)::sourcebpg
    type(realptr),dimension(:),pointer::tmpsource
    real(double),dimension(:,:),allocatable::tmpW
    type(realptr),dimension(:),pointer::tmpmarkerror
    type(realptr),dimension(:,:),pointer::tmpWT
    integer(long)::i,j,N,M,L

!    !check sourcebpg
!    if(check(sourcebpg).NE.0)return

    !get current size of sources and sourcebpg nodes
    N=size(this%source)
    M=size(sourcebpg%layer%node)
    L=size(this%layer%node)

    !allocate temp pointer array with size equal to the sum of the current source and new sourcebpg
    if(associated(tmpsource))nullify(tmpsource)
    allocate(tmpsource(0:N+M-1))

    !copy source pointer targets
    do i=0,N-1
       tmpsource(i)%ptr=>this%source(i)%ptr
    end do
    do i=1,M
       tmpsource(i+N-1)%ptr=>sourcebpg%layer%node(i)
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

    !temp store weight matrix
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

    !source bpg stuff below

    !get current number of sourcebpg marks and new marks (this bpg nodes)
    N=0
    if(associated(sourcebpg%markerror))N=size(sourcebpg%markerror)
    M=size(this%layer%node)
    L=size(sourcebpg%layer%node)

    !allocate temp pointer array for sourcebpg markerrors that point to this error array 
    if(associated(tmpmarkerror))nullify(tmpmarkerror)
    allocate(tmpmarkerror(N+M))

    !copy markerror pointer targets
    i=0
    do while(i.LT.N)
       i=i+1
       tmpmarkerror(i)%ptr=>sourcebpg%markerror(i)%ptr
    end do
    do i=N+1,M
       tmpmarkerror(i)%ptr=>this%error(i)
    end do

    !nullify markerror pointer array and allocate a new one with expanded size
    nullify(sourcebpg%markerror)
    allocate(sourcebpg%markerror(N+M))

    !copy markerror pointer targets
    do i=1,N+M
       sourcebpg%markerror(i)%ptr => tmpmarkerror(i)%ptr
    end do

    !update mark counter
    sourcebpg%nmark=N+M

    !nullify tmpsource
    nullify(tmpmarkerror)

    !create temp markweight matrix
    if(associated(tmpWT))nullify(tmpWT)

    !temp store markweight matrix targets if any
    if(N.GT.0)then
       allocate(tmpWT(N,L))
       do i=1,N
          do j=1,L
             tmpWT(i,j)%ptr=>sourcebpg%WT(i,j)%ptr
          end do
       end do
    end if

    !resize markweight matrix
    nullify(sourcebpg%WT)
    allocate(sourcebpg%WT(N+M,L))

    !return original markweights
    do j=1,L
       i=0
       do while (i.LT.N)
          i=i+1
          sourcebpg%WT(i,j)%ptr=>tmpWT(i,j)%ptr
       end do
       do i=1,M !new marks
          sourcebpg%WT(N+i,j)%ptr=>this%W(i,(this%nsource-1)-L+j)
       end do
    end do

    !nullify temp markweights
    if(associated(tmpWT))nullify(tmpWT)

  end subroutine BPG_BPG_link
  !======================================================================
  !> \brief Links the BPG object to a prior sourcelayer.
  !> \param this is the BPG object.
  !> \param sourcelayer is a layer object.
  !=====================================================================
  subroutine BPG_link(this,sourcelayer)
    use rand_class
    type(BPG),intent(inout)::this
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

    !temp store weight matrix
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

  end subroutine BPG_link

  !======================================================================
  !> \brief Creates and initializes the BPG object.
  !> \param this is the BPG object to be initialized.
  !> \param[in] file is an optional string containing the name of a previously stored BPG file.
!!$  !> \remark If no input file is provided the user must manually initialize THIS using stout.
  !=====================================================================
  subroutine BPG_init(this,N,activation)
    use rand_class
    type(BPG),intent(inout)::this
    integer(long),intent(in)::N
    character(len=*),optional,intent(in)::activation

    integer::i

    if(present(activation))then
       call make(this%layer,N=N,activation=activation)
    else
       call make(this%layer,N=N)
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

  end subroutine BPG_init

  !======================================================================
  !> \brief Destroys the BPG object.
  !> \param this is the BPG object to be destroyed.
  !====================================================================
  subroutine BPG_kill(this)
    type(BPG),intent(inout)::this
 
    !kill the layer primitive
    call kill(this%layer)

    if(associated(this%source))nullify(this%source)
    if(associated(this%zerothsource))nullify(this%zerothsource)
    if(associated(this%W))nullify(this%W)
    if(associated(this%error))nullify(this%error)
    if(associated(this%markerror))nullify(this%markerror)
    if(associated(this%WT))nullify(this%WT)


    !un-initialize bpg object
    this%initialized=.false.

  end subroutine BPG_kill

  !======================================================================
  !> \brief Computes the current state of BPG object.
  !> \param this is the BPG  object to be updated.
  !======================================================================
  subroutine BPG_update(this,derivative)
    type(BPG),intent(inout)::this
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

    !update bpg layer input according to weight
    this%layer%input(:)=matmul(this%W(:,:),source(:))

    !clean up copy of source state
    if(allocated(source))deallocate(source)

    !update bpg layer state according to input
    call update(this%layer,df)

  end subroutine BPG_update

  !======================================================================
  !> \brief Re-initiallizes the BPG object.
  !> \param this is the BPG  object to be re-initialized.
  !======================================================================
  subroutine BPG_reset(this)
    type(BPG),intent(inout)::this
!!$    integer(long)::istate,jstate
!!$
!!$    call Note('Begin BPG_reset.')
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
!!$    if(check(this).EQ.1)call stop('BPG_reset: failed final check!')

  end subroutine BPG_reset

  !======================================================================
  !> \brief Stores the current state of the BPG object to file.
  !> \param[in] this is the BPG  object to be updated.
  !> \param[in] file is a string containing the location of the store file.
  !======================================================================
  subroutine BPG_store(this,file)
    type(BPG),intent(in)::this
    character*(*),intent(in)::file

!!$    integer(short)::unit
!!$    logical::usedunit      
!!$
!!$    call note('Begin BPG_store.')
!!$    call Note('input file= '//file)
!!$    if(check(this).NE.0)then
!!$       call warn('BPG_store: failed check.','not saving object.')
!!$    else
!!$
!!$       !assign a unique unit label
!!$       unit=newunit()
!!$
!!$       !open store file
!!$       open(unit,file=file)
!!$
!!$       !always write the data type on the first line
!!$       write(unit,*)'BPG'
!!$
!!$       !store the layer
!!$       call store(this%layer,file//'.layer')
!!$
!!$       !write the location of the layer
!!$       write(unit,*)quote(file//'.layer')
!!$
!!$       !******      Store below all the derived type's attributes       ******!
!!$       !******         in the order the MAKE command reads them         ******!
!!$
!!$
!!$
!!$
!!$
!!$       !*********************************************************************!
!!$       !=====================================================================!
!!$       !******      Example - Store an attribute called 'var  '    ***********!
!!$       ! write(unit,*)this%var                                               !
!!$       !*********************************************************************!
!!$       !=====================================================================!
!!$       !***  Example - Store an NxM matrix attribute called 'matrix'  ********!
!!$       ! write(unit,*)((this%matrix(i,j),j=1,M),i=1,N)                       !
!!$       !*********************************************************************!
!!$
!!$
!!$       !finished saving all attributes - now close store file
!!$       close(unit)
!!$    end if
  end subroutine BPG_store


  !======================================================================
  !> \brief Displays the BPG object.
  !> \param[in] this is the BPG object.
  !> \param[in] msg is an optional string message to preface the displayed object.
  !======================================================================
  subroutine BPG_display(this,msg)
    type(BPG),intent(in)::this
    character*(*),intent(in),optional::msg

!!$    call Note('Begin BPG_display.')
!!$
!!$    if(check(this).NE.0)then
!!$       call warn('BPG_display: failed check','displaying nothing.')
!!$       return
!!$    end if
!!$
!!$    write(Dunit,*)'____________________________________________________'
!!$    write(Dunit,*)'-------------------   BPG   -------------------'
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
!!$    call display(this%layer,msg='BPG layer')
!!$    write(Dunit,*)'===================================================='

  end subroutine BPG_display

  !======================================================================
  !> \brief Checks the BPG object.
  !> \param[in] this is the BPG object to be checked.
  !> \return Nothing if all checks pass or 1 and a warn for the first failed check.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function BPG_check(this)
    use testing_class
    type(BPG),intent(in)::this

    integer(long)::i,j

    !initiate with no problems found 
    BPG_check=0

    !check that object is initialized
    call assert(this%initialized,msg='BPG_check: BPG object not initialized.',iostat=BPG_check)
    if(BPG_check.NE.0)return

    !check the layer
    call assert(check(this%layer).EQ.0,msg='BPG_check: failed layer check!',iostat=BPG_check)
    if(BPG_check.NE.0)return

    !check zerothsource
    call assert(associated(this%zerothsource),msg='BPG_check: zerothsource is not associated',iostat=BPG_check)
    if(BPG_check.NE.0)return
    call assert(this%zerothsource.EQ.1.0_double,msg='BPG_check: zerothsource is not 1.0',iostat=BPG_check)
    if(BPG_check.NE.0)return

    !check source
    call assert(associated(this%source),msg='BPG_check: source is not associated',iostat=BPG_check)
    if(BPG_check.NE.0)return
    call assert(size(this%source).GE.1,msg='BPG_check: source size is less than 1',iostat=BPG_check)
    if(BPG_check.NE.0)return
    do i=0,size(this%source)-1
       call assert(this%source(i)%ptr.EQ.this%source(i)%ptr&
            ,msg='BPG_check: source elements point to NaN values',iostat=BPG_check)
       if(BPG_check.NE.0)return
    end do

    !check source counter
    call assert(size(this%source).EQ.this%nsource,msg='BPG_check: source counter does not equal size of source pointer array',&
         iostat=BPG_check)
    if(BPG_check.NE.0)return

    !check weight matrix
    call assert(associated(this%W),msg='BPG_check: weight matrix W is not associated',iostat=BPG_check)
    if(BPG_check.NE.0)return
    call assert(size(this%W,1).EQ.this%layer%N&
         ,msg='BPG_check: W dim 1 size does not equal number of layer nodes',iostat=BPG_check)
    if(BPG_check.NE.0)return
    call assert(size(this%W,2).EQ.this%nsource&
         ,msg='BPG_check: W dim 1 size does not equal number of sources',iostat=BPG_check)
    if(BPG_check.NE.0)return
    call assert(all(this%W.EQ.this%W)&
         ,msg='BPG_check: Weight elements have NaN values',iostat=BPG_check)
    if(BPG_check.NE.0)return

    !check error vector
    call assert(associated(this%error),msg='BPG_check: error vector is not associated',iostat=BPG_check)
    if(BPG_check.NE.0)return
    call assert(size(this%error).EQ.this%layer%N+1&
         ,msg='BPG_check: error vector size does not equal number of layer nodes +1',iostat=BPG_check)
    if(BPG_check.NE.0)return
    call assert(all(this%error.EQ.this%error)&
         ,msg='BPG_check: error vector elements have NaN values',iostat=BPG_check)
    if(BPG_check.NE.0)return

    !check mark counter
    call assert(this%nmark.GE.0,msg='mark counter is less negative.',iostat=BPG_check)
    if(BPG_check.NE.0)return

    !if marks are present
    if(this%nmark.GT.0)then
       
       !check markerror
       call assert(associated(this%markerror),msg='mark error is not associated while mark counter is positive.',iostat=BPG_check)
       if(BPG_check.NE.0)return
       call assert(size(this%markerror).EQ.this%nmark,msg='BPG_check: mark counter does not equal size of markerror pointer array',&
            iostat=BPG_check)
       if(BPG_check.NE.0)return
       do i=1,this%nmark
          call assert(this%markerror(i)%ptr.EQ.this%markerror(i)%ptr&
               ,msg='BPG_check: markerror elements point to NaN values',iostat=BPG_check)
          if(BPG_check.NE.0)return
       end do

       !check markweight matrix
       call assert(size(this%WT,1).EQ.this%nmark&
            ,msg='BPG_check: WT dim 1 size does not equal number of marks',iostat=BPG_check)
       if(BPG_check.NE.0)return
       call assert(size(this%WT,2).EQ.this%layer%N&
            ,msg='BPG_check: WT dim 2 size does not equal number of layer nodes',iostat=BPG_check)
       if(BPG_check.NE.0)return
       do i=1,this%nmark
          do j=1,this%layer%N
             call assert(this%WT(i,j)%ptr.EQ.this%WT(i,j)%ptr&
                  ,msg='BPG_check: mark weight elements have NaN values',iostat=BPG_check)
          end do
       end do
       if(BPG_check.NE.0)return
    end if

  end function BPG_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the BPG methods.
  !> \param[in] this is the BPG object whose methods will be excercised.
  !> \return Nothing if all tests pass or 1 and a stop for the first failed test.
  !> \remark Will stop after first failed check.
  !======================================================================
  subroutine BPG_test
    use testing_class
    type(BPG)::this,sourcebpg,sourcebpg2,markbpg,markbpg2
    type(layer)::sourcelayer,sourcelayer2

    integer i,j

    write(*,*)'test bpg can be created with 4 nodes.'
    call make(this,N=4)
    call assert(check(this).EQ.0,msg='bpg object was not created properly.')
    call kill(this)

    write(*,*)'test bpg kill method sets initiallization flag to false.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.this%initialized,msg='bpg object remains initialized after killed.')

    write(*,*)'test bpg creates layer primitive properly with 4 nodes.'
    call make(this,N=4)
    call assert(check(this%layer).EQ.0,msg='bpg layer primitive was not created properly')
    call kill(this)

    write(*,*)'test bpg kill method sets layer primitive initiallization flag to false.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.this%layer%initialized,msg='bpg layer primitive object remains initialized after killed.')

    write(*,*)'test bpg make method allocates source array.'
    call make(this,N=4)
    call assert(associated(this%source),msg='bpg make mehtod doe snot allocate source array.')
    call kill(this)

    write(*,*)'test bpg make method allocates zerothsource.'
    call make(this,N=4)
    call assert(associated(this%zerothsource),msg='bpg zerothsource array is not associated after make.')
    call kill(this)

    write(*,*)'test bpg kill method nullifies zerothsource.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%zerothsource),msg='bpg zerothsource array remains associated after killed.')

    write(*,*)'test bpg make method allocates source array.'
    call make(this,N=4)
    call assert(associated(this%source),msg='bpg make method does not allocate source array.')
    call kill(this)

    write(*,*)'test bpg kill method nullifies source array.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%source),msg='bpg source array remains associated after killed.')

    write(*,*)'test source nodes can be added from a layer.'
    call make(sourcelayer,N=4)
    call make(this,N=2)
    call link(this,sourcelayer)
    call assert(size(this%source).EQ.4+1,msg='bpg source nodes from layer does not equal 4+1.')
    call kill(this)
    call kill(sourcelayer)

    write(*,*)'test source nodes can be added from a bpg.'
    call make(sourcebpg,N=4)
    call make(this,N=2)
    call link(this,sourcebpg)
    call assert(size(this%source).EQ.4+1,msg='bpg source nodes from bpg does not equal 4+1.')
    call kill(this)
    call kill(sourcelayer)

    write(*,*)'test changing sourcelayer nodes will change source nodes.'
    call make(sourcelayer,N=1)
    sourcelayer%node=0._double
    call make(this,N=2)
    call link(this,sourcelayer)
    sourcelayer%node=0.5_double
    call assert(this%source(1)%ptr.EQ.0.5,msg='bpg source pointer does not point to source node value.')
    call kill(this)
    call kill(sourcelayer)

    write(*,*)'test bpg make method allocates weight matrix.'
    call make(this,N=4)
    call assert(associated(this%W),msg='bpg weight matrix is not allocated after make.')
    call kill(this)

    write(*,*)'test bpg kill method nullifies weight matrix.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%W),msg='bpg weight matrix remains associated after killed.')

    write(*,*)'test bpg object passes check after link to layer method'
    call make(sourcelayer,N=1)
    call make(this,N=2)
    call link(this,sourcelayer)
    call assert(check(this).EQ.0,msg='bpg object failed check after link method.')
    call kill(this)
    call kill(sourcelayer)

    write(*,*)'test bpg object passes check after link to bpg method'
    call make(sourcebpg,N=1)
    call make(this,N=2)
    call link(this,sourcebpg)
    call assert(check(this).EQ.0,msg='bpg object failed check after link to bpg method.')
    call kill(this)
    call kill(sourcelayer)

    write(*,*)'test source bpg object passes check after link to bpg method'
    call make(sourcebpg,N=1)
    call make(this,N=2)
    call link(this,sourcebpg)
    call assert(check(sourcebpg).EQ.0,msg='source bpg object failed check after link to bpg method.')
    call kill(this)
    call kill(sourcelayer)

    write(*,*)'test original weights are preserved after link to layer method'
    call make(sourcelayer,N=1)
    call make(this,N=2)
    call link(this,sourcelayer)
    this%W=0.5_double
    call make(sourcelayer2,N=1)
    call link(this,sourcelayer2)
    call assert(this%W(1,1).EQ.0.5_double,msg='bpg original weights are not preserved after link to layer method.')
    call kill(this)
    call kill(sourcelayer)
    call kill(sourcelayer2)

    write(*,*)'test original weights are preserved after link to bpg method'
    call make(sourcebpg,N=1)
    call make(this,N=2)
    call link(this,sourcebpg)
    this%W=0.5_double
    call make(sourcebpg2,N=1)
    call link(this,sourcebpg2)
    call assert(this%W(1,1).EQ.0.5_double,msg='bpg original weights are not preserved after link to bpg method.')
    call kill(this)
    call kill(sourcebpg)
    call kill(sourcebpg2)

!!$    write(*,*)'test new weights are 0.0 after link to layer method'
!!$    call make(sourcelayer,N=1)
!!$    call make(this,N=2)
!!$    call link(this,sourcelayer)
!!$    this%W=0.5_double
!!$    call make(sourcelayer2,N=1)
!!$    call link(this,sourcelayer2)
!!$    call assert(this%W(1,2).EQ.0.0_double,msg='bpg new weights are 0.0 after link to layer method.')
!!$    call kill(this)
!!$    call kill(sourcelayer)
!!$    call kill(sourcelayer2)
!!$
!!$    write(*,*)'test new weights are 0.0 after link to bpg method'
!!$    call make(sourcebpg,N=1)
!!$    call make(this,N=2)
!!$    call link(this,sourcebpg)
!!$    this%W=0.5_double
!!$    call make(sourcebpg2,N=1)
!!$    call link(this,sourcebpg2)
!!$    call assert(this%W(1,2).EQ.0.0_double,msg='bpg new weights are 0.0 after link to bpg method.')
!!$    call kill(this)
!!$    call kill(sourcebpg)
!!$    call kill(sourcebpg2)

    write(*,*)'test bpg make method allocates error vector.'
    call make(this,N=4)
    call assert(associated(this%error),msg='bpg error vector is not associated after make.')
    call kill(this)

    write(*,*)'test bpg kill method nullifies error vector.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%error),msg='bpg error vector remains associated after killed.')

    write(*,*)'test bpg kill method nullifies markerror vector.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%markerror),msg='bpg markerror vector remains associated after killed.')

    write(*,*)'test sourcebpg markerror vector is allocated after link to bpg method.'
    call make(this,N=5)
    call make(markbpg,N=2)
    call link(markbpg,this)
    call assert(associated(this%markerror),msg='source bpg markerror vector not allocated after link.')
    call kill(markbpg)
    call kill(this)

    write(*,*)'test bpg link method reallocates source bpg markerror vector correctly.'
    call make(this,N=5)
    call make(markbpg,N=2)
    call link(markbpg,this)
    call assert(size(this%markerror).EQ.2,msg='size of markerror vector not equal to markbpg nodes after link.')
    call kill(markbpg)
    call kill(this)

    write(*,*)'test bpg link method correctly sets markerror pointer to markbpg error.'
    call make(this,N=2)
    call make(markbpg,N=3)
    markbpg%error=0.5
    call link(markbpg,this)
    call assert(this%markerror(3)%ptr.EQ.0.5_double,msg='bpg markerror pointer value does not equal markbpg error after link.')
    call kill(markbpg)
    call kill(this)

    write(*,*)'test bpg kill method nullifies WT matrix.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%WT),msg='bpg WT matrix remains associated after killed.')

    write(*,*)'test bpg link method allocates WT matrix.'
    call make(this,N=5)
    call make(markbpg,N=2)
    call link(markbpg,this)
    call assert(associated(this%WT),msg='bpg WT matrix is not allocated after link.')
    call kill(markbpg)
    call kill(this)

    write(*,*)'test bpg link method reallocates WT matrix correctly.'
    call make(this,N=5)
    call make(markbpg,N=2)
    call link(markbpg,this)
    call assert(size(this%WT,1).EQ.2,msg='size of WT matrix dim 1 not equal to markbpg nodes after link.')
    call assert(size(this%WT,2).EQ.5,msg='size of WT matrix dim 2 not equal to bpg nodes after link.')
    call kill(markbpg)
    call kill(this)

    write(*,*)'test bpg link method correctly sets points to WT.'
    call make(this,N=3)
    call make(markbpg,N=4)
    call link(markbpg,this)
    markbpg%W=0.5
    call assert(this%WT(4,3)%ptr.EQ.0.5_double,msg='WT matrix elements do not equal bpg Weights after link.')
    call kill(markbpg)
    call kill(this)

    write(*,*)'test bpg link method correctly sets points to WT given two sources.'
    call make(this,N=3)
    call make(sourcebpg,N=3)
    call make(markbpg,N=4)
    markbpg%W(:,:)=0.1_double
    call link(markbpg,sourcebpg)
    markbpg%W(:,1:3)=0.5_double
    call link(markbpg,this)
    markbpg%W(:,4:6)=0.2_double
!!$    do i=1,4
!!$       write(*,*)(sourcebpg%WT(i,j)%ptr,j=1,2)
!!$    end do
!!$    write(*,*)
!!$    do i=1,4
!!$       write(*,*)(this%WT(i,j)%ptr,j=1,3)
!!$    end do
!!$    write(*,*)
!!$    do i=1,4
!!$       write(*,*)(markbpg%W(i,j),j=0,5)
!!$    end do
    call assert(sourcebpg%WT(4,3)%ptr.EQ.0.5_double,msg='sourcebpg WT elements do not equal bpg Weights after link.')
    call assert(this%WT(4,3)%ptr.EQ.0.2_double,msg='bpg WT elements do not equal bpg Weights after link.')
    call kill(markbpg)
    call kill(sourcebpg)
    call kill(this)

    write(*,*)'test bpg link method correctly sets points to WT given two marks.'
    call make(this,N=3)
    call make(markbpg,N=2)
    call make(markbpg2,N=4)
    markbpg%W(:,:)=0.1_double
    markbpg2%W(:,:)=0.2_double
    call link(markbpg,this)
    markbpg%W(:,1:3)=0.3_double
    call link(markbpg2,this)
    markbpg2%W(:,1:3)=0.4_double
!!$    do i=1,6
!!$       write(*,*)(this%WT(i,j)%ptr,j=1,3)
!!$    end do
!!$    write(*,*)
!!$    do i=1,2
!!$       write(*,*)(markbpg%W(i,j),j=0,3)
!!$    end do
!!$    write(*,*)
!!$    do i=1,4
!!$       write(*,*)(markbpg2%W(i,j),j=0,3)
!!$    end do
    call assert(this%WT(2,1)%ptr.EQ.0.3_double,msg='bpg WT elements do not equal markbpg Weights after link.')
    call assert(this%WT(6,1)%ptr.EQ.0.4_double,msg='bpg WT elements do not equal markbpg2 Weights after link.')
    call kill(markbpg)
    call kill(markbpg2)
    call kill(this)

    write(*,*)'test update method sums previous layer input according weight'
    call make(sourcelayer,N=2)
    sourcelayer%node=3.0_double
    call make(this,N=2)
    call link(this,sourcelayer)
    this%W(:,:)=0._double
    this%W(:,1:2)=1._double
    call update(this)
    call assert(this%layer%input(1),6.0_double,msg='bpg update method does not update layer input correctly.')
    call kill(this)
    call kill(sourcelayer)

    write(*,*)'test bpg can be created with tanh activation and properly sets layer activation.'
    call make(this,N=2,activation='tanh')
    call assert(trim(this%layer%activation).EQ.'tanh',msg='Setting bpg activation does not set layer activation.')
    call kill(this)

    write(*,*)'test bpg can be updated with derivative flag and dnode is set properly'
    call make(this,N=2)
    this%layer%dnode=0.3_double
    call update(this,derivative=.true.)
    call assert(this%layer%dnode(1).EQ.1.0_double,msg='update with derivative flag does not update layer dnode properly.')
    call kill(this)

  end subroutine BPG_test
  !-----------------------------------------

end module BPG_class

