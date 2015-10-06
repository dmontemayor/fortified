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
  public::make, kill, display, store, update, reset, check, link

  type FFN
     logical::initialized=.false.
     type(layer)::layer

     integer(long)::M
     real(double),pointer::zerothsource
     type(realptr),dimension(:),pointer::source
     real(double),dimension(:,:),pointer::W
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

  !> Stores the current state of the FFN object.
  interface store
     module procedure FFN_store
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
  end interface

contains

  !======================================================================
  !> \brief Links the FFN object to a prior sourcelayer.
  !> \param this is the FFN object.
  !> \param sourcelayer is a layer object.
  !=====================================================================
  subroutine FFN_link(this,sourcelayer)
    type(FFN),intent(inout)::this
    type(layer),intent(in)::sourcelayer
    type(realptr),dimension(:),pointer::tmpsource
    real(double),dimension(:,:),allocatable::tmpW
    integer(long)::N,M,i,L

    !check sourcelayer
    if(check(sourcelayer).NE.0)return

    !get sizes of current sources and new sourcelayer
    N=size(this%source)
    M=size(sourcelayer%node)

    !allocate temp pointer array with size equal to the sum of the current source and new sourcelayer
    if(associated(tmpsource))nullify(tmpsource)
    allocate(tmpsource(0:N+M-1))

    !copy targets
    do i=0,N-1
       tmpsource(i)%ptr=>this%source(i)%ptr
    end do
    do i=N,N+M-1
       tmpsource(i)%ptr=>sourcelayer%node(i)
    end do

    !nullify source pointer array and allocate a new one with expanded size
    nullify(this%source)
    allocate(this%source(0:N+M-1))

    !copy targets
    do i=0,N+M-1
       this%source(i)%ptr => tmpsource(i)%ptr
    end do

    !update source counter
    this%M=N+M

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
    this%W=0.0_double

    !return original weights
    this%W(:,0:N-1)=tmpW

    !nullify temp weights
    deallocate(tmpW)

  end subroutine FFN_link

  !======================================================================
  !> \brief Creates and initializes the FFN object.
  !> \param this is the FFN object to be initialized.
  !> \param[in] file is an optional string containing the name of a previously stored FFN file.
!!$  !> \remark If no input file is provided the user must manually initialize THIS using stout.
  !=====================================================================
  subroutine FFN_init(this,N)!,file)
    type(FFN),intent(inout)::this
    integer(long),intent(in)::N


    call make(this%layer,N=N)

    !create zeroth source target
    if(associated(this%zerothsource))nullify(this%zerothsource)
    allocate(this%zerothsource)
    this%zerothsource=1.0_double

    !start source list with zerothsource
    if(associated(this%source))nullify(this%source)
    allocate(this%source(0:0))
    this%source(0)%ptr => this%zerothsource

    !set source counter
    this%M=1

    !create weight matrix
    if(associated(this%W))nullify(this%W)
    allocate(this%W(size(this%layer%node),0:0))

    !declare initialization complete
    this%initialized=.true.




!!$    call Note('Begin FFN_init.')
!!$
!!$    !check if input file is present and valid
!!$    if(present(file))then
!!$
!!$       !check if file is there
!!$       if(check(file).EQ.1)call stop('FFN_init: cannot find input file '//file)
!!$
!!$       !assign a unique unit label
!!$       unit=newunit()
!!$
!!$       !open the file
!!$       open(unit,file=file)
!!$
!!$       !read the file type - should always be on the first line
!!$       read(unit,*)filetype
!!$       filetype=adjustl(filetype)
!!$
!!$       !check if input file is the right kind of file
!!$       if(trim(filetype).NE.'FFN')& 
!!$            call Stop('FFN_init: input file is not valid.')
!!$
!!$    end if
!!$    
!!$    ! prepare primative
!!$    if(present(file))then
!!$       read(unit,*)infile
!!$       infile=adjustl(infile)
!!$       call make(this%layer,trim(infile))
!!$    else
!!$       call make(this%layer)
!!$    end if
!!$
!!$    !******     Initiate all your derived type's attributes below      ******!
!!$    !************************************************************************!
!!$
!!$
!!$
!!$
!!$    !************************************************************************!
!!$    !========================================================================!
!!$    !********           Example - Setup attribute 'var'            **********!
!!$    !                                                                        !
!!$    ! if(present(file))then                                                  !
!!$    !    read(unit,*)this%var          !Read from input file if present      !
!!$    ! else                                                                   !
!!$    !    write(*,*)'Please enter VAR.' !manually enter information           !
!!$    !    read(*,*)this%var                                                   !
!!$    ! end if                                                                 !
!!$    !                                                                        !
!!$    !************************************************************************!
!!$    !========================================================================!
!!$    !********   Example - Setup a matrix attribute called 'matrix'   ********!
!!$    !                                                                        !
!!$    ! if(associated(this%matrix))nullify(this%matrix) !cleanup memory        !
!!$    ! allocate(this%matirx(N,M))                      !create an NxM matrix  !
!!$    !                                                                        !
!!$    ! !Now fill that matrix with data from the input file or manually        !
!!$    ! if(present(file))then                                                  !
!!$    !    read(unit,*)((this%matrix(i,j),j=1,M),i=1,N) !read M loop first     !
!!$    ! else                                                                   !
!!$    !    write(*,*)'Enter value of matrix element:'                          !
!!$    !    do i=1,N                                                            !
!!$    !       do j=1,M                                                         !
!!$    !          write(*,*)'index=',i,j      !prompt user the matrix index     !
!!$    !          read(*,*)this%matrix(i,j)   !read data manually from keyboard !
!!$    !       end do                                                           !
!!$    !    end do                                                              !
!!$    ! end if                                                                 !
!!$    !                                                                        !
!!$    !************************************************************************!
!!$
!!$
!!$    !finished reading data - now close input file
!!$    if(present(file))close(unit)
!!$
!!$    !declare that initialization is complete
!!$    this%initialized=.true.
!!$
!!$    !update (or reset if brand new) this derived type
!!$    if(present(file))then
!!$       call update(this)
!!$    else
!!$       call reset(this)
!!$    end if
!!$
!!$    !do a final check before we exit
!!$    if(check(file).EQ.1)call stop('FFN_init: failed final check!')

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

!!$
!!$
!!$    !*******************       Nullify all pointers    **********************!
!!$
!!$ 
!!$
!!$
!!$
!!$    !************************************************************************!
!!$    !========================================================================!
!!$    !******        Example - cleanup matrix attribute 'matrix'       ********!
!!$    !                                                                        !
!!$    ! if(associated(this%matrix))nullify(this%matrix)                        !
!!$    !                                                                        !
!!$    !************************************************************************!
!!$
!!$
    !un-initialize ffn object
    this%initialized=.false.

  end subroutine FFN_kill

  !======================================================================
  !> \brief Computes the current state of FFN object.
  !> \param this is the FFN  object to be updated.
  !======================================================================
  subroutine FFN_update(this)
    type(FFN),intent(inout)::this

!!$    call Note('Begin FFN_update.')
!!$    
!!$    !Layers usually dont get updated
!!$    !    call update(this%layer)    
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
!!$    !                   of the layer's denisity                          !
!!$    !                                                                        !
!!$    ! this%var=0._double                                                     !
!!$    ! do istate=1,this%layer%nstate                                             !
!!$    !    this%var=this%var+this%layer%den(istate,istate)                        !
!!$    ! end do                                                                 !
!!$    !                                                                        !
!!$    !************************************************************************!
!!$
!!$    !Usually leave out final check before we exit
!!$    !if(check(this).EQ.1)call stop('FFN_update: failed final check!')

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
  !> \brief Stores the current state of the FFN object to file.
  !> \param[in] this is the FFN  object to be updated.
  !> \param[in] file is a string containing the location of the store file.
  !======================================================================
  subroutine FFN_store(this,file)
    type(FFN),intent(in)::this
    character*(*),intent(in)::file

!!$    integer(short)::unit
!!$    logical::usedunit      
!!$
!!$    call note('Begin FFN_store.')
!!$    call Note('input file= '//file)
!!$    if(check(this).NE.0)then
!!$       call warn('FFN_store: failed check.','not saving object.')
!!$    else
!!$
!!$       !assign a unique unit label
!!$       unit=newunit()
!!$
!!$       !open store file
!!$       open(unit,file=file)
!!$
!!$       !always write the data type on the first line
!!$       write(unit,*)'FFN'
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
  end subroutine FFN_store


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

    integer(long)::i

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
    call assert(size(this%source).EQ.this%M,msg='FFN_check: source counter M does not equal size of source pointer array',&
         iostat=FFN_check)
    if(FFN_check.NE.0)return

    !check weight matrix
    call assert(associated(this%W),msg='FFN_check: weight matrix W is not associated',iostat=FFN_check)
    if(FFN_check.NE.0)return
    call assert(size(this%W,1).EQ.size(this%layer%node)&
         ,msg='FFN_check: W dim 1 size does not equal number of layer nodes',iostat=FFN_check)
    if(FFN_check.NE.0)return
    call assert(size(this%W,2).EQ.size(this%source)&
         ,msg='FFN_check: W dim 1 size does not equal number of sources',iostat=FFN_check)
    if(FFN_check.NE.0)return
    call assert(all(this%W.EQ.this%W)&
         ,msg='FFN_check: Weight elements have NaN values',iostat=FFN_check)
    if(FFN_check.NE.0)return

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
    type(FFN)::this
    type(layer)::sourcelayer,sourcelayer2

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

    write(*,*)'test ffn kill method nullifies source array.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%source),msg='ffn source array remains associated after killed.')

    write(*,*)'test ffn kill method nullifies zerothsource.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%zerothsource),msg='ffn zerothsource array remains associated after killed.')

    write(*,*)'test source layer nodes can be added.'
    call make(sourcelayer,N=4)
    call make(this,N=2)
    call link(this,sourcelayer)
    call assert(size(this%source).EQ.4+1,msg='ffn source layer nodes does not equal 4+1.')
    call kill(this)
    call kill(sourcelayer)

    write(*,*)'test changing sourcelayer nodes will change source nodes.'
    call make(sourcelayer,N=1)
    sourcelayer%node=0._double
    call make(this,N=2)
    call link(this,sourcelayer)
    sourcelayer%node=0.5_double
    call assert(this%source(1)%ptr.EQ.0.5,msg='ffn source node do not equal source layer nodes.')
    call kill(this)
    call kill(sourcelayer)

    write(*,*)'test ffn kill method nullifies weight matrix.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%W),msg='ffn weight matrix remains associated after killed.')

    write(*,*)'test ffn object passes check after link method'
    call make(sourcelayer,N=1)
    call make(this,N=2)
    call link(this,sourcelayer)
    call assert(check(this).EQ.0,msg='ffn object failed check after link method.')
    call kill(this)
    call kill(sourcelayer)

    write(*,*)'test original weights are preserved after link method'
    call make(sourcelayer,N=1)
    call make(this,N=2)
    call link(this,sourcelayer)
    this%W=0.5_double
    call make(sourcelayer2,N=1)
    call link(this,sourcelayer)
    call assert(this%W(1,1).EQ.0.5_double,msg='ffn original weights are not preserved after link.')
    call kill(this)
    call kill(sourcelayer)
    call kill(sourcelayer2)

    write(*,*)'test new weights are 0.0 after link method'
    call make(sourcelayer,N=1)
    call make(this,N=2)
    call link(this,sourcelayer)
    this%W=0.5_double
    call make(sourcelayer2,N=1)
    call link(this,sourcelayer)
    call assert(this%W(1,2).EQ.0.0_double,msg='ffn new weights are 0.0 after link.')
    call kill(this)
    call kill(sourcelayer)
    call kill(sourcelayer2)

  end subroutine FFN_test
  !-----------------------------------------

end module FFN_class

