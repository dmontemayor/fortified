!>\brief
!! layer abstract class
!!\details
!! A layer is a network of N independent real valued nodes. Node states are
!! depend on the input signal I a 1xN vector where element I_k determines
!! the activity (state) of node k. The layer abstract class forms the basis
!! for all artificial neural network models.
 !!\author
!! Daniel Montemayor Sept 2015
!<------------------------------------------------------------------------
module layer_class
  use type_kinds
  use testing_class
  use functions
  implicit none
  private

  public::layer, layer_test
  public::make, kill, display, store, update, reset, check

  type layer
     logical::initialized=.false.
     !**********     Enter your derived type's attributes here     **********!
     integer(long)::N=huge(1_long)
     real(double),dimension(:),pointer::node
     real(double),dimension(:),pointer::input

     !***********************************************************************!
     !=======================================================================!
     !***   Here are some example attributes your derived type may have   ***!
     ! type(primitive)::primitive2                     !an extra primitive   !
     ! integer(short)::label                           !a short integer      !
     ! integer(long)::ndim                             !a long integer       !
     ! real(double)::var                               !a real variable      !
     ! complex(double)::zed                            !a complex variable   !
     ! real(double),dimension(:,:),pointer::matrix     !a real matrix        !
     ! complex(double),dimension(:,:),pointer::Zmatrix !a complex matrix     !
     !***********************************************************************!
  end type layer

  !> Creates the layer object.
  interface make
     module procedure layer_init
  end interface

  !> Destroys the layer object.
  interface kill
     module procedure layer_kill
  end interface

  !> Displays the current state of the layer object.
  interface display
     module procedure layer_display
  end interface

  !> Stores the current state of the layer object.
  interface store
     module procedure layer_store
  end interface

  !> Recaluclates the layer object.
  interface update
     module procedure layer_update
  end interface

  !> Reinitializes the layer object.
  interface reset
     module procedure layer_reset
  end interface

  !> Checks that the layer object.
  interface check
     module procedure layer_check
  end interface


contains

  !======================================================================
  !> \brief Creates and initializes the layer object.
  !> \param this is the layer object to be initialized.
  !> \param[in] file is an optional string containing the name of a previously stored layer file.
!!$  !> \remark If no input file is provided the user must manually initialize THIS using stout.
  !=====================================================================
  subroutine layer_init(this,N)
    type(layer),intent(inout)::this
    integer(long),intent(in)::N

    !declare initialization in progress
    this%initialized=.false.

    this%N=N !assign number of nodes from input
    if(N.LE.0)return

    if(associated(this%node))nullify(this%node) !cleanup up memory
    allocate(this%node(0:this%N))               !allocate nodes (zeroth node is for bias and is always 1)
    this%node(0)=1_double

    if(associated(this%input))nullify(this%input) !cleanup up memory
    allocate(this%input(this%N))                  !allocate input array

    !declare initialization complete
    this%initialized=.true.

  end subroutine layer_init

  !======================================================================
  !> \brief Destroys the layer object.
  !> \param this is the layer object to be destroyed.
  !====================================================================
  subroutine layer_kill(this)
    type(layer),intent(inout)::this
 
    if(associated(this%node))nullify(this%node)
    if(associated(this%input))nullify(this%input)

    !un-initialize object
    this%initialized=.false.

  end subroutine layer_kill

  !======================================================================
  !> \brief Computes the current state of layer object.
  !> \param this is the layer  object to be updated.
  !======================================================================
  subroutine layer_update(this)
    type(layer),intent(inout)::this

!!$    call Note('Begin layer_update.')
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
!!$    !if(check(this).EQ.1)call stop('layer_update: failed final check!')

  end subroutine layer_update

  !======================================================================
  !> \brief Re-initiallizes the layer object.
  !> \param this is the layer  object to be re-initialized.
  !======================================================================
  subroutine layer_reset(this)
    type(layer),intent(inout)::this
!!$    integer(long)::istate,jstate
!!$
!!$    call Note('Begin layer_reset.')
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
!!$    if(check(this).EQ.1)call stop('layer_reset: failed final check!')

  end subroutine layer_reset

  !======================================================================
  !> \brief Stores the current state of the layer object to file.
  !> \param[in] this is the layer  object to be updated.
  !> \param[in] file is a string containing the location of the store file.
  !======================================================================
  subroutine layer_store(this,file)
    type(layer),intent(in)::this
    character*(*),intent(in)::file

!!$    integer(short)::unit
!!$    logical::usedunit      
!!$
!!$    call note('Begin layer_store.')
!!$    call Note('input file= '//file)
!!$    if(check(this).NE.0)then
!!$       call warn('layer_store: failed check.','not saving object.')
!!$    else
!!$
!!$       !assign a unique unit label
!!$       unit=newunit()
!!$
!!$       !open store file
!!$       open(unit,file=file)
!!$
!!$       !always write the data type on the first line
!!$       write(unit,*)'layer'
!!$
!!$       !store the primitive
!!$       call store(this%primitive,file//'.primitive')
!!$
!!$       !write the location of the primitive
!!$       write(unit,*)quote(file//'.primitive')
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
  end subroutine layer_store


  !======================================================================
  !> \brief Displays the layer object.
  !> \param[in] this is the layer object.
  !> \param[in] msg is an optional string message to preface the displayed object.
  !======================================================================
  subroutine layer_display(this,msg)
    type(layer),intent(in)::this
    character*(*),intent(in),optional::msg

!!$    call Note('Begin layer_display.')
!!$
!!$    if(check(this).NE.0)then
!!$       call warn('layer_display: failed check','displaying nothing.')
!!$       return
!!$    end if
!!$
!!$    write(Dunit,*)'____________________________________________________'
!!$    write(Dunit,*)'-------------------   layer   -------------------'
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
!!$    call display(this%primitive,msg='layer primitive')
!!$    write(Dunit,*)'===================================================='

  end subroutine layer_display

  !======================================================================
  !> \brief Checks the layer object.
  !> \param[in] this is the layer object to be checked.
  !> \return Nothing if all checks pass or 1 and a warn for the first failed check.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function layer_check(this)
    type(layer),intent(in)::this
    integer(short)::errsum=0

    !initiate with no problems found 
    layer_check=0

    !check that layer object is initialized
    call assert(this%initialized,msg='layer_check: layer not initialized.',iostat=layer_check)
    if(layer_check.NE.0)return

    !********    Check all attributes are within acceptable values    *******!
    !check number of nodes N is NaN or Huge
    call assert(this%N,this%N,msg='layer_check: N is NaN or Huge',iostat=layer_check)
    if(layer_check.NE.0)return

    !check number of nodes is positive
    call assert(this%N.GT.0,msg='layer_check: N is not a positive integer',iostat=layer_check)
    if(layer_check.NE.0)return

    !check node array points to something
    call assert(associated(this%node),msg='layer_check: node memory not associated',iostat=layer_check)
    if(layer_check.NE.0)return

    !check node array is correct size
    call assert(size(this%node),this%N+1,msg='layer_check: node array size is not N+1',iostat=layer_check)
    if(layer_check.NE.0)return

    !check zeroth node array element is has value 1.0
    call assert(this%node(0).EQ.1_double,msg='layer_check: zeroth node array element is not 1.0',iostat=layer_check)
    if(layer_check.NE.0)return

    !check input array points to something
    call assert(associated(this%input),msg='layer_check: input memory not associated',iostat=layer_check)
    if(layer_check.NE.0)return

    !check input array is correct size
    call assert(size(this%input),this%N,msg='layer_check: input array size is not N',iostat=layer_check)
    if(layer_check.NE.0)return



  end function layer_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the layer methods.
  !> \param[in] this is the layer object whose methods will be excercised.
  !> \remark Will stop after first failed test
  !======================================================================
  subroutine layer_test
    type(layer)::this
    integer(short)::ierr

    write(*,*)'test layer can be created with 4 nodes.'
    call make(this,N=4)
    call assert(check(this).EQ.0,msg='layer object was not created properly.')
    call kill(this)

    write(*,*)'test layer kill method sets initiallization flag to false.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.this%initialized,msg='layer object remains initialized after killed.')

    write(*,*)'test layer kill method cleans up node memory.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%node),msg='kill layer method does not nullify node array pointer.')

    write(*,*)'test layer kill method cleans up input memory.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%input),msg='kill layer method does not nullify input array pointer.')

    write(*,*)'test layer remains uninitiallized when created with 0 nodes.'
    call make(this,N=0)
    call assert(.not.this%initialized,msg='layer object remains initialized with zero node number.')
    call kill(this)

    write(*,*)'test layer remains uninitiallized when created with -4 nodes.'
    call make(this,N=-4)
    call assert(.not.this%initialized,msg='layer object remains initialized with -4 node number.')
    call kill(this)

    write(*,*)'test layer node array remains unassociated when created with 0 nodes.'
    call make(this,N=0)
    call assert(.not.associated(this%node),msg='layer node array pointer is associated when created with 0 nodes.')
    call kill(this)

    write(*,*)'test layer input array remains unassociated when created with -4 nodes.'
    call make(this,N=-4)
    call assert(.not.associated(this%input),msg='layer input array pointer is associated when created with -4 nodes.')
    call kill(this)

    !test for node memory leaks with double make
    write(*,*)'Test that succesive makes destroy previously held data.'
    call make(this,N=5)
    this%node(5)=99.9
    call make(this,N=4)
    call assert(this%node(5).NE.99.9,msg='second make did not prevent accesing previously stored node data.')

  end subroutine layer_test
  !-----------------------------------------

end module layer_class

