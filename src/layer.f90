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
  implicit none
  private

  public::layer, layer_test
  public::make, kill, display, backup, update, reset, check

  type layer
     logical::initialized=.false.
     character(len=label)::name='layer'
     integer(long)::N=huge(1_long)
     character(len=label)::activation
     real(double),dimension(:),pointer::node,dnode
     real(double),dimension(:),pointer::input
     !for future: RMBs and RNNs
     !real(double),dimension(:),pointer::mirror,dmirror
     !real(double),dimension(:),pointer::minput
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

  !> Backups the current state of the layer object.
  interface backup
     module procedure layer_backup
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
  !> \param[in] N is the number of nodes in layer.
  !> \param[in] activation is an optional string determining layer activation function.
  !> \remark If no activation is provided linear activity is assumed.
  !=====================================================================
  subroutine layer_init(this,N,activation,file)
    type(layer),intent(inout)::this
    integer(long),optional,intent(in)::N
    character(len=*),optional,intent(in)::activation
    character*(*),intent(in),optional::file

    !declare initialization in progress
    this%initialized=.false.

    !declare layer size
    if(present(N))then
       if(N.LE.0)return
       this%N=N !assign number of nodes from input
    else
       this%N=1 !assign default number of nodes
    end if

    if(associated(this%node))nullify(this%node) !cleanup up memory
    allocate(this%node(this%N))               !allocate nodes

    if(associated(this%input))nullify(this%input) !cleanup up memory
    allocate(this%input(this%N))                  !allocate input array

    if(associated(this%dnode))nullify(this%dnode) !cleanup up memory
    allocate(this%dnode(this%N))               !allocate nodes

    this%activation='linear'!default activation
    if(present(activation))this%activation=activation

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
    if(associated(this%dnode))nullify(this%dnode)

    !un-initialize object
    this%initialized=.false.

  end subroutine layer_kill

  !======================================================================
  !> \brief Computes the current state of layer object.
  !> \param this is the layer  object to be updated.
  !> \param derivative is an optional boolean that when true will calcuate
  !> the derivative of the activation function wrt to input at current input
  !======================================================================
  subroutine layer_update(this,derivative)
    use math
    use functions 
    type(layer),intent(inout)::this
    logical,intent(in),optional::derivative

    !consider: ReLU=max(0,x), and softplus=ln(1+e^x),dsoftplu=logistic
    ! also leaky ReLU=max(x,ax) for 1>a>0
    ! alo noisy units f(x)=acitvation(x)+sigma*gran()

    select case(trim(this%activation))
    case ('linear')
       this%node(1:this%N)=this%input
       if(present(derivative).and.derivative)this%dnode=1.0_double
    case ('softplus')
       this%node(1:this%N)=softplus(this%input)
       if(present(derivative).and.derivative)this%dnode=logistic(this%input)
    case ('logistic')
       this%node(1:this%N)=logistic(this%input)
       if(present(derivative).and.derivative)this%dnode=dlogistic(this%input)
    case ('tanh')
       this%node(1:this%N)=tanh(this%input)
       if(present(derivative).and.derivative)this%dnode=dtanh(this%input)
    case ('gaussian')
       this%node(1:this%N)=gaussian(this%input)
       if(present(derivative).and.derivative)&
            this%dnode=-this%input*gaussian(this%input)
    case ('bernoulli')
       this%node(1:this%N)=bernoulli(this%input)
       if(present(derivative).and.derivative)this%dnode=dlogistic(this%input)
    case ('oscillator')
       this%node(1:this%N)=sin(this%input)
       if(present(derivative).and.derivative)this%dnode=cos(this%input)
    case ('poisson')
       this%node(1:this%N)=poisson(this%input)
       if(present(derivative).and.derivative)this%dnode=-exp(-this%input)
    case ('softmax')
       this%node(1:this%N)=softmax(this%input)
       if(present(derivative).and.derivative)&
            this%dnode=diag(dsoftmax(-this%input))
    case default
       this%node(1:this%N)=this%input
       if(present(derivative).and.derivative)this%dnode=1.0_double
    end select

  end subroutine layer_update

  !======================================================================
  !> \brief Re-initiallizes the layer object.
  !> \param this is the layer  object to be re-initialized.
  !======================================================================
  subroutine layer_reset(this)
    type(layer),intent(inout)::this

  end subroutine layer_reset

  !======================================================================
  !> \brief Backups the current state of the layer object to file.
  !> \param[in] this is the layer  object to be updated.
  !> \param[in] file is a string containing the location of the backup file.
  !======================================================================
  subroutine layer_backup(this,file)
    type(layer),intent(in)::this
    character*(*),intent(in)::file

!!$    integer(short)::unit
!!$    logical::usedunit      
!!$
!!$    call note('Begin layer_backup.')
!!$    call Note('input file= '//file)
!!$    if(check(this).NE.0)then
!!$       call warn('layer_backup: failed check.','not saving object.')
!!$    else
!!$
!!$       !assign a unique unit label
!!$       unit=newunit()
!!$
!!$       !open backup file
!!$       open(unit,file=file)
!!$
!!$       !always write the data type on the first line
!!$       write(unit,*)'layer'
!!$
!!$       !backup the primitive
!!$       call backup(this%primitive,file//'.primitive')
!!$
!!$       !write the location of the primitive
!!$       write(unit,*)quote(file//'.primitive')
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
  end subroutine layer_backup

  !======================================================================
  !> \brief Retrun the layer object as a single line record entry.
  !> \param[in] this is the layer object.
  !> \param[in] msg is an optional string message to annotate the displayed object.
  !======================================================================
  character(len=line) function layer_display(this,msg)
    type(layer),intent(in)::this
    character*(*),intent(in),optional::msg
    character(len=7)::FMT='(A10)'

    write(layer_display,FMT)'helloworld'

   
  end function layer_display

  !======================================================================
  !> \brief Checks the layer object.
  !> \param[in] this is the layer object to be checked.
  !> \return Nothing if all checks pass or 1 and a warn for the first failed check.
  !> \remark Will return after first failed check.
  !======================================================================
  integer(short)function layer_check(this)
    use testing_class
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
    call assert(size(this%node),this%N,msg='layer_check: node array size is not N',iostat=layer_check)
    if(layer_check.NE.0)return

    !check input array points to something
    call assert(associated(this%input),msg='layer_check: input memory not associated',iostat=layer_check)
    if(layer_check.NE.0)return

    !check input array is correct size
    call assert(size(this%input),this%N,msg='layer_check: input array size is not N',iostat=layer_check)
    if(layer_check.NE.0)return

    !check dnode array points to something
    call assert(associated(this%dnode),msg='layer_check: dnode memory not associated',iostat=layer_check)
    if(layer_check.NE.0)return

    !check dnode array is correct size
    call assert(size(this%dnode),this%N,msg='layer_check: dnode array size is not N',iostat=layer_check)
    if(layer_check.NE.0)return

  end function layer_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the layer methods.
  !> \param[in] this is the layer object whose methods will be excercised.
  !> \remark Will stop after first failed test
  !======================================================================
  subroutine layer_test
    use testing_class
    type(layer)::this
    integer(short)::ierr
    character(len=label)::string

    !verify layer is compatible with current version
    include 'verification'
    
    !additional make tests

    !additional kill tests

    !additional display tests

    !additional backup tests

    !additional update tests

    !additional reset tests

    
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

    write(*,*)'Test that succesive makes destroy previously held data.'
    call make(this,N=5)
    this%node(5)=99.9
    call make(this,N=4)
    call assert(this%node(4).NE.99.9,msg='second make did not prevent accesing previously backuped node data.')
    call kill(this)

    write(*,*)'Test that default activation sets node states equal to input.'
    call make(this,N=2)
    this%input(1)=0.12345
    call update(this)
    call assert(this%node(1).EQ.0.12345,msg='default activation does not set node state equal to input signal.')
    call kill(this)

    write(*,*)'Test layer object can be initiallized with linear activation'
    call make(this,N=2,activation='linear')
    call assert(check(this).EQ.0,msg='layer cannot be initiallized with linear activation.')
    call kill(this)

    write(*,*)'Test update layer with linear activation node state equals input=0.12345.'
    call make(this,N=2,activation='linear')
    this%input=0.12345
    call update(this)
    call assert(this%node(1).EQ.0.12345,msg='linear activation layer node does not equal input=0.12345')
    call kill(this)

    write(*,*)'Test update layer with logistic activation node=0.5 when input=0.0 .'
    call make(this,N=2,activation='logistic')
    this%input=0_double
    call update(this)
    call assert(this%node(1).EQ.0.5_double,msg='logistic activation layer node does not equal 0.5 when input=0.0')
    call kill(this)

    write(*,*)'Test update layer with logistic activation node=1/(1+e) when input=-1.0 .'
    call make(this,N=2,activation='logistic')
    this%input=-1_double
    call update(this)
    call assert(this%node(1).EQ.1/(1+exp(1.0_double)),&
         msg='logistic activation layer node does not equal 1/(1+e) when input=-1.0')
    call kill(this)

    write(*,*)'Test update layer with tanh activation node=0.0 when input=0.0 .'
    call make(this,N=2,activation='tanh')
    this%input=0.0_double
    call update(this)
    call assert(this%node(1).EQ.0.0_double,msg='tanh activation layer node does not equal 0.0 when input=0.0')
    call kill(this)

    write(*,*)'Test update layer with tanh activation node=tanh(-1) when input=-1.0 .'
    call make(this,N=2,activation='tanh')
    this%input=-1_double
    call update(this)
    call assert(this%node(1).EQ.tanh(-1.0_double),&
         msg='tanh activation layer node does not equal tanh(-1) when input=-1.0')
    call kill(this)

    write(*,*)'Test update layer with gaussian activation node=1.0 when input=0.0 .'
    call make(this,N=2,activation='gaussian')
    this%input=0.0_double
    call update(this)
    call assert(this%node(1).EQ.1.0_double,msg='gaussian activation layer node does not equal 1.0 when input=0.0')
    call kill(this)

    write(*,*)'Test update layer with gaussian activation node=exp(-0.5) when input=-1.0 .'
    call make(this,N=2,activation='gaussian')
    this%input=-1_double
    call update(this)
    call assert(this%node(1).EQ.exp(-0.5_double),&
         msg='gaussian activation layer node does not equal exp(-0.5) when input=-1.0')
    call kill(this)

    write(*,*)'test layer kill method cleans up dnode memory.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%dnode),msg='kill layer method does not nullify node array pointer.')

    write(*,*)'test update layer does not update dnode without derivative flag'
    call make (this,N=4)
    this%dnode=3.0_double
    this%input=0.5_double
    call update(this)
    call assert(all(this%dnode.EQ.3.0_double),msg='update layer has changed dnode without derivative flag.')
    call kill(this)

    write(*,*)'test update layer does update dnode with true derivative flag'
    call make (this,N=4)
    this%dnode=3.0_double
    this%input=0.5_double
    call update(this,derivative=.true.)
    call assert(all(this%dnode.EQ.1.0_double),msg='update layer has not changed dnode with true derivative flag.')
    call kill(this)

    !test layer update for all activation types
    !linear
    !softplus
    !logistic
    !tanh
    !gaussian
    !bernoulli
    !oscillator
    !poisson
    !softmax
    
    
  end subroutine layer_test
  !-----------------------------------------

end module layer_class

