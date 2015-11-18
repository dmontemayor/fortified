!>\brief
!! device class
!!\details
!! A linked list of ffn layers begining with a visible layer and followed by
!! a series of ffn objects. The device class can train the network by
!! backpropagation.
!<------------------------------------------------------------------------
module DEVICE_class
  use type_kinds
  use layer_class
  use ffn_class
  implicit none
  private

  public::DEVICE, DEVICE_test
  public::make, kill, display, store, update, reset, check
  public::link, forwardpass, backprop, measure

  type ffnptr
     type(ffn),pointer::ffn
     integer(long)::rank
     real(double),dimension(:,:),pointer::gW
  end type ffnptr

  type DEVICE
     logical::initialized=.false.
     type(layer)::vislayer
     integer(long)::nlayer=huge(1)
     integer(long)::rank=huge(1)
     type(ffnptr),dimension(:),pointer::hiddenlayer
     real(double)::learningrate
  end type DEVICE

  !> calculate device perfomance
  interface measure
     module procedure DEVICE_measure
  end interface

  !> calculate error gradient by backpropagation
  interface backprop
     module procedure DEVICE_backprop
  end interface

  !> propagate input signal up to change top layerstates
  interface forwardpass
     module procedure DEVICE_forwardpass
  end interface

  !> link ffn to DEVICE
  interface link
     module procedure DEVICE_link
  end interface

  !> Creates the DEVICE object.
  interface make
     module procedure DEVICE_init
  end interface

  !> Destroys the DEVICE object.
  interface kill
     module procedure DEVICE_kill
  end interface

  !> Displays the current state of the DEVICE object.
  interface display
     module procedure DEVICE_display
  end interface

  !> Stores the current state of the DEVICE object.
  interface store
     module procedure DEVICE_store
  end interface

  !> Recaluclates the DEVICE object.
  interface update
     module procedure DEVICE_update
  end interface

  !> Reinitializes the DEVICE object.
  interface reset
     module procedure DEVICE_reset
  end interface

  !> Checks that the DEVICE object.
  interface check
     module procedure DEVICE_check
  end interface

contains
  !======================================================================
  !> \brief calucate Device performance given given dataset
  !> \param this is the DEVICE object to be measured.
  !> \param[in] dataset is a realmatrix of size  Ninput + Noutput x Nsample
  !>  containing a batch of Nsample training data vetor pairs.
  !> \param[in] file is an optional string name of output file to save analysis
  !> \param[inout] mse is an optional real output for mean square error
  !> \param[inout] bce is an optional real output for binary cross entropy
  !!> \param[inout] xrange is an optional real output for dataset xrange 
  !!> \param[inout] xmin is an optional real output for dataset minimum x value
  !> \param[in] nmeasure is an optional integer to repeat the measurement
  !>            ,average measurement is returned
  !=====================================================================
  !subroutine DEVICE_measure(this,dataset,mse,bce,file,nmeasure,xrange,xmin)
  subroutine DEVICE_measure(this,dataset,mse,bce,file,nmeasure)
    use testing_class
    type(DEVICE),intent(inout)::this
    real(double),intent(in)::dataset(:,:)
    real(double),intent(inout),optional::mse,bce
    character(len=*),intent(in),optional::file
    integer(long),intent(in),optional::nmeasure
    !real(double),intent(inout),optional::xrange,xmin

    integer(short)::ierr
    integer(long)::i,j,nsample,nx,ny,nmeas
    real(double),allocatable::prob(:)
    !real(double)::xL,x0

    !setup io
    nx=this%vislayer%N
    do i=1,this%nlayer
       if(this%hiddenlayer(i)%ffn%nmark.EQ.0)&
            ny=this%hiddenlayer(i)%ffn%layer%N
    end do
    
    !check dataset
    call assert(size(dataset,1).EQ.nx+ny,msg='DEVICE_measure: dataset dim 1 is incongruent with device io',iostat=ierr)
    if(ierr.NE.0)return

    nmeas=1
    if(present(nmeasure).and.nmeasure.GT.0)nmeas=nmeasure
    nsample=size(dataset,2)

    if(present(mse))then
       !calculate mse
       mse=0._double
       do i=1,nsample
          !load input
          this%vislayer%input(1:nx)=dataset(1:nx,i)

          !repeat measurement
          do j=1,nmeas
             !forward pass input
             call forwardpass(this)
             !accumulate mse
             mse=mse+sum((dataset(nx+1:nx+ny,i)&
                  -this%hiddenlayer(this%nlayer)%ffn%layer%node(1:ny))**2)
          end do
       end do
       mse=mse/real(nsample*nmeas,double)
    end if

    if(present(bce))then
       !calculate bce
       if(allocated(prob))deallocate(prob)
       allocate(prob(ny))
       bce=0._double
       do i=1,nsample
          !load input
          this%vislayer%input(1:nx)=dataset(1:nx,i)

          !repeat measurement
          prob=0._double
          do j=1,nmeas
             !forward pass input
             call forwardpass(this)
             !accumulate prob
             prob=prob+this%hiddenlayer(this%nlayer)%ffn%layer%node(1:ny)
          end do
          prob=prob/real(nmeas,double)
          !accumulate bce
          bce=bce+sum(log(prob)*dataset(nx+1:nx+ny,i)&
               +log(1-prob)*(1-dataset(nx+1:nx+ny,i)))
       end do
       !normalize bce
       bce=bce/real(nsample*nmeas,double)
       if(allocated(prob))deallocate(prob)
    end if
    
    if(present(file))then
       !xL=1.0_double
       !if(present(xrange))xL=xrange
       !x0=0._double
       !if(present(xmin))x0=xmin
       open(2232,file=file)
       do i=1,nsample
          !load input
          this%vislayer%input(1:nx)=dataset(1:nx,i)
          !forward pass input
          call forwardpass(this)
          !dump output
          write(2232,*)i,this%vislayer%input(1:nx),this%hiddenlayer(this%nlayer)%ffn%layer%node(1:ny)
!!$          write(2232,*)i,this%vislayer%input(1:nx)*xL+x0,this%hiddenlayer(this%nlayer)%ffn%layer%node(1:ny)
       end do
       close(2232)
    end if
    
  end subroutine DEVICE_measure

  !======================================================================
  !> \brief back propagate error graident given dataset.
  !> \param this is the DEVICE object to be measured.
  !> \param[in] dataset is a realmatrix of size  Ninput + Noutput x Nsample
  !>  containing a batch of Nsample training data vetor pairs.
  !=====================================================================
  subroutine DEVICE_backprop(this,dataset)
    !use testing_class
    type(DEVICE),intent(inout)::this
    real(double),intent(in)::dataset(:,:)

    integer(short)::ierr
    integer(long)::isample,nsample,nx,ny,irank,ilayer,inode,imark,isource
    real(double)::zsum

    !!check device
    !call assert(check(this).EQ.0,msg='DEVICE_backprop: input device did not pass check',iostat=ierr)
    !if(ierr.NE.0)return

    !setup io
    nx=this%vislayer%N
    do ilayer=1,this%nlayer
       if(this%hiddenlayer(ilayer)%ffn%nmark.EQ.0)&
            ny=this%hiddenlayer(ilayer)%ffn%layer%N
    end do

    !!check dataset
    !call assert(size(dataset,1).EQ.nx+ny,msg='DEVICE_backprop: dataset dim 1 is incongruent with device io',iostat=ierr)
    !if(ierr.NE.0)return

    nsample=size(dataset,2)

    !initialize derivate weights
    do ilayer=1,this%nlayer
       this%hiddenlayer(ilayer)%gW=0._double
    end do

    !loop over samples
    do isample=1,nsample
       !load input
       this%vislayer%input(1:nx)=dataset(1:nx,isample)
       !forward pass input
       call forwardpass(this,derivative=.true.)

       !calculate gradient starting with highest ranked layers to lowest
       do irank=this%rank,1,-1
          do ilayer=1,this%nlayer
             if(this%hiddenlayer(ilayer)%rank.EQ.irank)then
                if(this%hiddenlayer(ilayer)%ffn%nmark.EQ.0)then

                   !ffn is an output layer
                   this%hiddenlayer(ilayer)%ffn%error(1:ny)=&
                        (dataset(nx+1:nx+ny,isample)-this%hiddenlayer(ilayer)%ffn%layer%node(1:ny))&
                        *this%hiddenlayer(ilayer)%ffn%layer%dnode(1:ny)
                else!ffn is an internal layer
                   do inode=1,this%hiddenlayer(ilayer)%ffn%layer%N

                      zsum=0.0_double
                      do imark=1,this%hiddenlayer(ilayer)%ffn%nmark
                         zsum=zsum+this%hiddenlayer(ilayer)%ffn%markerror(imark)%ptr&
                              *this%hiddenlayer(ilayer)%ffn%WT(imark,inode)%ptr
                      end do
                      this%hiddenlayer(ilayer)%ffn%error(inode)=this%hiddenlayer(ilayer)%ffn%error(inode)&
                           +zsum*this%hiddenlayer(ilayer)%ffn%layer%dnode(inode)
      
                   end do
                   
                end if
             end if
          end do
       end do

       !accumulate derivate weights
       do ilayer=1,this%nlayer
          do inode=1,this%hiddenlayer(ilayer)%ffn%layer%N
             do isource=0,this%hiddenlayer(ilayer)%ffn%nsource-1
                this%hiddenlayer(ilayer)%gW(inode,isource)=this%hiddenlayer(ilayer)%gW(inode,isource)&
                     +this%hiddenlayer(ilayer)%ffn%error(inode)&
                     *this%hiddenlayer(ilayer)%ffn%source(isource)%ptr
             end do
          end do
       end do
    end do

    !update weights for all hidden layers
    do ilayer=1,this%nlayer
       this%hiddenlayer(ilayer)%ffn%W=this%hiddenlayer(ilayer)%ffn%W&
            +this%learningrate*this%hiddenlayer(ilayer)%gW/real(nsample,double)
    end do


  end subroutine DEVICE_backprop
  !======================================================================
  !> \brief propagate any input signal up the DEVICE object.
  !> \param this is the DEVICE object to be updated.
  !=====================================================================
  subroutine DEVICE_forwardpass(this,derivative)
    type(DEVICE),intent(inout)::this
    logical,intent(in),optional::derivative

    logical::df
    integer(long)::i,j

    df=.false.
    if(present(derivative))df=derivative

    !update visible layer
    call update(this%vislayer,df)

    !update hiddenlayers from lowest rank to highest rank
    do i=1,this%rank
       do j=1,this%nlayer
          if(this%hiddenlayer(j)%rank.EQ.i)call update(this%hiddenlayer(j)%ffn,df)
       end do
    end do

  end subroutine DEVICE_forwardpass
 
 !======================================================================
  !> \brief links ffn to the DEVICE object.
  !> \param this is the DEVICE object to be initialized.
  !> \param ffn is an ffn to add to the device
  !=====================================================================
  subroutine DEVICE_link(this,newffn)
    type(DEVICE),intent(inout)::this
    type(ffn),intent(inout),target::newffn
    type(ffnptr),dimension(:),pointer::tmphiddenlayer

    integer(long)::i,maxrank=0

    if(this%nlayer.GT.0)then
       !create temphiddenlayer
       if(associated(tmphiddenlayer))nullify(tmphiddenlayer)
       allocate(tmphiddenlayer(this%nlayer))
       !copy hiddenlayer targets to temphiddenlayer and store maxrank
       do i=1,this%nlayer
          tmphiddenlayer(i)%ffn => this%hiddenlayer(i)%ffn
          tmphiddenlayer(i)%rank = this%hiddenlayer(i)%rank
          if(this%hiddenlayer(i)%rank.GE.maxrank)maxrank=this%hiddenlayer(i)%rank
       end do
    end if

    !reallocate hiddenlayer array to accomodate newffn
    nullify(this%hiddenlayer)
    allocate(this%hiddenlayer(this%nlayer+1))

    if(this%nlayer.GT.0)then
       !return original hiddenlayer targets and ranks
       do i=1,this%nlayer
          this%hiddenlayer(i)%ffn => tmphiddenlayer(i)%ffn
          this%hiddenlayer(i)%rank = tmphiddenlayer(i)%rank
       end do
       !nullify temp hidden layer array
       if(associated(tmphiddenlayer))nullify(tmphiddenlayer)
    end if

    !add newffn target to hiddenlayer pointer
    this%hiddenlayer(this%nlayer+1)%ffn => newffn

    if(this%nlayer.GT.0)then
       !link newffn to last ffn in hiddenlayer array
       call link(newffn,this%hiddenlayer(this%nlayer)%ffn)
       !set newffn rank to last ffn +1
       this%hiddenlayer(this%nlayer+1)%rank = this%hiddenlayer(this%nlayer)%rank+1
    else
       !link newffn to visible layer
       call link(newffn,this%vislayer)
       !set newffn rank to 1
       this%hiddenlayer(this%nlayer+1)%rank = 1
    end if

    !check if new ffn rank is greater than maxrank
    if(this%hiddenlayer(this%nlayer+1)%rank.GE.maxrank)maxrank=this%hiddenlayer(this%nlayer+1)%rank

    !increment layer counter
    this%nlayer=this%nlayer+1

    !set derivative weight matricies
    do i=1,this%nlayer
       if(associated(this%hiddenlayer(i)%gW))nullify(this%hiddenlayer(i)%gW)
       allocate(this%hiddenlayer(i)%gW(size(this%hiddenlayer(i)%ffn%W,1),0:size(this%hiddenlayer(i)%ffn%W,2)-1))
    end do

    !set device rank
    this%rank=maxrank

  end subroutine DEVICE_link

  !======================================================================
  !> \brief Creates and initializes the DEVICE object.
  !> \param this is the DEVICE object to be initialized.
  !> \param[in] file is an optional string containing the name of a previously stored DEVICE file.
  !=====================================================================
  subroutine DEVICE_init(this,N,activation,loss)
    type(DEVICE),intent(inout)::this
    integer(long),intent(in)::N
    character(len=*),optional,intent(in)::activation,loss

!softmax=
!square =1/2N sum_i^N(y_i-yhat_i)
!binarycrossentropy=-1/N sum_i^N(log(p_i)*yhat_i+log(1-p_i)*(1-yhat_i)


    if(present(activation))then
       call make(this%vislayer,N=N,activation=activation)
    else
       call make(this%vislayer,N=N)
    end if

    this%nlayer=0
    this%rank=0
    this%learningrate=0.001

    !declare initialization complete
    this%initialized=.true.
  end subroutine DEVICE_init

  !======================================================================
  !> \brief Destroys the DEVICE object.
  !> \param this is the DEVICE object to be destroyed.
  !====================================================================
  subroutine DEVICE_kill(this)
    type(DEVICE),intent(inout)::this
 
    call kill(this%vislayer)
    if(associated(this%hiddenlayer))nullify(this%hiddenlayer)

    !un-initialized DEVICE object
    this%initialized=.false.
  end subroutine DEVICE_kill

  !======================================================================
  !> \brief Computes the current state of DEVICE object.
  !> \param this is the DEVICE  object to be updated.
  !======================================================================
  subroutine DEVICE_update(this)
    type(DEVICE),intent(inout)::this


  end subroutine DEVICE_update

  !======================================================================
  !> \brief Re-initiallizes the DEVICE object.
  !> \param this is the DEVICE  object to be re-initialized.
  !======================================================================
  subroutine DEVICE_reset(this)
    type(DEVICE),intent(inout)::this
  end subroutine DEVICE_reset

  !======================================================================
  !> \brief Stores the current state of the DEVICE object to file.
  !> \param[in] this is the DEVICE  object to be updated.
  !> \param[in] file is a string containing the location of the store file.
  !======================================================================
  subroutine DEVICE_store(this,file)
    type(DEVICE),intent(in)::this
    character*(*),intent(in)::file

  end subroutine DEVICE_store


  !======================================================================
  !> \brief Displays the DEVICE object.
  !> \param[in] this is the DEVICE object.
  !> \param[in] msg is an optional string message to preface the displayed object.
  !======================================================================
  subroutine DEVICE_display(this,msg)
    type(DEVICE),intent(in)::this
    character*(*),intent(in),optional::msg

  end subroutine DEVICE_display

  !======================================================================
  !> \brief Checks the DEVICE object.
  !> \param[in] this is the DEVICE object to be checked.
  !> \return Nothing if all checks pass or 1 and a warn for the first failed check.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function DEVICE_check(this)
    use testing_class
    type(DEVICE),intent(in)::this

    integer(long)::i

    !initiate with no problems found 
    DEVICE_check=0
    !call Note('Checking DEVICE.')

    !check that object is initialized
    call assert(this%initialized,msg='DEVICE_check: DEVICE object not initialized.',iostat=DEVICE_check)
    if(DEVICE_check.NE.0)return

    !check vislayer
    call assert(check(this%vislayer).EQ.0,msg='DEVICE_check: vislayer failed check.',iostat=DEVICE_check)
    if(DEVICE_check.NE.0)return

    !check nlayer
    call assert(this%nlayer,this%nlayer,msg='DEVICE_check: nlayer has bad values.',iostat=DEVICE_check)
    if(DEVICE_check.NE.0)return
    call assert(this%nlayer.GE.0,msg='DEVICE_check: nlayer is negative.',iostat=DEVICE_check)
    if(DEVICE_check.NE.0)return

    !check ffn hiddenlayer pointer array points to something and is correct size
    if(this%nlayer.GT.0)then
       call assert(associated(this%hiddenlayer)&
            ,msg='DEVICE_check: is not allocated hiddenlayer pointer array.',iostat=DEVICE_check)
       if(DEVICE_check.NE.0)return
       call assert(size(this%hiddenlayer).EQ.this%nlayer&
            ,msg='DEVICE_check: is not allocated hiddenlayer pointer array.',iostat=DEVICE_check)
       if(DEVICE_check.NE.0)return
    end if

    !check ffn targets in hiddenlayer pointer array
    if(this%nlayer.GT.0)then
       do i=1,this%nlayer
          call assert(check(this%hiddenlayer(i)%ffn).EQ.0&
               ,msg='DEVICE_check: bad ffn found in hiddenlayer pointer array.',iostat=DEVICE_check)
          if(DEVICE_check.NE.0)return
       end do
    end if

    !check rank
    call assert(this%rank,this%rank,msg='DEVICE_check: device rank has bad value.',iostat=DEVICE_check)
    if(DEVICE_check.NE.0)return
    call assert(this%rank.GE.0,msg='DEVICE_check: rank is negative.',iostat=DEVICE_check)
    if(DEVICE_check.NE.0)return

    !check hiddenlayer rank assignments
    if(this%nlayer.GT.0)then
       do i=1,this%nlayer
          call assert(this%hiddenlayer(i)%rank.EQ.this%hiddenlayer(i)%rank&
               ,msg='DEVICE_check: bad ffn rank value found in hiddenlayer array.',iostat=DEVICE_check)
          if(DEVICE_check.NE.0)return
          call assert(this%hiddenlayer(i)%rank.LE.this%rank&
               ,msg='DEVICE_check: hiddenlayer ffn rank is greater than device rank.',iostat=DEVICE_check)
          if(DEVICE_check.NE.0)return
       end do
    end if

    !check hiddenlayer derivate weight assignment
    if(this%nlayer.GT.0)then
       do i=1,this%nlayer
          call assert(associated(this%hiddenlayer(i)%gW)&
               ,msg='DEVICE_check: derivative weights not allocated.',iostat=DEVICE_check)
          if(DEVICE_check.NE.0)return
          call assert(size(this%hiddenlayer(i)%gW,1).EQ.size(this%hiddenlayer(i)%ffn%W,1)&
               ,msg='DEVICE_check: derivative weight dimension 1 does not equal ffn weight dim 1.',iostat=DEVICE_check)
          if(DEVICE_check.NE.0)return
          call assert(size(this%hiddenlayer(i)%gW,2).EQ.size(this%hiddenlayer(i)%ffn%W,2)&
               ,msg='DEVICE_check: derivative weight dimension 2 does not equal ffn weight dim 2.',iostat=DEVICE_check)
          if(DEVICE_check.NE.0)return
          call assert(all(this%hiddenlayer(i)%gW.EQ.this%hiddenlayer(i)%gW)&
               ,msg='DEVICE_check: hiddenlayer gW has bad values.',iostat=DEVICE_check)
          if(DEVICE_check.NE.0)return
       end do
    end if

    !check learningrate
    call assert(this%learningrate,this%learningrate,msg='DEVICE_check: learning rate has bad values.'&
         ,iostat=DEVICE_check)
    if(DEVICE_check.NE.0)return
    call assert(this%learningrate.GT.epsilon(1.0_double),msg='learningrate is tiny.',iostat=DEVICE_check)
    if(DEVICE_check.NE.0)return
          
  end function DEVICE_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the DEVICE methods.
  !> \return Will print to stdio result of tests
  !======================================================================
  subroutine DEVICE_test
    use testing_class
    use progressbar
    use MDutils

    type(DEVICE)::this
    type(ffn)::ffn1,ffn2

    integer(short)::ierr
    integer(long)::epoch,nepoch,ndump,batch,nbatch,sample,nsample,nx,resid
    real(double),allocatable::dataset(:,:),dataset2(:,:,:),histogram(:)
    real(double)::MSE0,MSE,MSEtemp,accuracy
    integer(long),allocatable::batchlen(:)
    character(len=1)::char1,char1out
    logical::bit1,bit2
    integer(long)::i,j,ncorrect

    write(*,*)'test DEVICE can be created with vislayer size specification.'
    call make(this,N=4)
    call assert(check(this).EQ.0,msg='DEVICE object was not created properly.')
    write(*,*)'test DEVICE kill method sets initiallization flag to false.'
    call kill(this)
    call assert(.not.this%initialized,msg='DEVICE object remains initialized after killed.')

    write(*,*)'test make method creates vislayer of correct size.'
    call make(this,N=4)
    call assert(this%vislayer%N.EQ.4,msg='device make method does not create vislayer of correct size')
    call kill(this)

    write(*,*)'test kill method kills vislayer'
    call assert(check(this%vislayer).NE.0,msg='device kill method does not kill vislayer')

    write(*,*)'test make method sets layer counter to zero'
    call make(this,N=4)
    call assert(this%nlayer.EQ.0,msg='make mehtod does not set layer counter to zero')
    call kill(this)

    write(*,*)'test make method does not allocate hiddenlayer.'
    call make(this,N=4)
    call assert(.not.associated(this%hiddenlayer),msg='device hidden layer array is allocated after make.')
    call kill(this)

    write(*,*)'test kill method deallocates hiddenlayer.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.associated(this%hiddenlayer),msg='device kill method does not deallocate hidden layer array.')

    write(*,*)'test link method allocates hiddenlayer.'
    call make(ffn1,N=2)
    call make(this,N=4)
    call link(this,ffn1)
    call assert(associated(this%hiddenlayer),msg='device link method does not allocate hidden layer array.')
    call kill(this)
    call kill(ffn1)

    write(*,*)'test device object and newffn pass check after link method.'
    call make(ffn1,N=2)
    call make(this,N=4)
    call link(this,ffn1)
    call assert(check(this).EQ.0,msg='device does not pass check after link method')
    call assert(check(ffn1).EQ.0,msg='new ffn not pass check after link method')
    call kill(this)
    call kill(ffn1)

    write(*,*)'test link method sets rank assigment is correct.'
    call make(ffn1,N=4)
    call make(this,N=3)
    call link(this,ffn1)
    call assert(this%rank.EQ.1,msg='link method does not set device rank properly.')
    call assert(this%hiddenlayer(1)%rank.EQ.1,msg='link method does not set hiddenlayer rank properly')
    call kill(this)
    call kill(ffn1)

    write(*,*)'test forward pass with linked single node ffn.'
    call make(this,N=1)
    call make(ffn1,N=1)
    call link(this,ffn1)
    ffn1%W=0._double !set all weights to zero
    ffn1%W(1,1)=1.0_double !set ffn weight to input node to 1
    this%vislayer%input=5.0_double
    call forwardpass(this)
    call assert(ffn1%layer%node(1),5.0_double,msg='forward pass does not propagate input signal properly.')
    call kill(ffn1)
    call kill(this)

    write(*,*)'test forward pass updated layer dnode when passed with derivative flag.'
    call make(this,N=1)
    call make(ffn1,N=1)
    ffn1%layer%dnode=0.3_double
    call link(this,ffn1)
    call forwardpass(this,derivative=.true.)
    call assert(ffn1%layer%dnode(1),1.0_double,msg='forward pass does set dnode properly with true derivative flag.')
    call kill(ffn1)
    call kill(this)

    write(*,*)'test forward pass with linked two layers with derivative flag.'
    call make(this,N=1)
    call make(ffn1,N=1)
    call make(ffn2,N=1)
    call link(this,ffn1)
    ffn1%W=0._double !set all weights to zero
    ffn1%W(1,1)=1.0_double !set ffn weight to input node to 1
    call link(this,ffn2)
    ffn2%W=0._double !set all weights to zero
    ffn2%W(1,1)=1.0_double !set ffn weight to input node to 1
    this%vislayer%input=5.0_double
    call forwardpass(this,derivative=.true.)
    call assert(ffn2%layer%node(1),5.0_double,msg='forward pass does not propagate input signal properly with derivative flag.')
    call kill(ffn1)
    call kill(this)

    write(*,*)'link method allocates derivative weights in hiddenlayer array'
    call make(this,N=3)
    call make(ffn1,N=2)
    call link(this,ffn1)
    call assert(associated(this%hiddenlayer(1)%gW),msg='link mehtod does not allocate derivative weights.')
    call kill(this)
    call kill(ffn1)

    write(*,*)'link method allocates derivative weights of correct size'
    call make(this,N=3)
    call make(ffn1,N=2)
    call link(this,ffn1)
    call assert(size(this%hiddenlayer(1)%gW,1).EQ.size(ffn1%W,1)&
         ,msg='link mehtod does not allocate derivative weights of correct size.')
    call assert(size(this%hiddenlayer(1)%gW,2).EQ.size(ffn1%W,2)&
         ,msg='link mehtod does not allocate derivative weights of correct size.')
    call kill(this)
    call kill(ffn1)


!!$    write(*,*)'test training by back propagation reduces MSE with zpot dataset on 2 layer device.'
!!$    call make(this,N=1)!for x value
!!$    this%learningrate=1E-8
!!$    nepoch=50000
!!$    call make(ffn1,N=30,activation='tanh')
!!$    call make(ffn2,N=1) !for y value
!!$    call link(this,ffn1)
!!$    call link(this,ffn2)
!!$    !allocate dataset
!!$    if(allocated(dataset))deallocate(dataset)
!!$    allocate(dataset(2,500)) 
!!$    !read training set
!!$    open(111,file='data/zpot.dat')
!!$    do i=1,500
!!$       read(111,*)dataset(:,i)
!!$    end do
!!$    close(111)
!!$    !get starting MSE
!!$    call measure(this,dataset,MSE=MSE0)
!!$    open(123,file='testbackprop_zpot.error')
!!$    do epoch=0,nepoch
!!$       call backprop(this,dataset)
!!$       if(mod(epoch,500).EQ.0)then
!!$          call assert(check(this).EQ.0,msg='backpropagation failed')
!!$          call progress(epoch,nepoch)
!!$          this%learningrate=this%learningrate*.95
!!$          !get new MSE
!!$          call measure(this,dataset,MSE=MSE,file='testbackprop_zpot.out')
!!$          write(123,*)epoch,MSE
!!$          flush(123)
!!$       end if
!!$    end do
!!$    close(123)
!!$    !get final MSE
!!$    call measure(this,dataset,MSE=MSE,file='testbackprop_zpot.out')
!!$    !cleanup data set 
!!$    if(allocated(dataset))deallocate(dataset)
!!$    call system('gnuplot -persist data/zpot.plt')
!!$    call assert(MSE.LT.MSE0,msg='MSE is not reduced by training with backpropagation.')
!!$    call kill(this)
!!$    call kill(ffn2)
!!$    call kill(ffn1)
!!$
!!$    write(*,*)'test training by back propagation reduces MSE with cubic dataset on 2 layer device.'
!!$    call make(this,N=1)
!!$    this%learningrate=1E-11
!!$    nepoch=50000
!!$    call make(ffn1,N=20,activation='tanh')
!!$    call make(ffn2,N=1)
!!$    call link(this,ffn1)
!!$    call link(this,ffn2)
!!$    ffn2%W(1,:)=-60.0
!!$    !allocate dataset
!!$    if(allocated(dataset))deallocate(dataset)
!!$    allocate(dataset(2,500)) 
!!$    !read training set
!!$    open(111,file='data/cubic.dat')
!!$    do i=1,500
!!$       read(111,*)dataset(:,i)
!!$    end do
!!$    close(111)
!!$    !get starting MSE
!!$    call measure(this,dataset,MSE=MSE0)
!!$    open(123,file='testbackprop_cubic.error')
!!$    do epoch=0,nepoch
!!$       call backprop(this,dataset)
!!$       if(mod(epoch,1000).EQ.0)then
!!$          call assert(check(this).EQ.0,msg='backpropagation failed')
!!$          call progress(epoch,nepoch)
!!$          this%learningrate=this%learningrate*.95
!!$          !get new MSE
!!$          call measure(this,dataset,MSE=MSE,file='testbackprop_cubic.out')
!!$          write(123,*)epoch,MSE
!!$          flush(123)
!!$       end if
!!$    end do
!!$    close(123)
!!$    !get final MSE
!!$    call measure(this,dataset,MSE=MSE,file='testbackprop_cubic.out')
!!$    !cleanup data set 
!!$    if(allocated(dataset))deallocate(dataset)
!!$    call system('gnuplot -persist data/cubic.plt')
!!$    call assert(MSE.LT.MSE0,msg='cubic dataset MSE is not reduced by training with backpropagation.')
!!$    call kill(this)
!!$    call kill(ffn2)
!!$    call kill(ffn1) 
!!$
!!$    write(*,*)'test backprop training reduces MSE for xor dataset with sigmoid hidden layer.'
!!$    call make(this,N=2)
!!$    this%learningrate=1E-5
!!$    nepoch=10000
!!$    call make(ffn1,N=100,activation='tanh')
!!$    call make(ffn2,N=1,activation='bernoulli')
!!$    call link(this,ffn1)
!!$    call link(this,ffn2)
!!$    !allocate dataset
!!$    if(allocated(dataset))deallocate(dataset)
!!$    allocate(dataset(3,500)) 
!!$    !read training set
!!$    open(111,file='data/xor.dat')
!!$    do i=1,500
!!$       read(111,*)dataset(:,i)
!!$    end do
!!$    close(111)
!!$    !get starting MSE
!!$    call measure(this,dataset,MSE=MSE0)
!!$    open(123,file='testbackprop_xor.error')
!!$    do epoch=0,nepoch
!!$       call backprop(this,dataset)
!!$       if(mod(epoch,nepoch/100).EQ.0)then
!!$          call assert(check(this).EQ.0,msg='backpropagation failed')
!!$          call progress(epoch,nepoch)
!!$          this%learningrate=this%learningrate*.95
!!$          !get new MSE
!!$          call measure(this,dataset,MSE=MSE,file='testbackprop_xor.out')
!!$          write(123,*)epoch,MSE
!!$          flush(123)
!!$       end if
!!$    end do
!!$    close(123)
!!$    !get final MSE
!!$    call measure(this,dataset,MSE=MSE,file='testbackprop_xor.out')
!!$    !cleanup data set 
!!$    if(allocated(dataset))deallocate(dataset)
!!$    call system('gnuplot -persist data/xor.plt')
!!$    call assert(MSE.LT.MSE0,msg='MSE is not reduced by training with backpropagation.')
!!$    call kill(this)
!!$    call kill(ffn2)
!!$    call kill(ffn1)
!!$
!!$     write(*,*)'test backprop training reduces MSE for xor dataset with sinusoid hidden layer.'
!!$    call make(this,N=2)
!!$    this%learningrate=1E-5
!!$    nepoch=10000
!!$    call make(ffn1,N=100,activation='oscillator')
!!$    call make(ffn2,N=1,activation='bernoulli')
!!$    call link(this,ffn1)
!!$    call link(this,ffn2)
!!$    !allocate dataset
!!$    if(allocated(dataset))deallocate(dataset)
!!$    allocate(dataset(3,500)) 
!!$    !read training set
!!$    open(111,file='data/xor.dat')
!!$    do i=1,500
!!$       read(111,*)dataset(:,i)
!!$    end do
!!$    close(111)
!!$    !get starting MSE
!!$    call measure(this,dataset,MSE=MSE0)
!!$    open(123,file='testbackprop_xorsinusoid.error')
!!$    do epoch=0,nepoch
!!$       call backprop(this,dataset)
!!$       if(mod(epoch,nepoch/100).EQ.0)then
!!$          call assert(check(this).EQ.0,msg='backpropagation failed')
!!$          call progress(epoch,nepoch)
!!$          this%learningrate=this%learningrate*.95
!!$          !get new MSE
!!$          call measure(this,dataset,MSE=MSE,file='testbackprop_xorsinusoid.out')
!!$          write(123,*)epoch,MSE
!!$          flush(123)
!!$       end if
!!$    end do
!!$    close(123)
!!$    !get final MSE
!!$    call measure(this,dataset,MSE=MSE,file='testbackprop_xorsinusoid.out')
!!$    !cleanup data set 
!!$    if(allocated(dataset))deallocate(dataset)
!!$    call system('gnuplot -persist data/xorsinusoid.plt')
!!$    call assert(MSE.LT.MSE0,msg='MSE is not reduced by training with backpropagation.')
!!$    call kill(this)
!!$    call kill(ffn2)
!!$    call kill(ffn1)
!!$

    write(*,*)'test backprop training reduces MSE for protstruct dataset with tanh hidden layer.'

    nx=5*9
    nbatch=111
    nsample=3
    nepoch=100
    ndump=100

    call make(this,N=nx)
    this%learningrate=1E-6
    call make(ffn1,N=1000,activation='oscillator')
    call make(ffn2,N=2,activation='bernoulli')
    call link(this,ffn1)
    call link(this,ffn2)

    !allocate batchlen
    if(allocated(batchlen))deallocate(batchlen)
    allocate(batchlen(nbatch)) 
    !allocate dataset2
    if(allocated(dataset2))deallocate(dataset2)
    allocate(dataset2(nx+2,498,nbatch)) 

    !record training set
    dataset2=0
    open(111,file='data/protstruct.seqdat.dat')
    read(111,*)!skip first seq break 
    do batch=1,nbatch
       ierr=0
       batchlen(batch)=0
       do while(ierr.EQ.0)
          batchlen(batch)=batchlen(batch)+1
          read(111,*,iostat=ierr)dataset2(:,batchlen(batch),batch)
       end do
    end do
    close(111)

    !get starting MSE
    MSE0=0.0_double
    do batch=1,nbatch
       call measure(this,dataset2(:,1:batchlen(batch),batch)&
            ,MSE=MSEtemp,nmeasure=nsample)
       !call measure(this,dataset2(:,1:batchlen(batch),batch)&
       !     ,BCE=MSEtemp,nmeasure=nsample)
       MSE0=MSE0+MSEtemp
    end do
    MSE0=MSE0/real(nbatch)

    !open output file
    open(123,file='testbackprop_protstruct.error')
    do epoch=0,nepoch
       do batch=1,nbatch
          call backprop(this,dataset2(:,1:batchlen(batch),batch))
       end do
       if(mod(epoch,nepoch/ndump).EQ.0.and.check(this).EQ.0)then
          call progress(epoch,nepoch)
          !this%learningrate=this%learningrate*.97
          !get new MSE
          MSE=0.0_double
          do batch=1,nbatch
             call measure(this,dataset2(:,1:batchlen(batch),batch)&
                  ,MSE=MSEtemp,nmeasure=nsample)
            ! call measure(this,dataset2(:,1:batchlen(batch),batch)&
            !      ,BCE=MSEtemp,nmeasure=nsample)
             MSE=MSE+MSEtemp
          end do
          MSE=MSE/real(nbatch)
          write(123,*)epoch,MSE
          flush(123)
       end if
    end do
    close(123)

    !write final prediction for batches 1-10
    nsample=10
    open(333,file='testbackprop_protstruct.out')
    write(333,*)'True State | Predicted States with Probability | Max Prob state | accuracy'
    if(allocated(histogram))deallocate(histogram)
    allocate(histogram(0:3))
    !batch=1
    accuracy=0
    do batch=1,10
       ncorrect=0
       do resid=1,batchlen(batch)
          !load input
          this%vislayer%input(1:nx)=dataset2(1:nx,resid,batch)
          
          !repeat measurement
          histogram=0.0_double
          do sample=1,nsample
             !forward pass input
             call forwardpass(this)
             
             !get output state and accumulate histogram
             i=int(this%hiddenlayer(this%nlayer)%ffn%layer%node(1))
             j=int(this%hiddenlayer(this%nlayer)%ffn%layer%node(2))
             histogram(i*2+j)=histogram(i*2+j)+1.0_double
             
          end do
          histogram=histogram/real(nsample)
          
          write(333,*)!start on newline
          
          !Translate input to human readable
          !write(333,fmt='(A1,1X)',advance='no') 'X'
          
          !Translate dataset label to human readable
          i=dataset2(nx+1,resid,batch)
          j=dataset2(nx+2,resid,batch)
          bit1=.false.
          bit2=.false.
          if(i.EQ.1)bit1=.true.
          if(j.EQ.1)bit2=.true.
          char1=secondarystruct(bit1,bit2)
          write(333,fmt='(2(A1,1X))',advance='no') char1,'|'
                    
          !write max prob predicted state
          i= (maxloc(histogram,1)-1)/2
          j= mod(maxloc(histogram,1)-1,2)
          bit1=.false.
          bit2=.false.
          if(i.EQ.1)bit1=.true.
          if(j.EQ.1)bit2=.true.
          char1out=secondarystruct(bit1,bit2)
          write(333,fmt='(2(A1,1X))',advance='no') char1out,'|'
          
          !write accuracy
          if(char1.EQ.char1out)ncorrect=ncorrect+1
          write(333,fmt='(F8.6,1X,A1,1X)',advance='no') ncorrect/real(resid),'|'

          !write predicted state with probabilities
          do i=0,1
             do j=0,1
                !Translate output state to human readable
                bit1=.false.
                bit2=.false.
                if(i.EQ.1)bit1=.true.
                if(j.EQ.1)bit2=.true.
                char1out=secondarystruct(bit1,bit2)
                write(333,fmt='(4(A1,1X,F8.6,1X))',advance='no')&
                     char1out,histogram(i*2+j)
             end do
          end do
          
       end do
       write(333,*)
       write(333,*)'batch=',batch, 'accuracy=',ncorrect/real(batchlen(batch))
       write(333,*)'---------------------------------------------------------'
       accuracy=accuracy+ncorrect/real(batchlen(batch))
    end do
    write(333,*)'Average accuracy=',accuracy/real(batch)
    close(333)
       
    call system('gnuplot -persist data/protstruct.plt')
    call assert(MSE.LT.MSE0,msg='MSE is not reduced by training with backpropagation.')
    
    !cleanup memory
    if(allocated(histogram))deallocate(histogram)
    if(allocated(batchlen))deallocate(batchlen)
    if(allocated(dataset2))deallocate(dataset2)
    call kill(this)
    call kill(ffn2)
    call kill(ffn1)
    

    write(*,*)'ALL DEVICE TESTS PASSED!'
  end subroutine DEVICE_test
  !-----------------------------------------

end module DEVICE_class

