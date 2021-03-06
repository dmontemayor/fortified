!>\brief
!! device class
!!\details
!! A linked list of ffn layers begining with a visible layer and followed by
!! a series of ffn objects. The device class can train the network by
!! backpropagation.
!<------------------------------------------------------------------------
module device_class
  use type_kinds
  use layer_class
  use ffn_class
  implicit none
  private

  public::device, device_test
  public::make, kill, status, backup, update, reset, check, describe
  public::link, forwardpass, backprop, measure

  type ffnptr
     type(ffn),pointer::ffn
     integer(long)::rank
     real(double),dimension(:,:),pointer::gW
  end type ffnptr

  type device
     logical::initialized=.false.
     character(len=label)::name='device'
     type(layer)::vislayer
     integer(long)::nlayer=huge(1)
     integer(long)::rank=huge(1)
     type(ffnptr),dimension(:),pointer::hiddenlayer
     real(double)::learningrate
  end type device

  !> calculate device perfomance
  interface measure
     module procedure device_measure
  end interface

  !> calculate error gradient by backpropagation
  interface backprop
     module procedure device_backprop
  end interface

  !> propagate input signal up to change top layerstates
  interface forwardpass
     module procedure device_forwardpass
  end interface

  !> link ffn to device
  interface link
     module procedure device_link
  end interface

  !> Creates the device object.
  interface make
     module procedure device_init
  end interface

  !> Destroys the device object.
  interface kill
     module procedure device_kill
  end interface

  !> Returns the current state of the device object.
  interface status
     module procedure device_status
  end interface

  !> Returns a plain text description of the device object.
  interface describe
     module procedure device_describe
  end interface
  
  !> Backups the current state of the device object.
  interface backup
     module procedure device_backup
  end interface

  !> Recaluclates the device object.
  interface update
     module procedure device_update
  end interface

  !> Reinitializes the device object.
  interface reset
     module procedure device_reset
  end interface

  !> Checks that the device object.
  interface check
     module procedure device_check
  end interface

contains
  !======================================================================
  !> \brief Retruns a description of device as a string.
  !> \param[in] this is the device object.
  !======================================================================
  character(len=comment) function device_describe(this)
    type(device),intent(in)::this
    character(len=5)::FMT='(A)'

    write(device_describe,FMT)'No description for device has been provided.'
   
  end function device_describe

  !======================================================================
  !> \brief calucate Device performance given given dataset
  !> \param this is the device object to be measured.
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
  !subroutine device_measure(this,dataset,mse,bce,file,nmeasure,xrange,xmin)
  subroutine device_measure(this,dataset,mse,bce,file,nmeasure)
    use testing_class
    type(device),intent(inout)::this
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
    call assert(size(dataset,1).EQ.nx+ny,msg='device_measure: dataset dim 1 is incongruent with device io',iostat=ierr)
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
    
  end subroutine device_measure

  !======================================================================
  !> \brief back propagate error graident given dataset.
  !> \param this is the device object to be measured.
  !> \param[in] dataset is a realmatrix of size  Ninput + Noutput x Nsample
  !>  containing a batch of Nsample training data vetor pairs.
  !=====================================================================
  subroutine device_backprop(this,dataset)
    !use testing_class
    type(device),intent(inout)::this
    real(double),intent(in)::dataset(:,:)

    integer(short)::ierr
    integer(long)::isample,nsample,nx,ny,irank,ilayer,inode,imark,isource
    real(double)::zsum

    !!check device
    !call assert(check(this).EQ.0,msg='device_backprop: input device did not pass check',iostat=ierr)
    !if(ierr.NE.0)return

    !setup io
    nx=this%vislayer%N
    do ilayer=1,this%nlayer
       if(this%hiddenlayer(ilayer)%ffn%nmark.EQ.0)&
            ny=this%hiddenlayer(ilayer)%ffn%layer%N
    end do

    !!check dataset
    !call assert(size(dataset,1).EQ.nx+ny,msg='device_backprop: dataset dim 1 is incongruent with device io',iostat=ierr)
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


  end subroutine device_backprop
  !======================================================================
  !> \brief propagate any input signal up the device object.
  !> \param this is the device object to be updated.
  !=====================================================================
  subroutine device_forwardpass(this,derivative)
    type(device),intent(inout)::this
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

  end subroutine device_forwardpass
 
 !======================================================================
  !> \brief links ffn to the device object.
  !> \param this is the device object to be initialized.
  !> \param ffn is an ffn to add to the device
  !=====================================================================
  subroutine device_link(this,newffn)
    type(device),intent(inout)::this
    type(ffn),intent(inout),target::newffn
    type(ffnptr),dimension(:),pointer::tmphiddenlayer

    integer(long)::i,maxrank=0

    if(this%nlayer.GT.0)then
       !create temphiddenlayer
       if(associated(tmphiddenlayer))nullify(tmphiddenlayer)
       allocate(tmphiddenlayer(this%nlayer))
       !copy hiddenlayer targets to temphiddenlayer and backup maxrank
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

  end subroutine device_link

  !======================================================================
  !> \brief Creates and initializes the device object.
  !> \param this is the device object to be initialized.
  !> \param[in] file is an optional string containing the name of a previously backupd device file.
  !=====================================================================
  subroutine device_init(this,file,N,activation,loss)
    type(device),intent(inout)::this
    character*(*),intent(in),optional::file
    integer(long),optional,intent(in)::N
    character*(*),optional,intent(in)::activation,loss

    integer(long)::M

    !softmax=
    !square =1/2N sum_i^N(y_i-yhat_i)
    !binarycrossentropy=-1/N sum_i^N(log(p_i)*yhat_i+log(1-p_i)*(1-yhat_i)

    !set default node size
    M=1
    if(present(N))M=N

    if(present(activation))then
       call make(this%vislayer,N=M,activation=activation)
    else
       call make(this%vislayer,N=M)
    end if

    this%nlayer=0
    this%rank=0
    this%learningrate=0.001

    !declare initialization complete
    this%initialized=.true.
  end subroutine device_init

  !======================================================================
  !> \brief Destroys the device object.
  !> \param this is the device object to be destroyed.
  !====================================================================
  subroutine device_kill(this)
    type(device),intent(inout)::this
 
    call kill(this%vislayer)
    if(associated(this%hiddenlayer))nullify(this%hiddenlayer)

    !un-initialized device object
    this%initialized=.false.
  end subroutine device_kill

  !======================================================================
  !> \brief Computes the current state of device object.
  !> \param this is the device  object to be updated.
  !======================================================================
  subroutine device_update(this)
    type(device),intent(inout)::this


  end subroutine device_update

  !======================================================================
  !> \brief Re-initiallizes the device object.
  !> \param THIS is the device  object to be re-initialized.
  !> \param STATE is an optional integer when 0 will deallocate all dynamic
  !>        memory and return the object in an un-initiallized state.
  !> \remark When STATE is not present dynamic memory will be re-allocated
  !======================================================================
  subroutine device_reset(this,state)
    type(device),intent(inout)::this
    integer(long),intent(in),optional::STATE
    logical::nullstate

    nullstate=.false.
    if(present(state).and.state.EQ.0)nullstate=.true.


    !memory management
    if(nullstate)then
       !nullify all pointers
       !******        Example - cleanup pointer attribute 'PPP'       ****
       !if(associated(this%PPP))nullify(this%PPP)
       !******************************************************************
       
       !kill all objects
       !**** example **********
       !call kill(this%object)
       !***********************
       
       !un-initialized device object
       this%initialized=.false.
    else
       !reallocate all dynamic memory
       !******        Example - cleanup pointer attribute 'PPP'       ****
       !******                  then reallocate memory                ****
       !if(associated(this%PPP))nullify(this%PPP)
       !allocate(this%PPP(0:this%npt-1))
       !******************************************************************
       
       !reset all objects
       !**** example **********
       !call reset(this%object)
       !***********************
    end if

    !Reset any attributes to satisfy re-initiation conditions
    !**********  Example - attribute 'var' is always initially a     ********
    !**********            Gaussian random number                    ******** 
    ! this%var=gran()
    !************************************************************************

  end subroutine device_reset

  !======================================================================
  !> \brief Backups the current state of the device object to file.
  !> \param[in] THIS is the device  object to be updated.
  !> \param[in] FILE is a string containing the location of the backup file.
  !======================================================================
  subroutine device_backup(this,file)
    use filemanager
    use string
    use testing_class
    type(device),intent(in)::this
    character*(*),intent(in)::file
    integer(short)::unit
    logical::fileisopen
    integer(long)::i,j
    
    !check input file
    inquire(file=file,opened=fileisopen,number=unit)
    if(unit.LT.0)unit=newunit()
    if(.not.fileisopen)open(unit,file=file)
    
    !check device object
    call assert(check(this).EQ.0,msg='device object does not pass check.')

    !always write the data type on the first line
    write(unit,*)'device'

    !******      Backup below all the derived type's attributes       ****
    !******         in the order the MAKE method reads them           ****


    !First, Scalar parameters
    !******          Example - Backup a scalar parameter            ******
    ! write(unit,*)this%var
    !*********************************************************************


    !Second, Dynamic arrays
    !***       Example - Backup an NxM matrix                          ***
    ! write(unit,*)((this%matrix(i,j),j=1,M),i=1,N)
    !*********************************************************************


    !Last,objects
    !******              Example - Backup an object            ***********
    ! call backup(this%object,file//'.object')
    ! write(unit,*)quote(file//'.object')!write the object location
    !*********************************************************************
    

    !finished writing all attributes - now close backup file
    close(unit)  
  end subroutine device_backup
  
  !======================================================================
  !> \brief Retrun the device object as a single line record entry.
  !> \param[in] this is the device object.
  !> \param[in] msg is an optional string message to annotate the status.
  !======================================================================
  character(len=line) function device_status(this,msg)
    type(device),intent(in)::this
    character*(*),intent(in),optional::msg
    character(len=7)::FMT='(A10)'

    write(device_status,FMT)'helloworld'

   
  end function device_status

  !======================================================================
  !> \brief Checks the device object.
  !> \param[in] this is the device object to be checked.
  !> \return Nothing if all checks pass or 1 and a warn for the first failed check.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function device_check(this)
    use testing_class
    type(device),intent(in)::this

    integer(long)::i
  
    !initiate with no problems found 
    device_check=0
    !call Note('Checking device.')

    !check that object is initialized
    call assert(this%initialized,msg='device_check: device object not initialized.',iostat=device_check)
    if(device_check.NE.0)return

    !check vislayer
    call assert(check(this%vislayer).EQ.0,msg='device_check: vislayer failed check.',iostat=device_check)
    if(device_check.NE.0)return

    !check nlayer
    call assert(this%nlayer,this%nlayer,msg='device_check: nlayer has bad values.',iostat=device_check)
    if(device_check.NE.0)return
    call assert(this%nlayer.GE.0,msg='device_check: nlayer is negative.',iostat=device_check)
    if(device_check.NE.0)return

    !check ffn hiddenlayer pointer array points to something and is correct size
    if(this%nlayer.GT.0)then
       call assert(associated(this%hiddenlayer)&
            ,msg='device_check: is not allocated hiddenlayer pointer array.',iostat=device_check)
       if(device_check.NE.0)return
       call assert(size(this%hiddenlayer).EQ.this%nlayer&
            ,msg='device_check: is not allocated hiddenlayer pointer array.',iostat=device_check)
       if(device_check.NE.0)return
    end if

    !check ffn targets in hiddenlayer pointer array
    if(this%nlayer.GT.0)then
       do i=1,this%nlayer
          call assert(check(this%hiddenlayer(i)%ffn).EQ.0&
               ,msg='device_check: bad ffn found in hiddenlayer pointer array.',iostat=device_check)
          if(device_check.NE.0)return
       end do
    end if

    !check rank
    call assert(this%rank,this%rank,msg='device_check: device rank has bad value.',iostat=device_check)
    if(device_check.NE.0)return
    call assert(this%rank.GE.0,msg='device_check: rank is negative.',iostat=device_check)
    if(device_check.NE.0)return

    !check hiddenlayer rank assignments
    if(this%nlayer.GT.0)then
       do i=1,this%nlayer
          call assert(this%hiddenlayer(i)%rank.EQ.this%hiddenlayer(i)%rank&
               ,msg='device_check: bad ffn rank value found in hiddenlayer array.',iostat=device_check)
          if(device_check.NE.0)return
          call assert(this%hiddenlayer(i)%rank.LE.this%rank&
               ,msg='device_check: hiddenlayer ffn rank is greater than device rank.',iostat=device_check)
          if(device_check.NE.0)return
       end do
    end if

    !check hiddenlayer derivate weight assignment
    if(this%nlayer.GT.0)then
       do i=1,this%nlayer
          call assert(associated(this%hiddenlayer(i)%gW)&
               ,msg='device_check: derivative weights not allocated.',iostat=device_check)
          if(device_check.NE.0)return
          call assert(size(this%hiddenlayer(i)%gW,1).EQ.size(this%hiddenlayer(i)%ffn%W,1)&
               ,msg='device_check: derivative weight dimension 1 does not equal ffn weight dim 1.',iostat=device_check)
          if(device_check.NE.0)return
          call assert(size(this%hiddenlayer(i)%gW,2).EQ.size(this%hiddenlayer(i)%ffn%W,2)&
               ,msg='device_check: derivative weight dimension 2 does not equal ffn weight dim 2.',iostat=device_check)
          if(device_check.NE.0)return
          call assert(all(this%hiddenlayer(i)%gW.EQ.this%hiddenlayer(i)%gW)&
               ,msg='device_check: hiddenlayer gW has bad values.',iostat=device_check)
          if(device_check.NE.0)return
       end do
    end if

    !check learningrate
    call assert(this%learningrate,this%learningrate,msg='device_check: learning rate has bad values.'&
         ,iostat=device_check)
    if(device_check.NE.0)return
    call assert(this%learningrate.GT.epsilon(1.0_double),msg='learningrate is tiny.',iostat=device_check)
    if(device_check.NE.0)return
          
  end function device_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the device methods.
  !> \return Will print to stdio result of tests
  !======================================================================
  subroutine device_test
    use testing_class
    use filemanager
    use progressbar
    use MDutils

    type(device)::this
    type(ffn)::ffn1,ffn2

    integer(short)::ierr
    integer(long)::epoch,nepoch,ndump,batch,nbatch,sample,nsample,nx,resid
    real(double),allocatable::dataset(:,:),dataset2(:,:,:),histogram(:)
    real(double)::MSE0,MSE,MSEtemp,accuracy
    integer(long),allocatable::batchlen(:)
    character(len=1)::char1,char1out
    logical::bit1,bit2
    integer(long)::i,j,ncorrect

    real(double)::dummy1,dummy2,dummy3,dummy4
    real(double),allocatable::datasetvar(:),datasetbar(:)
    integer(long)::pcaworksizemax,pcaworksize
    real(double),allocatable::PCAeigenvector(:,:),PCAeigenvalue(:),PCAwork(:)


    character(len=label)::string
    integer(long)::unit
    
    !verify device is compatible with current version
    include 'verification'

    write(*,*)'test device can be created with vislayer size specification.'
    call make(this,N=4)
    call assert(check(this).EQ.0,msg='device object was not created properly.')
    write(*,*)'test device kill method sets initiallization flag to false.'
    call kill(this)
    call assert(.not.this%initialized,msg='device object remains initialized after killed.')

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


    write(*,*)'test training by back propagation reduces MSE with zpot dataset with tanh regression.'
    call make(this,N=1)!for x value
    this%learningrate=1E-8
    nepoch=50000
    call make(ffn1,N=30,activation='tanh')
    call make(ffn2,N=1) !for y value
    call link(this,ffn1)
    call link(this,ffn2)
    !allocate dataset
    if(allocated(dataset))deallocate(dataset)
    allocate(dataset(2,500)) 
    !read training set
    open(111,file='data/zpot.dat')
    do i=1,500
       read(111,*)dataset(:,i)
    end do
    close(111)
    !get starting MSE
    call measure(this,dataset,MSE=MSE0)
    open(123,file='testbackprop_zpot.error')
    do epoch=0,nepoch
       call backprop(this,dataset)
       if(mod(epoch,500).EQ.0)then
          call assert(check(this).EQ.0,msg='backpropagation failed')
          call progress(epoch,nepoch)
          this%learningrate=this%learningrate*.95
          !get new MSE
          call measure(this,dataset,MSE=MSE,file='testbackprop_zpot.out')
          write(123,*)epoch,MSE
          flush(123)
       end if
    end do
    close(123)
    !get final MSE
    call measure(this,dataset,MSE=MSE,file='testbackprop_zpot.out')
    !cleanup data set 
    if(allocated(dataset))deallocate(dataset)
    call system('gnuplot -persist data/zpot.plt')
    call assert(MSE.LT.MSE0,msg='MSE is not reduced by training with backpropagation.')
    call kill(this)
    call kill(ffn2)
    call kill(ffn1)

    write(*,*)'test training by back propagation reduces MSE with cubic dataset with tanh regression.'
    call make(this,N=1)
    this%learningrate=1E-11
    nepoch=50000
    call make(ffn1,N=20,activation='tanh')
    call make(ffn2,N=1)
    call link(this,ffn1)
    call link(this,ffn2)
    ffn2%W(1,:)=-60.0
    !allocate dataset
    if(allocated(dataset))deallocate(dataset)
    allocate(dataset(2,500)) 
    !read training set
    open(111,file='data/cubic.dat')
    do i=1,500
       read(111,*)dataset(:,i)
    end do
    close(111)
    !get starting MSE
    call measure(this,dataset,MSE=MSE0)
    open(123,file='testbackprop_cubic.error')
    do epoch=0,nepoch
       call backprop(this,dataset)
       if(mod(epoch,1000).EQ.0)then
          call assert(check(this).EQ.0,msg='backpropagation failed')
          call progress(epoch,nepoch)
          this%learningrate=this%learningrate*.95
          !get new MSE
          call measure(this,dataset,MSE=MSE,file='testbackprop_cubic.out')
          write(123,*)epoch,MSE
          flush(123)
       end if
    end do
    close(123)
    !get final MSE
    call measure(this,dataset,MSE=MSE,file='testbackprop_cubic.out')
    !cleanup data set 
    if(allocated(dataset))deallocate(dataset)
    call system('gnuplot -persist data/cubic.plt')
    call assert(MSE.LT.MSE0,msg='cubic dataset MSE is not reduced by training with backpropagation.')
    call kill(this)
    call kill(ffn2)
    call kill(ffn1) 

    write(*,*)'test backprop training reduces MSE for xor dataset with tanh bernoulli regression.'
    call make(this,N=2)
    this%learningrate=1E-5
    nepoch=10000
    call make(ffn1,N=100,activation='tanh')
    call make(ffn2,N=1,activation='bernoulli')
    call link(this,ffn1)
    call link(this,ffn2)
    !allocate dataset
    if(allocated(dataset))deallocate(dataset)
    allocate(dataset(3,500)) 
    !read training set
    open(111,file='data/xor.dat')
    do i=1,500
       read(111,*)dataset(:,i)
    end do
    close(111)
    !get starting MSE
    call measure(this,dataset,MSE=MSE0)
    open(123,file='testbackprop_xor.error')
    do epoch=0,nepoch
       call backprop(this,dataset)
       if(mod(epoch,nepoch/100).EQ.0)then
          call assert(check(this).EQ.0,msg='backpropagation failed')
          call progress(epoch,nepoch)
          this%learningrate=this%learningrate*.95
          !get new MSE
          call measure(this,dataset,MSE=MSE,file='testbackprop_xor.out')
          write(123,*)epoch,MSE
          flush(123)
       end if
    end do
    close(123)
    !get final MSE
    call measure(this,dataset,MSE=MSE,file='testbackprop_xor.out')
    !cleanup data set 
    if(allocated(dataset))deallocate(dataset)
    call system('gnuplot -persist data/xor.plt')
    call assert(MSE.LT.MSE0,msg='MSE is not reduced by training with backpropagation.'&
         ,iostat=ierr)
    call kill(this)
    call kill(ffn2)
    call kill(ffn1)

     write(*,*)'test backprop training reduces MSE for xor dataset with oscillator bernoulli regression.'
    call make(this,N=2)
    this%learningrate=1E-5
    nepoch=10000
    call make(ffn1,N=100,activation='oscillator')
    call make(ffn2,N=1,activation='bernoulli')
    call link(this,ffn1)
    call link(this,ffn2)
    !allocate dataset
    if(allocated(dataset))deallocate(dataset)
    allocate(dataset(3,500)) 
    !read training set
    open(111,file='data/xor.dat')
    do i=1,500
       read(111,*)dataset(:,i)
    end do
    close(111)
    !get starting MSE
    call measure(this,dataset,MSE=MSE0)
    open(123,file='testbackprop_xorsinusoid.error')
    do epoch=0,nepoch
       call backprop(this,dataset)
       if(mod(epoch,nepoch/100).EQ.0)then
          call assert(check(this).EQ.0,msg='backpropagation failed')
          call progress(epoch,nepoch)
          this%learningrate=this%learningrate*.95
          !get new MSE
          call measure(this,dataset,MSE=MSE,file='testbackprop_xorsinusoid.out')
          write(123,*)epoch,MSE
          flush(123)
       end if
    end do
    close(123)
    !get final MSE
    call measure(this,dataset,MSE=MSE,file='testbackprop_xorsinusoid.out')
    !cleanup data set 
    if(allocated(dataset))deallocate(dataset)
    call system('gnuplot -persist data/xorsinusoid.plt')
    call assert(MSE.LT.MSE0,msg='MSE is not reduced by training with backpropagation.'&
         ,iostat=ierr)
    call kill(this)
    call kill(ffn2)
    call kill(ffn1)
    
    write(*,*)'test backprop training reduces MSE for protstruct dataset.'

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
    



    
    write(*,*)'ALL device TESTS PASSED!'
  end subroutine device_test
  !-----------------------------------------

end module device_class
