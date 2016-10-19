!>\brief
!! Class som
!!\details
!! Self organizing map som is a type of feed foward network distinguished
!! by a map of interconnected nodes in the map space.
!! weights to external input nodes. 
!<------------------------------------------------------------------------
module som_class
  use type_kinds
  use ffn_class
  implicit none
  private

  public::som, som_test
  public::make, kill, status, backup, update, reset, check, describe, trainstep

  type som
     logical::initialized=.false.
     character(len=label)::name='som'
     type(ffn)::ffn
     integer(long)::ND1,ND2,ND3
     real(double)::D1halfcell,D2halfcell,D3halfcell
     real(double),dimension(:,:),pointer::nodeCoord
     logical::mexhat,PBC
     real(double)::sigmastart,sigmadecay
     real(double)::learningstart,learningdecay
  end type som

  !> Train som object by one step.
  interface trainstep
     module procedure som_trainstep
  end interface

  !> Creates the som object.
  interface make
     module procedure som_init
  end interface

  !> Destroys the som object.
  interface kill
     module procedure som_kill
  end interface

  !> Returns the current state of the som object.
  interface status
     module procedure som_status
  end interface

  !> Returns a plain text description of the som object.
  interface describe
     module procedure som_describe
  end interface
  
  !> Backups the current state of the som object.
  interface backup
     module procedure som_backup
  end interface

  !> Recaluclates the som object.
  interface update
     module procedure som_update
  end interface

  !> Reinitializes the som object.
  interface reset
     module procedure som_reset
  end interface

  !> Checks that the som object.
  interface check
     module procedure som_check
  end interface

contains
  !======================================================================
  !> \brief Retruns a description of som as a string.
  !> \param[in] this is the som object.
  !======================================================================
  character(len=comment) function som_describe(this)
    type(som),intent(in)::this
    character(len=5)::FMT='(A)'

    write(som_describe,FMT)'No description for som has been provided.'
   
  end function som_describe

  !======================================================================
  !> \brief Evolves the som object by one training step.
  !> \param this is the som object to be evolved.
  !> \param[in] S is the current training step
  !=====================================================================
  subroutine som_trainstep(this,S)
    type(som),intent(inout)::this
    integer(long),intent(in)::S
    real(double)::dW(this%ffn%layer%N,this%ffn%nsource),NF(this%ffn%layer%N)

    integer(long)::i,j,winningnode,N,M
    real(double)::xmax,x,LF

    N=this%ffn%layer%N
    M=this%ffn%nsource

    !calculate distances
    do i=1,N
       do j=1,M
          dW(i,j)=this%ffn%source(j-1)%ptr-this%ffn%W(i,j-1)
       end do
    end do

    !find winning node
    winningnode=0
    do i=1,N
       x=sqrt(sum(dW(i,:)**2))
       if(i.EQ.1)then
          xmax=x
          winningnode=i
       end if
       if(x.LE.xmax)then
          xmax=x
          winningnode=i
       end if
    end do

    !evolve nodes parameterized by winning node
    LF=Lfunc(this,S)
    NF=NFunc(this,winningnode,S)
    do j=1,M
       this%ffn%W(:,j-1)=this%ffn%W(:,j-1)+LF*NF*dW(:,j)
    end do

  end subroutine som_trainstep

  !======================================================================
  !> \brief som learning function
  !> \param this is the som object to measure.
  !> \param[in] S is the current training step.
  !=====================================================================
  function Lfunc(this,S)
    use functions
    type(som),intent(inout)::this
    integer(long),intent(in)::S
    real(double)::Lfunc
    Lfunc=this%learningstart*exp(-this%learningdecay*S)
    return
  end function Lfunc

  !======================================================================
  !> \brief som neigborhood function
  !> \param this is the som object to measure.
  !> \param[in] A is the winning node.
  !> \param[in] S is the current training step.
  !=====================================================================
  function Nfunc(this,A,S)
    use functions
    type(som),intent(inout)::this
    integer(long),intent(in)::A,S
    real(double),dimension(this%ffn%layer%N)::Nfunc,x
    integer(long)::i
    real(double)::sigma
    do i=1,this%ffn%layer%N
       x(i)=cartdist(this,i,A)
    end do
    sigma=this%sigmastart*exp(-this%sigmadecay*S)
    !put in units of standard deviation sigma
    x=x/sigma
    if(this%mexhat)then
       Nfunc=mexhat(x)
    else
       Nfunc=gaussian(x)
    end if
  end function Nfunc

  !======================================================================
  !> \brief cartesian (L2norm) distance between two nodes in som.
  !> \param this is the som object to measure.
  !> \param[in] A is the starting node.
  !> \param[in] B is the ending node.
  !=====================================================================
  real(double) function cartdist(this,A,B)
    type(som),intent(inout)::this
    integer(long),intent(in)::A,B
    real(double)::x,y,z

    x=this%nodecoord(A,1)-this%nodecoord(B,1)
    y=this%nodecoord(A,2)-this%nodecoord(B,2)
    z=this%nodecoord(A,3)-this%nodecoord(B,3)
    if(this%PBC)then
       if(x.gt.this%D1halfcell)x=x-this%D1halfcell*2.0_double
       if(x.le.-this%D1halfcell)x=x+this%D1halfcell*2.0_double
       if(y.gt.this%D2halfcell)y=y-this%D2halfcell*2.0_double
       if(y.le.-this%D2halfcell)y=y+this%D2halfcell*2.0_double
       if(z.gt.this%D3halfcell)z=z-this%D3halfcell*2.0_double
       if(z.le.-this%D3halfcell)z=z+this%D3halfcell*2.0_double
    end if
    !cartdist=sqrt(sum((this%nodecoord(A,:)-this%nodecoord(B,:))**2))
    cartdist=sqrt(x*x+y*y+z*z)

    return
  end function cartdist

  !======================================================================
  !> \brief Creates and initializes the som object.
  !> \param this is the som object to be initialized.
  !> \param[in] file is an optional string containing the name of a previously backupd som file.
  !> \remark If no input file is provided the user must manually initialize THIS using stout.
  !=====================================================================
  subroutine som_init(this,N,ND1,ND2,ND3,hcp,file,activation)
    type(som),intent(inout)::this
    integer(long),intent(in),optional::N
    integer(long),intent(in),optional::ND1,ND2,ND3
    logical,intent(in),optional::hcp
    character*(*),intent(in),optional::file
    character*(*),optional,intent(in)::activation

    real(double)::root3,third,twothirdsroot6
    integer(long)::i,j,k,M

    !set default classifier size
    M=1
    if(present(N))M=N
    
    !set number of nodes on dimension 1
    this%ND1=M
    if(present(ND1))this%ND1=ND1

    !set number of nodes on dimension 2
    this%ND2=1
    if(present(ND2))this%ND2=ND2

    !set number of nodes on dimension 3
    this%ND3=1
    if(present(ND3))this%ND3=ND3

    if(present(activation))then
       call make(this%ffn,N=this%ND1*this%ND2*this%ND3,activation=activation)
    else
       call make(this%ffn,N=this%ND1*this%ND2*this%ND3)
    end if

    !set node coordinates
    root3=sqrt(3.0_double)
    third=1/3.0_double
    twothirdsroot6=2*sqrt(6.0_double)/3.0_double
    if(associated(this%nodecoord))nullify(this%nodeCoord)
    allocate(this%nodecoord(this%ffn%layer%N,3))
    M=0 !node counter
    do k=0,this%ND3-1
       do j=0,this%ND2-1
          do i=0,this%ND1-1
             M=M+1 !increment node counter
             if(present(hcp).and.hcp)then
                this%nodecoord(M,1)=2*i+mod(j+k,2)
                this%nodecoord(M,2)=root3*(j+third*mod(k,2))
                this%nodecoord(M,3)=twothirdsroot6*k
                this%nodecoord(M,:)=0.5_double*this%nodecoord(M,:)
             else
                this%nodecoord(M,1)=i
                this%nodecoord(M,2)=j
                this%nodecoord(M,3)=k
             end if
          end do
       end do
    end do
    !cell length
    if(present(hcp).and.hcp)then
       this%D1halfcell=this%ND1*.5_double
       this%D2halfcell=this%ND2*root3*.25_double
       this%D3halfcell=this%ND3*twothirdsroot6*.25_double
    else
       this%D1halfcell=this%ND1*.5_double
       this%D2halfcell=this%ND2*.5_double
       this%D3halfcell=this%ND3*.5_double
    end if

    !set default value for priodic boundary condition option
    this%PBC=.false.

    !set default value for mexhat neighborhood function option
    this%mexhat=.false.

    !set default value of sigmastart manual option
    this%sigmastart=1.0_double

    !set default value of sigmadecay manual option
    this%sigmadecay=0.0_double

    !set default value of learning factor manual option
    this%learningstart=0.01_double

    !set default value of learning decay manual option
    this%learningdecay=0.0_double

    !declare initialization complete
    this%initialized=.true.
  end subroutine som_init

  !======================================================================
  !> \brief Destroys the som object.
  !> \param this is the som object to be destroyed.
  !====================================================================
  subroutine som_kill(this)
    type(som),intent(inout)::this
 
    !kill the layer primitive
    call kill(this%ffn)

    !kill node coordinates
    if(associated(this%nodecoord))nullify(this%nodecoord)

    !un-initialize som object
    this%initialized=.false.

  end subroutine som_kill

  !======================================================================
  !> \brief Computes the current state of som object.
  !> \param this is the som  object to be updated.
  !======================================================================
  subroutine som_update(this)
    type(som),intent(inout)::this

  end subroutine som_update

  !======================================================================
  !> \brief Re-initiallizes the som object.
  !> \param THIS is the som  object to be re-initialized.
  !> \param STATE is an optional integer when 0 will deallocate all dynamic
  !>        memory and return the object in an un-initiallized state.
  !> \remark When STATE is not present dynamic memory will be re-allocated
  !======================================================================
  subroutine som_reset(this,state)
    type(som),intent(inout)::this
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
       
       !un-initialized som object
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

  end subroutine som_reset

  !======================================================================
  !> \brief Backups the current state of the som object to file.
  !> \param[in] THIS is the som  object to be updated.
  !> \param[in] FILE is a string containing the location of the backup file.
  !======================================================================
  subroutine som_backup(this,file)
    use filemanager
    use string
    use testing_class
    type(som),intent(in)::this
    character*(*),intent(in)::file
    integer(short)::unit
    logical::fileisopen
    integer(long)::i,j
    
    !check input file
    inquire(file=file,opened=fileisopen,number=unit)
    if(unit.LT.0)unit=newunit()
    if(.not.fileisopen)open(unit,file=file)
    
    !check som object
    call assert(check(this).EQ.0,msg='som object does not pass check.')

    !always write the data type on the first line
    write(unit,*)'som'

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
  end subroutine som_backup
  
  !======================================================================
  !> \brief Retrun the som object as a single line record entry.
  !> \param[in] this is the som object.
  !> \param[in] msg is an optional string message to annotate the status.
  !======================================================================
  character(len=line) function som_status(this,msg)
    type(som),intent(in)::this
    character*(*),intent(in),optional::msg
    character(len=7)::FMT='(A10)'

    write(som_status,FMT)'helloworld'

   
  end function som_status

  !======================================================================
  !> \brief Checks the som object.
  !> \param[in] this is the som object to be checked.
  !> \return Nothing if all checks pass or 1 and a warn for the first failed check.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function som_check(this)
    use testing_class
    type(som),intent(in)::this

    !initiate with no problems found 
    som_check=0

    !check that object is initialized
    call assert(this%initialized,msg='som_check: som object not initialized.',iostat=som_check)
    if(som_check.NE.0)return

    !check the ffn
    call assert(check(this%ffn).EQ.0,msg='som_check: failed ffn check!',iostat=som_check)
    if(som_check.NE.0)return

    !check map dimensions
    call assert(this%ND1*this%ND2*this%ND3.EQ.this%ffn%layer%N,msg='som_check: all nodes must fit in rectangular 3D grid.'&
         ,iostat=som_check)
    if(som_check.NE.0)return

    !check node coordinates
    call assert(associated(this%nodeCoord),msg='som_check: node coordinates not associated',iostat=som_check)
    if(som_check.NE.0)return
    call assert(size(this%nodecoord).EQ.this%ffn%layer%N*3,msg='som_check: size of node coordinates not equal to number of nodes.',&
         iostat=som_check)
    if(som_check.NE.0)return
    call assert(all(this%nodecoord.EQ.this%nodecoord),msg='som_check: node coordinates have NaN values.',iostat=som_check)
    if(som_check.NE.0)return

    !check sigmastart
    call assert(this%sigmastart,this%sigmastart,msg='som_check: sigmastart has NaN or huge vaule.',iostat=som_check)
    if(som_check.NE.0)return
    call assert(this%sigmastart.GT.epsilon(1.0_double),msg='som_check: sigmastart has tiny or negative value.',iostat=som_check)
    if(som_check.NE.0)return

    !check sigmadecay
    call assert(this%sigmadecay,this%sigmadecay,msg='som_check: sigmadecay has NaN or huge vaule.',iostat=som_check)
    if(som_check.NE.0)return
    call assert(this%sigmadecay.GE.0.0_double,msg='som_check: sigmadecay has negative value.',iostat=som_check)
    if(som_check.NE.0)return

    !check learningstart
    call assert(this%learningstart,this%learningstart,msg='som_check: learningstart has NaN or huge vaule.',iostat=som_check)
    if(som_check.NE.0)return
    call assert(this%learningstart.GE.0.0_double,msg='som_check: learningstart has negative value.',iostat=som_check)
    if(som_check.NE.0)return

    !check learningdecay
    call assert(this%learningdecay,this%learningdecay,msg='som_check: learningdecay has NaN or huge vaule.',iostat=som_check)
    if(som_check.NE.0)return
    call assert(this%learningdecay.GE.0.0_double,msg='som_check: learningdecay has negative value.',iostat=som_check)
    if(som_check.NE.0)return

  end function som_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the som methods.
  !> \param[in] this is the som object whose methods will be excercised.
  !> \return Nothing if all tests pass or 1 and a stop for the first failed test.
  !> \remark Will stop after first failed check.
  !======================================================================
  subroutine som_test
    use testing_class
    use filemanager
    use layer_class
    use progressbar
    use rand_class
    type(som)::this
    type(layer)::sourcelayer

    real(double),allocatable::vector(:),mat(:,:)
    real(double)::x,y,z
    integer(long)::i,j,k,kprime
    integer(long)::ierr,epoch,nres,count,nepoch
    character::residue,resname(20),classifier
    character(len=1000)::seq,labseq
    character(len=label)::string
    integer(long)::unit

    integer(long),parameter::nmaxsubject=100,nmaxCID=500,nmaxrecord=1000
    integer(long)::subjectlist(nmaxsubject),CIDlist(nmaxCID),nsubject,nCID
    integer(long),allocatable::nrecordpersubject(:)
    real(double)::dummy1,dummy2,dummy3,dummy4
    real(double),allocatable::dataset(:,:),datasetvar(:),datasetbar(:)
    real(double),allocatable::trainingbatch(:,:,:),trainingbatchlabel(:,:,:)
    integer(long)::nx,pcaworksizemax,pcaworksize

    real(double)::R0,RA,RB,RC,RD,norm,R0var,RDvar

    real(double),allocatable::PCAeigenvector(:,:),PCAeigenvalue(:),PCAwork(:)
    real(double)::RGBmin(3),RGBL(3)

    real(double)::fitness,overlap,Finit,Ffinal
    real(double),allocatable::activity(:)

    !verify som is compatible with current version
    include 'verification'

    write(*,*)'test som can be created with 4 nodes.'
    call make(this,N=4)
    call assert(check(this).EQ.0,msg='som object was not created properly.')
    call kill(this)

    write(*,*)'test som kill method sets initiallization flag to false.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.this%initialized,msg='som object remains initialized after killed.')

    write(*,*)'test som creates ffn primitive properly with 4 nodes.'
    call make(this,N=4)
    call assert(check(this%ffn).EQ.0,msg='som ffn primitive was not created properly')
    call kill(this)

    write(*,*)'test som kill method sets ffn primitive flag to false.'
    call make(this,N=4)
    call kill(this)
    call assert(.not.this%ffn%initialized,msg='som ffn primitive remains initialized after killed.')

    write(*,*)'test som default with all nodes on dimension 1'
    call make(this,N=4)
    call assert(this%ND1.EQ.4,msg='som did not default with all nodes on dimension.')
    call kill(this)

    write(*,*)'test number of nodes is overwritten to fit rectangular map space.'
    call make(this,N=5,ND1=2,ND2=2)
    call assert(this%ffn%layer%N.EQ.4,msg='som number of nodes does not fit rectangular map space.')
    call kill(this)

    write(*,*)'test node coordinates are nullified after kill'
    call make(this,N=9,ND1=3,ND2=3)
    call kill(this)
    call assert(.not.associated(this%nodeCoord),msg='som nodecoord remains associated after killed.')

    write(*,*)'test node coordinate for node 5 is (1,1,0) in 3x3 map'
    call make(this,N=9,ND1=3,ND2=3)
    call assert(this%nodecoord(5,1).EQ.1.0,msg='node 5 x coordinate is not 1.0 in 3x3 map')
    call assert(this%nodecoord(5,2).EQ.1.0,msg='node 5 y coordinate is not 1.0 in 3x3 map')
    call assert(this%nodecoord(5,3).EQ.0.0,msg='node 5 z coordinate is not 0.0 in 3x3 map')
    call kill(this)

    write(*,*)'test distance between 1st and 9th node in 3x3 map is sqrt(8).'
    call make(this,N=9,ND1=3,ND2=3)
    call assert(cartdist(this,1,9).EQ.sqrt(8.0_double),msg='distance between 1st and 9th node in 3x3 map is not sqrt(8).')
    call kill(this)

    write(*,*)'test distance between 8th and 21st node in 3x3x3 map is 3.'
    call make(this,N=27,ND1=3,ND2=3,ND3=3)
    call assert(cartdist(this,8,21).EQ.3.0_double,msg='distance between 8th and 21st node in 3x3x3 map is not 3.')
    call kill(this)

    write(*,*)'test som can be created with hexagonal close packed geometry'
    call make(this,N=27,ND1=3,ND2=3,ND3=3,hcp=.true.)
    call assert(check(this).EQ.0,msg='som cannot be created with hexagonal close packed geometry.')
    call kill(this)

    write(*,*)'test distance between 2nd and 3rd node in 2x2 hcp map is 1.'
    call make(this,N=4,ND1=2,ND2=2,hcp=.true.)
    call assert(cartdist(this,2,3),1.0_double,msg='distance between 2nd and 3rd node in 2x2 hcp map is not 1.')
    call kill(this)

    write(*,*)'test distance between 3rd and 8th node in 2x2x2 hcp map is 1.'
    call make(this,N=8,ND1=2,ND2=2,ND3=2,hcp=.true.)
    call assert(cartdist(this,3,8),1.0_double,msg='distance between 3rd and 8th node in 2x2x2 hcp map is not 1.')
    call kill(this)

    write(*,*)'test neigborhood function is 1 when node is same as winning node.'
    call make(this,N=10)
    if(allocated(vector))deallocate(vector)
    allocate(vector(10))
    vector=NFunc(this,5,0)
    call assert(vector(5).EQ.1_double,msg='neighborhood function is not 1 when node is same as winning node.')
    deallocate(vector)
    call kill(this)

    write(*,*)'test default neigborhood function is exp(-1/2) when node is nearest neighbor of winning node.'
    call make(this,N=10)
    if(allocated(vector))deallocate(vector)
    allocate(vector(10))
    vector=NFunc(this,5,0)
    call assert(vector(4).EQ.exp(-0.5_double),&
         msg='default neighborhood function is not exp(-1/2) when node is nearest neighbor of winning node.')
    deallocate(vector)
    call kill(this)

    write(*,*)'test som manual feature mexhat neighborhood function is 1 when node is winning node.'
    call make(this,N=10)
    this%mexhat=.true.
    if(allocated(vector))deallocate(vector)
    allocate(vector(10))
    vector=NFunc(this,5,0)
    call assert(vector(5).EQ.1_double,msg='mexhat neighborhood function is not 1 when node is same as winning node.')
    deallocate(vector)
    call kill(this)

    write(*,*)'test som manual feature mexhat neigborhood function is 0 when node is nearest neighbor of winning node.'
    call make(this,N=10)
    this%mexhat=.true.
    if(allocated(vector))deallocate(vector)
    allocate(vector(10))
    vector=NFunc(this,5,0)
    call assert(vector(4).EQ.0._double,&
         msg='mexhat neighborhood function is not 0 when node is nearest neighbor of winning node.')
    deallocate(vector)
    call kill(this)

    write(*,*)'test som manual feature sigmastart cannot be zero.'
    call make(this,N=10)
    this%sigmastart=0.0_double
    call assert(check(this).NE.0,msg='som manual feature should not be able to be set to zero.')
    call kill(this)

    write(*,*)'test som manual feature sigmastart set to 2 gives default Nfunc=exp(-1/8) for nearest neighbor.'
    call make(this,N=10)
    this%sigmastart=2.0_double
    if(allocated(vector))deallocate(vector)
    allocate(vector(10))
    vector=NFunc(this,5,0)
    call assert(vector(4).EQ.exp(-1/8.0_double)&
         ,msg='sigmastart set to 2 does not give default Nfunc=exp(-1/8) for nearest neighbor.')
    deallocate(vector)
    call kill(this)

    write(*,*)'test som manual feature sigmastart set to 2 gives mexhat Nfunc=exp(-1/8)*3/4 for nearest neighbor.'
    call make(this,N=10)
    this%sigmastart=2.0_double
    this%mexhat=.true.
    if(allocated(vector))deallocate(vector)
    allocate(vector(10))
    vector=NFunc(this,5,0)
    call assert(vector(4).EQ.exp(-1/8.0_double)*3/4.0_double&
         ,msg='sigmastart set to 2 does not give default Nfunc=exp(-1/8)*3/4 for nearest neighbor.')
    deallocate(vector)
    call kill(this)

    write(*,*)'test som manual feature sigmadecay cannot be negative.'
    call make(this,N=10)
    this%sigmadecay=-1.0_double
    call assert(check(this).NE.0,msg='som sigmadecay should not be able to be set to a negative value.')
    call kill(this)

    write(*,*)'test som manual feature sigmadecay set to 0.5 gives default Nfunc=exp(-0.5/exp(-1))&
         & for nearest neighbor after 1 step.'
    call make(this,N=10)
    this%sigmadecay=0.5_double
    if(allocated(vector))deallocate(vector)
    allocate(vector(10))
    vector=NFunc(this,5,1)
    call assert(vector(4),exp(-0.5_double/exp(-1.0_double))&
         ,msg='sigmadecay set to 0.5 does not give default Nfunc=exp(-0.5/exp(-1)) for nearest neighbor after 1 step.')
    deallocate(vector)
    call kill(this)

    write(*,*)'test som learning factor is 0.01 when created'
    call make(this,N=9,ND1=3,ND2=3)
    call assert(this%Learningstart.EQ.0.01_double,msg='som learningstart is not 0.01 when created.')
    call kill(this)

    write(*,*)'test som learning factor cannot be manually set to negative'
    call make(this,N=9,ND1=3,ND2=3)
    this%learningstart=-0.01
    call assert(check(this).NE.0,msg='som learningstart should not be able to be set to negative.')
    call kill(this)

    write(*,*)'test som learning decay is 0.0 when created'
    call make(this,N=9,ND1=3,ND2=3)
    call assert(this%Learningdecay.EQ.0.0_double,msg='som learningdecay is not 0.0 when created.')
    call kill(this)

    write(*,*)'test som learning decay cannot be manually set to negative'
    call make(this,N=9,ND1=3,ND2=3)
    this%learningdecay=-0.01
    call assert(check(this).NE.0,msg='som learningdecay should not be able to be set to negative.')
    call kill(this)

    write(*,*)'test learning function equals learningstart at zeroth learning step.'
    call make(this,N=10)
    this%learningstart=2.0_double
    call assert(LFunc(this,0).EQ.2.0,msg='learning function does not equal learningstart at zeroth learning step.')
    call kill(this)

    write(*,*)'test learning function equals exp(-1)/100 at first learning step when learningdecay is set to 1.'
    call make(this,N=10)
    this%learningdecay=1.0_double
    call assert(LFunc(this,1),.01_double*exp(-1.0_double)&
         ,msg='learning function does not equal exp(-1)/100 at first learning step when learningdecay is set to 1.')
    call kill(this)


    write(*,*)'test som manual feature PBC finds cardist of first and last node equal to 1.0 when true.'
    call make(this,N=4)
    this%PBC=.true.
    call assert(cartdist(this,1,4).EQ.1.0,msg='PBC=.true. does not return cartdist&
         & of 1.0 between first and last nodes for 1D SOM')
    call kill(this)

    write(*,*)'test som manual feature PBC finds cardist of first and last node equal to& 
         & sqrt(2.0) when true for 3X3 surf.'
    call make(this,N=9,ND1=3,ND2=3)
    this%PBC=.true.
    x=sqrt(2.0_double)
    call assert(cartdist(this,1,9).EQ.x,msg='PBC=.true. does not return cartdist&
         & of sqrt(2.0) between first and last nodes for 1D SOM')
    call kill(this)

    write(*,*)'test som training step evolves node toward source node at 1.0'
    call make (sourcelayer,N=1)
    sourcelayer%node=1.0_double
    call make(this,N=1)
    call link(this%ffn,sourcelayer)
    this%ffn%W=0.0_double
    x=0._double
    do i=0,this%ffn%nsource-1
       x=x+(this%ffn%W(1,i)-this%ffn%source(i)%ptr)**2
    end do
    x=sqrt(x)
    call trainstep(this,0)
    y=0_double
    do i=0,this%ffn%nsource-1
       y=y+(this%ffn%W(1,i)-this%ffn%source(i)%ptr)**2
    end do
    y=sqrt(y)
    call assert(y.LT.x,msg='som node did not evolve node toward source node at 1.0.')
    call kill(this)
    call kill(sourcelayer)


    write(*,*)'test som training step evolves node toward source node vector at (-1.0,-1.0)'
    call make (sourcelayer,N=2)
    sourcelayer%node(1:2)=-1.0_double
    call make(this,N=1)
    call link(this%ffn,sourcelayer)
    x=0.0_double
    do i=0,this%ffn%nsource-1
       x=x+(this%ffn%W(1,i)-this%ffn%source(i)%ptr)**2
    end do
    x=sqrt(x)
    call trainstep(this,0)
    y=0_double
    do i=0,this%ffn%nsource-1
       y=y+(this%ffn%W(1,i)-this%ffn%source(i)%ptr)**2
    end do
    y=sqrt(y)
    call assert(y.LT.x,msg='som node did not evolve node toward source node vector at (-1.0,-1.0).')
    call kill(this)
    call kill(sourcelayer)


    write(*,*)'test training som with doublewell data set evolves 2 classifiers to points (1,2) and (-1,-2).'
    call make (sourcelayer,N=2)
    call make(this,N=2)
    call link(this%ffn,sourcelayer)
    open(111,file='data/doublewell.dat')
    open(222,file='doublewelltraining.dat')
    this%sigmadecay=0.01
    this%learningdecay=0.001
    do epoch=0,500
       ierr=0
       do while(ierr.EQ.0)
          read(111,*,iostat=ierr)sourcelayer%input
          if(ierr.EQ.0)then
             call update(sourcelayer)
             call trainstep(this,epoch)
          end if
       end do
       rewind(111)
       write(222,*)epoch,this%ffn%W(1,1:2),this%ffn%W(2,1:2)
    end do
    close(111)
    close(222)
    x=sqrt((this%ffn%W(1,1)+1)**2+(this%ffn%W(1,2)+2)**2)
    y=sqrt((this%ffn%W(1,1)-1)**2+(this%ffn%W(1,2)-2)**2)
    if(x.LT.y)then !node 1 closer to point (-1,-2)
       y=sqrt((this%ffn%W(2,1)-1)**2+(this%ffn%W(2,2)-2)**2) !calc distance of node 2 to point (1,2)
    else !node 1 closer to point (1,2)
       x=sqrt((this%ffn%W(2,1)+1)**2+(this%ffn%W(2,2)+2)**2) !calc distance of node 2 to point (-1,-2)
    end if
    call assert(x.LT.0.16.and.y.LT.0.05,msg='som did not evolve well 2 nodes toward the centers of the doublepotential.')
    call kill(this)
    call kill(sourcelayer)


    write(*,*)'test training som with zpotential data set evolves 10 classifiers on 1D map toward z potential.'
    call make (sourcelayer,N=2)
    call make(this,N=10)
    call link(this%ffn,sourcelayer)
    open(111,file='data/zpot.dat')
    this%sigmadecay=0.01
    this%learningdecay=0.001
    do epoch=0,500
       ierr=0
       do while(ierr.EQ.0)
          read(111,*,iostat=ierr)sourcelayer%input
          if(ierr.EQ.0)then
             call update(sourcelayer)
             call trainstep(this,epoch)
          end if
       end do
       rewind(111)
    end do
    close(111)
    open(222,file='zpottraining.dat')
    do i=1,this%ffn%layer%N
       write(222,*)this%ffn%W(i,1:2)
    end do
    close(222)
    !calculate linear error with zpotential ridge equation y=x(x+1)(x-1)
    z=0.0_double
    do i=1,this%ffn%layer%N
       x=this%ffn%W(i,1)
       y=x*(x+1)*(x-1)
       z=z+(this%ffn%W(i,2)-y)**2
    end do
    z=z/500.0
    call assert(z.LT.1.5E-3,msg='som zpotential training MSE is greater than 1.5E-3')
    call kill(this)
    call kill(sourcelayer)

!!$stop
!!$    write(*,*)'test training som with protstruct dataset to recover 3 classifiers&
!!$         & alphahelix, betasheet, and randomcoil with acceptable accuracy.'
!!$    resname(1)='A'!Alanine ALA
!!$    resname(2)='R'!Arginine ARG
!!$    resname(3)='N'!Asparagine ASN
!!$    resname(4)='D'!Aspartic acid ASP
!!$    resname(5)='C'!Cysteine CYS
!!$    resname(6)='Q'!Glutamine GLN
!!$    resname(7)='E'!Glutamic acid GLU
!!$    resname(8)='G'!Glycine GLY
!!$    resname(9)='H'!Histidine HIS
!!$    resname(10)='I'!Isoleucine ILE
!!$    resname(11)='L'!Leucine LEU
!!$    resname(12)='K'!Lysine LYS
!!$    resname(13)='M'!Methionine MET
!!$    resname(14)='F'!Phenylalanine PHE
!!$    resname(15)='P'!Proline PRO
!!$    resname(16)='S'!Serine SER
!!$    resname(17)='T'!Threonine THR
!!$    resname(18)='W'!Tryptophan TRP
!!$    resname(19)='Y'!Tyrosine TYR
!!$    resname(20)='V'!Valine VAL
!!$    call make (sourcelayer,N=60)
!!$    call make(this,N=3)
!!$    call link(this%ffn,sourcelayer)
!!$    this%sigmadecay=0.001
!!$    !this%learningdecay=0.001
!!$    open(111,file='data/protstructdataset.dat')
!!$    do epoch=0,5!000
!!$       write(*,*)'epoch=',epoch
!!$       ierr=0
!!$       do while(ierr.EQ.0)
!!$          read(111,*,iostat=ierr)residue
!!$          if(ierr.EQ.0)then
!!$             select case (residue)
!!$             case('#')!comment
!!$             case('<')!begin of sequence
!!$                !clear sequence register
!!$                nres=0
!!$                seq=''
!!$             case('e')!end of sequence 
!!$                !begin sequence analysis
!!$                !write(*,*)nres,trim(seq)
!!$                !loop over residues
!!$                do i=1,nres
!!$                   sourcelayer%input=0._double
!!$                   !get current residue
!!$                   do j=1,20
!!$                      if(seq(i:i).EQ.resname(j))sourcelayer%input(j)=1.0_double
!!$                   end do
!!$                   !get next residude distance
!!$                   do j=1,20
!!$                      k=nres
!!$                      do while (k.GT.i)
!!$                         if(seq(k:k).EQ.resname(j))sourcelayer%input(20+j)=k-i
!!$                         k=k-1
!!$                      end do
!!$                   end do
!!$                   !get previous residue distance
!!$                   do j=1,20
!!$                      k=1
!!$                      do while (k.LT.i)
!!$                         if(seq(k:k).EQ.resname(j))sourcelayer%input(40+j)=i-k
!!$                         k=k+1
!!$                      end do
!!$                   end do
!!$
!!$                   !load points into som and train
!!$                   write(*,*)trim(seq)
!!$                   write(*,*)sourcelayer%input(1:20)
!!$                   write(*,*)sourcelayer%input(21:39)
!!$                   write(*,*)sourcelayer%input(41:60)
!!$                   stop
!!$                   call update(sourcelayer)
!!$                   call trainstep(this,epoch)
!!$                end do
!!$             case default
!!$                if(any(resname.EQ.residue))then
!!$                   nres=nres+1
!!$                   !push to sequence
!!$                   seq=trim(seq)//residue
!!$                end if
!!$             end select
!!$          end if
!!$       end do
!!$       rewind(111)
!!$    end do
!!$    close(111)
!!$
!!$    !calculate class correlation matrix
!!$    if(allocated(mat))deallocate(mat)
!!$    allocate(mat(3,3))
!!$    mat=0.0_double
!!$    count=0
!!$    nres=0
!!$    open(111,file='data/protstructdataset.dat')
!!$    ierr=0
!!$    do while(ierr.EQ.0)
!!$       read(111,*,iostat=ierr)residue,classifier
!!$       if(ierr.EQ.0)then
!!$          select case (residue)
!!$          case('#')!comment
!!$          case('<')!begin of sequence
!!$             !clear sequence register
!!$             nres=0
!!$             seq=''
!!$             labseq=''
!!$          case('e')!end of sequence 
!!$             !begin sequence analysis
!!$             !write(*,*)nres,trim(seq)
!!$             !write(*,*)nres,trim(labseq)
!!$             !write(*,*)seq(nres:nres)//' '//labseq(nres:nres)
!!$             !stop
!!$             !loop over residues
!!$             do i=1,nres
!!$                sourcelayer%input=0._double
!!$                !get current residue
!!$                do j=1,20
!!$                   if(seq(i:i).EQ.resname(j))sourcelayer%input(j)=1.0_double
!!$                end do
!!$                !get next residude distance
!!$                do j=1,20
!!$                   k=nres
!!$                   do while (k.GT.i)
!!$                      if(seq(k:k).EQ.resname(j))sourcelayer%input(20+j)=k-i
!!$                      k=k-1
!!$                   end do
!!$                end do
!!$                !get previous residue distance
!!$                do j=1,20
!!$                   k=1
!!$                   do while (k.LT.i)
!!$                      if(seq(k:k).EQ.resname(j))sourcelayer%input(40+j)=i-k
!!$                      k=k+1
!!$                   end do
!!$                end do
!!$                !find nearest som classifier
!!$                do j=1,3
!!$                   !get distance of input coord with som node j
!!$                   x=sqrt(sum( (this%ffn%W(j,1:60)-sourcelayer%input(1:60))**2))
!!$                   if(j.EQ.1)then
!!$                      !autoamtic save mindistance z and classifier k if first classifier
!!$                      z=x
!!$                      k=j
!!$                   end if
!!$                   if(x.LE.z)then
!!$                      !save mindistance z and classifier k
!!$                      z=x
!!$                      k=j
!!$                   end if
!!$                end do
!!$                
!!$                !correlate som classifier k with training classifier j
!!$                select case(labseq(i:i))
!!$                case ('h')
!!$                   j=1
!!$                case ('e')
!!$                   j=2
!!$                case ('_')
!!$                   j=3
!!$                case default
!!$                   write(*,*)nres,i
!!$                   write(*,*)trim(seq)
!!$                   write(*,*)trim(labseq)
!!$                   write(*,*)'unkown classifier program will stop. ***'//labseq(i:i)//'***'
!!$                   stop
!!$                end select
!!$                mat(k,j)=mat(k,j)+1_double
!!$                count=count+1
!!$             end do
!!$          case default
!!$             if(any(resname.EQ.residue))then
!!$                nres=nres+1
!!$                !push to sequence
!!$                seq=trim(seq)//residue
!!$                labseq=trim(labseq)//classifier
!!$             end if
!!$          end select
!!$       end if
!!$    end do
!!$    close(111)
!!$    mat=mat/real(count,double)
!!$    write(*,*)'-----------correlation matrix-----------'
!!$    write(*,*)'              h                         e                          _'
!!$    do k=1,3
!!$       write(*,*)k,(mat(k,j),j=1,3)
!!$    end do
!!$    write(*,*)'----------------------------------------'
!!$
!!$    open(222,file='protstructtraining.dat')
!!$    write(222,*)'-----------correlation matrix-----------'
!!$    write(222,*)'              h                         e                          _'
!!$    do k=1,3
!!$       write(222,*)k,(mat(k,j),j=1,3)
!!$    end do
!!$    write(222,*)'----------------------------------------'
!!$    close(222)
!!$!    call assert(z.LT.1.07E-3,msg='som zpotential training MSE is greater than 1.07E-3')
!!$    call kill(this)
!!$    call kill(sourcelayer)

!!!!!!!!!!!!!!!!!

    write(*,*)'choose training-validation and unknown data'
    !perception data column labels 
    !subjects #,CID,dilution,detection,intensity,pleasantness,familiarity,identification,edible,bakery,sweet,fruit,fish,garlic,spices,cold,sour,burnt,acid,warm,musky,sweaty,ammonia/urinous,decayed,wood,grass,flower,chemical
    !allocate dataset
    nx=24
    if(allocated(dataset))deallocate(dataset)
    allocate(dataset(nx,35090))

    !load available data
    dataset=0
    subjectlist=0
    CIDlist=0
    nsubject=0
    nCID=0
    open(111,file='data/allsubjectperception.dat')
    ierr=0
    i=0
    do while(ierr.EQ.0)
       i=i+1
       read(111,*,iostat=ierr)dummy1,dummy2,dummy3,dummy4,dataset(:,i)

       !count unique CIDs and number of subjects
       if(all(subjectlist.NE.int(dummy1)))then
          nsubject=nsubject+1
          subjectlist(nsubject)=int(dummy1)
       end if
       if(all(CIDlist.NE.int(dummy2)))then
          nCID=nCID+1
          CIDlist(nCID)=int(dummy2)
       end if

    end do
    close(111)

    write(*,*)nsubject,'subjects found'!,subjectlist(1:nsubject)
    write(*,*)nCID,'CIDs found'!,CIDlist(1:nCID)

    !randomly pick 70% of subjects 
    do k=1,1000 !swap values in subjectlist 1000 times
       i=ceiling(uran()*nsubject) !choose random array position to swap
       j=ceiling(uran()*nsubject) !choose random array position to swap
       !save value in position i
       ierr=subjectlist(i)
       !swap j value to position i
       subjectlist(i)=subjectlist(j)
       !return saved value to position j
       subjectlist(j)=ierr
    end do
    open(123,file='knownsubjects.dat')
    do i=1,int(nsubject*.7)
       write(123,*)subjectlist(i)
    end do
    close(123)

    !randomly pick 70% of CIDs 
    do k=1,1000 !swap values in CIDlist 1000 times
       i=ceiling(uran()*nCID) !choose random array position to swap
       j=ceiling(uran()*nCID) !choose random array position to swap
       !save value in position i
       ierr=CIDlist(i)
       !swap j value to position i
       CIDlist(i)=CIDlist(j)
       !return saved value to position j
       CIDlist(j)=ierr
    end do
    open(123,file='knownCIDs.dat')
    do i=1,int(nCID*.7)
       write(123,*)CIDlist(i)
    end do
    close(123)

    !clean up dataset memory
    if(allocated(dataset))deallocate(dataset)

!!!!!
    write(*,*)'test training som with frist 3PCA perception data&
         & set evolves 49 classifiers on 2D hcp map toward perception data.'
    call make (sourcelayer,N=3)
    call make(this,N=49,ND1=7,ND2=7,hcp=.true.)
    call link(this%ffn,sourcelayer)


    !get list of known subjects
    open(123,file='knownsubjects.dat')
    ierr=0
    nsubject=0
    subjectlist=0
    do while (ierr.EQ.0)
       nsubject=nsubject+1
       read(123,*,iostat=ierr) subjectlist(nsubject)
    end do
    close(123)
    !get list of known CIDs
    open(123,file='knownCIDs.dat')
    ierr=0
    nCID=0
    CIDlist=0
    do while (ierr.EQ.0)
       nCID=nCID+1
       read(123,*,iostat=ierr) CIDlist(nCID)
    end do
    close(123)
    !allocate number of records per known subjects array
    if(allocated(nrecordpersubject))deallocate(nrecordpersubject)
    allocate(nrecordpersubject(nsubject))
    nrecordpersubject=0
    !allocate dataset mean and variance 
    if(allocated(datasetbar))deallocate(datasetbar)
    allocate(datasetbar(nx))
    datasetbar=0.0_double
    if(allocated(datasetvar))deallocate(datasetvar)
    allocate(datasetvar(nx))
    datasetvar=0.0_double
    !allocate training batches
    nx=24
    if(allocated(trainingbatch))deallocate(trainingbatch)
    allocate(trainingbatch(nx,nmaxrecord,nsubject))
    !get training data as batches of known subjects
    open(111,file='data/allsubjectperception.dat')
    ierr=0
    nrecordpersubject=0
    do while(ierr.EQ.0)
       read(111,*,iostat=ierr)dummy1,dummy2,dummy3,dummy4,datasetbar
       do i=1,nsubject
          if (subjectlist(i).EQ.dummy1.and.any(CIDlist.EQ.int(dummy2)))then
             nrecordpersubject(i)=nrecordpersubject(i)+1
             call assert(nrecordpersubject(i).LE.nmaxrecord,msg='nrecord per subject is greater than nmaxrecord')
             trainingbatch(:,nrecordpersubject(i),i)=datasetbar
          end if
       end do
    end do
    close(111)

    !calculate dataset mean and variance per column
    datasetbar=0._double
    datasetbar=0._double
    do i=1,nsubject
       do j=1,nrecordpersubject(i)
          datasetbar=datasetbar+trainingbatch(:,j,i)
          datasetvar=datasetvar+trainingbatch(:,j,i)**2
       end do
    end do
    norm=real(sum(nrecordpersubject),double)
    datasetbar=datasetbar/norm
    datasetvar=sqrt(datasetvar/norm-datasetbar**2)

    write(*,*)"Perception label         mean          and          stdev"    
    write(*,*)"intensity      ",datasetbar(1),datasetvar(1)  
    write(*,*)"pleasantness   ",datasetbar(2),datasetvar(2)  
    write(*,*)"familiarity    ",datasetbar(3),datasetvar(3)  
    write(*,*)"identification ",datasetbar(4),datasetvar(4)  
    write(*,*)"edible         ",datasetbar(5),datasetvar(5)  
    write(*,*)"bakery         ",datasetbar(6),datasetvar(6)  
    write(*,*)"sweet          ",datasetbar(7),datasetvar(7)  
    write(*,*)"fruit          ",datasetbar(8),datasetvar(8)  
    write(*,*)"fish           ",datasetbar(9),datasetvar(9)  
    write(*,*)"garlic         ",datasetbar(10),datasetvar(10)
    write(*,*)"spices         ",datasetbar(11),datasetvar(11)
    write(*,*)"cold           ",datasetbar(12),datasetvar(12)
    write(*,*)"sour           ",datasetbar(13),datasetvar(13)
    write(*,*)"burnt          ",datasetbar(14),datasetvar(14)
    write(*,*)"acid           ",datasetbar(15),datasetvar(15)
    write(*,*)"warm           ",datasetbar(16),datasetvar(16)
    write(*,*)"musky          ",datasetbar(17),datasetvar(17)
    write(*,*)"sweaty         ",datasetbar(18),datasetvar(18)
    write(*,*)"ammonia/urinous",datasetbar(19),datasetvar(19)
    write(*,*)"decayed        ",datasetbar(20),datasetvar(20)
    write(*,*)"wood           ",datasetbar(21),datasetvar(21)
    write(*,*)"grass          ",datasetbar(22),datasetvar(22)
    write(*,*)"flower         ",datasetbar(23),datasetvar(23)
    write(*,*)"chemical       ",datasetbar(24),datasetvar(24)
    write(*,*)

    !allocate PCA eigenvecotor and eignevalue memory
    if(allocated(PCAeigenvector))deallocate(PCAeigenvector)
    allocate(PCAeigenvector(nx,nx))
    if(allocated(PCAeigenvalue))deallocate(PCAeigenvalue)
    allocate(PCAeigenvalue(nx))
    if(allocated(PCAwork))deallocate(PCAwork)
    pcaworksizemax=3*nx
    allocate(PCAwork(pcaworksizemax))

    !calculate correlation matrix
    !write(*,*)'calculating perception correlation matrix'
    PCAeigenvector=0._double
    do k=1,nsubject
       do i=1,nx
          do j=i,nx
             PCAeigenvector(i,j)=PCAeigenvector(i,j)&
                  +sum((trainingbatch(i,1:nrecordpersubject(k),k)-datasetbar(i))&
                  *(trainingbatch(j,1:nrecordpersubject(k),k)-datasetbar(j)))&
                  /(datasetvar(i)*datasetvar(j))
             !PCAeigenvector(j,i)=PCAeigenvector(i,j)
          end do
       end do
    end do
    norm=real(sum(nrecordpersubject),double)
    PCAeigenvector=PCAeigenvector/norm

    !query the optimal workspace
    pcaworksize=-1
    call dsyev('V','U',nx,PCAeigenvector,nx,PCAeigenvalue,pcawork,pcaworksize,i)
    pcaworksize = min(pcaworksizemax,int(pcawork(1)))
    !compute eigenvalues and eigenvectors of a real symmetric matrix
    call dsyev('V','U',nx,PCAeigenvector,nx,PCAeigenvalue,pcawork,pcaworksize,i)
    !check for convergence
    call assert(i.EQ.0,msg='lapack routine dsyev failed to compute eigenvalues')

    write(*,*)'PCA eigen values            percentage,                cumulative percentange,  1- cum percentage'
    do i=1,nx
       write(*,*)pcaeigenvalue(i),pcaeigenvalue(i)/sum(pcaeigenvalue)&
            ,sum(pcaeigenvalue(1:i))/sum(pcaeigenvalue),1.0-sum(pcaeigenvalue(1:i))/sum(pcaeigenvalue)
    end do
    write(*,*)

    write(*,*)"Perception label  PCA1 vector,               PCA2 vector,              PCA3 vector"    
    write(*,*)"intensity      ",pcaeigenvector(nx,1),pcaeigenvector(nx-1,1),pcaeigenvector(nx-2,1)
    write(*,*)"pleasantness   ",pcaeigenvector(nx,2),pcaeigenvector(nx-1,2),pcaeigenvector(nx-2,2)
    write(*,*)"familiarity    ",pcaeigenvector(nx,3),pcaeigenvector(nx-1,3),pcaeigenvector(nx-2,3)
    write(*,*)"identification ",pcaeigenvector(nx,4),pcaeigenvector(nx-1,4),pcaeigenvector(nx-2,4)
    write(*,*)"edible         ",pcaeigenvector(nx,5),pcaeigenvector(nx-1,5),pcaeigenvector(nx-2,5)
    write(*,*)"bakery         ",pcaeigenvector(nx,6),pcaeigenvector(nx-1,6),pcaeigenvector(nx-2,6)
    write(*,*)"sweet          ",pcaeigenvector(nx,7),pcaeigenvector(nx-1,7),pcaeigenvector(nx-2,7)
    write(*,*)"fruit          ",pcaeigenvector(nx,8),pcaeigenvector(nx-1,8),pcaeigenvector(nx-2,8)
    write(*,*)"fish           ",pcaeigenvector(nx,9),pcaeigenvector(nx-1,9),pcaeigenvector(nx-2,9)
    write(*,*)"garlic         ",pcaeigenvector(nx,10),pcaeigenvector(nx-1,10),pcaeigenvector(nx-2,10)
    write(*,*)"spices         ",pcaeigenvector(nx,11),pcaeigenvector(nx-1,11),pcaeigenvector(nx-2,11)
    write(*,*)"cold           ",pcaeigenvector(nx,12),pcaeigenvector(nx-1,12),pcaeigenvector(nx-2,12)
    write(*,*)"sour           ",pcaeigenvector(nx,13),pcaeigenvector(nx-1,13),pcaeigenvector(nx-2,13)
    write(*,*)"burnt          ",pcaeigenvector(nx,14),pcaeigenvector(nx-1,14),pcaeigenvector(nx-2,14)
    write(*,*)"acid           ",pcaeigenvector(nx,15),pcaeigenvector(nx-1,15),pcaeigenvector(nx-2,15)
    write(*,*)"warm           ",pcaeigenvector(nx,16),pcaeigenvector(nx-1,16),pcaeigenvector(nx-2,16)
    write(*,*)"musky          ",pcaeigenvector(nx,17),pcaeigenvector(nx-1,17),pcaeigenvector(nx-2,17)
    write(*,*)"sweaty         ",pcaeigenvector(nx,18),pcaeigenvector(nx-1,18),pcaeigenvector(nx-2,18)
    write(*,*)"ammonia/urinous",pcaeigenvector(nx,19),pcaeigenvector(nx-1,19),pcaeigenvector(nx-2,19)
    write(*,*)"decayed        ",pcaeigenvector(nx,20),pcaeigenvector(nx-1,20),pcaeigenvector(nx-2,20)
    write(*,*)"wood           ",pcaeigenvector(nx,21),pcaeigenvector(nx-1,21),pcaeigenvector(nx-2,21)
    write(*,*)"grass          ",pcaeigenvector(nx,22),pcaeigenvector(nx-1,22),pcaeigenvector(nx-2,22)
    write(*,*)"flower         ",pcaeigenvector(nx,23),pcaeigenvector(nx-1,23),pcaeigenvector(nx-2,23)
    write(*,*)"chemical       ",pcaeigenvector(nx,24),pcaeigenvector(nx-1,24),pcaeigenvector(nx-2,24)
    write(*,*)

    !rotate perception data
    do k=1,nsubject
       do j=1,nrecordpersubject(k)
          trainingbatch(:,j,k)=matmul(transpose(pcaeigenvector),trainingbatch(:,j,k))
       end do
    end do

    !dump som evolution in rgb coordinates
    open(222,file='SOMperceptionRGB.dat')
    !get rgb domain of ffn weights
    RGBmin(1)=minval(this%ffn%W(:,1))
    RGBL(1)=maxval(this%ffn%W(:,1))
    RGBmin(2)=minval(this%ffn%W(:,2))
    RGBL(2)=maxval(this%ffn%W(:,2))
    RGBmin(3)=minval(this%ffn%W(:,3))
    RGBL(3)=maxval(this%ffn%W(:,3))
    RGBL=RGBL-RGBmin
    do j=1,50 !dump 50 frames of initial state
       do i=1,this%ffn%layer%N
          !transform weights into rgb coords from 0-255
          write(222,*)this%nodeCoord(i,1:2),int(256*(this%ffn%W(i,1:3)-RGBmin(1:3))/RGBL(1:3))
       end do
       write(222,*)!new data block for every epoch
    end do
    !train som
    this%sigmastart=2.0
    this%sigmadecay=0.0002
    this%learningstart=1E-4
    this%learningdecay=0.00005
    nepoch=10000!0
!!$    this%sigmastart=2.0
!!$    this%sigmadecay=0.00001
!!$    this%learningstart=5E-1
!!$    this%learningdecay=0.0000001
!!$    nepoch=1000!0
    do epoch=0,nepoch
       call progress(epoch,nepoch)
       !k=ceiling(uran()*nsubject)
       k=mod(epoch,nsubject)+1
       do i=1,nrecordpersubject(k)
          sourcelayer%input=trainingbatch(nx-2:nx,i,k)
          !sourcelayer%input=dataset(nx-2:nx,i)
          !j=ceiling(uran()*size(dataset,2))
          !sourcelayer%input=dataset(nx-2:nx,j)
          call update(sourcelayer)
          call trainstep(this,epoch)
       end do
       if(mod(epoch,10).EQ.0)then !dump data for animation
          !get rgb domain of ffn weights
          RGBmin(1)=minval(this%ffn%W(:,1))
          RGBL(1)=maxval(this%ffn%W(:,1))
          RGBmin(2)=minval(this%ffn%W(:,2))
          RGBL(2)=maxval(this%ffn%W(:,2))
          RGBmin(3)=minval(this%ffn%W(:,3))
          RGBL(3)=maxval(this%ffn%W(:,3))
          RGBL=RGBL-RGBmin
          do k=1,this%ffn%layer%N
             !transform weights into rgb coords from 0-255
             write(222,*)this%nodeCoord(k,1:2),int(256*(this%ffn%W(k,1:3)-RGBmin(1:3))/RGBL(1:3))
          end do
          write(222,*)!new data block for every animation frame
       end if
    end do
    !get rgb domain of ffn weights
    RGBmin(1)=minval(this%ffn%W(:,1))
    RGBL(1)=maxval(this%ffn%W(:,1))
    RGBmin(2)=minval(this%ffn%W(:,2))
    RGBL(2)=maxval(this%ffn%W(:,2))
    RGBmin(3)=minval(this%ffn%W(:,3))
    RGBL(3)=maxval(this%ffn%W(:,3))
    RGBL=RGBL-RGBmin
    do j=1,50 !dump 50 frames of final state
       do i=1,this%ffn%layer%N
          !transform weights into rgb coords from 0-255
          write(222,*)this%nodeCoord(i,1:2),int(256*(this%ffn%W(i,1:3)-RGBmin(1:3))/RGBL(1:3))
       end do
       write(222,*)!new data block for every epoch
    end do
    !close SOM evolution output file
    close(222)
    !animate SOM evolution
    write(*,*)'animating SOM evolution and saving to file anim.gif'
    call system('gnuplot SOMperceptionanim.plt')

    !check neighboring gcells point to nearby coordinates in perception space
    !compared to average distance of corrdiantes by any pair of gcells
    !
    !calculate expected cell pair pointer distance
    R0=0.0_double
    R0var=0.0_double
    do i=1,this%ffn%layer%N
       do j=i,this%ffn%layer%N
          if(i.NE.j)then
             !accumulate distance
             x=sum((this%ffn%W(i,:)-this%ffn%W(j,:))**2)
             R0=R0+sqrt(x)
             R0var=R0var+x
          end if
       end do
    end do
    norm=.5_double*(this%ffn%layer%N**2-this%ffn%layer%N)
    R0=R0/norm
    R0var=R0var/norm-R0*R0
    write(*,*)'average gcell pointer distance and stdev',R0,sqrt(R0var)
    !Group A -compare cells around cell 09: 02,03,08,09,10,16,17 
    RA=0.0_double
    RA=RA+sqrt(sum((this%ffn%W(2,:)-this%ffn%W(9,:))**2))
    RA=RA+sqrt(sum((this%ffn%W(3,:)-this%ffn%W(9,:))**2))
    RA=RA+sqrt(sum((this%ffn%W(8,:)-this%ffn%W(9,:))**2))
    RA=RA+sqrt(sum((this%ffn%W(10,:)-this%ffn%W(9,:))**2))
    RA=RA+sqrt(sum((this%ffn%W(16,:)-this%ffn%W(9,:))**2))
    RA=RA+sqrt(sum((this%ffn%W(17,:)-this%ffn%W(9,:))**2))
    RA=RA/6.0_double
    write(*,*)'average group A gcell pointer distance',RA,RA/R0
    !Group B -compare cells around cell 20: 12,13,19,20,21,26,27
    RB=0.0_double
    RB=RB+sqrt(sum((this%ffn%W(12,:)-this%ffn%W(20,:))**2))
    RB=RB+sqrt(sum((this%ffn%W(13,:)-this%ffn%W(20,:))**2))
    RB=RB+sqrt(sum((this%ffn%W(19,:)-this%ffn%W(20,:))**2))
    RB=RB+sqrt(sum((this%ffn%W(21,:)-this%ffn%W(20,:))**2))
    RB=RB+sqrt(sum((this%ffn%W(26,:)-this%ffn%W(20,:))**2))
    RB=RB+sqrt(sum((this%ffn%W(27,:)-this%ffn%W(20,:))**2))
    RB=RB/6.0_double
    write(*,*)'average group B gcell pointer distance',RB,RB/R0
    !Group C -compare cells around cell 38: 31,32,37,38,39,45,46
    RC=0.0_double
    RC=RC+sqrt(sum((this%ffn%W(31,:)-this%ffn%W(38,:))**2))
    RC=RC+sqrt(sum((this%ffn%W(32,:)-this%ffn%W(38,:))**2))
    RC=RC+sqrt(sum((this%ffn%W(37,:)-this%ffn%W(38,:))**2))
    RC=RC+sqrt(sum((this%ffn%W(39,:)-this%ffn%W(38,:))**2))
    RC=RC+sqrt(sum((this%ffn%W(45,:)-this%ffn%W(38,:))**2))
    RC=RC+sqrt(sum((this%ffn%W(46,:)-this%ffn%W(38,:))**2))
    RC=RC/6.0_double
    write(*,*)'average group C gcell pointer distance',RC,RC/R0

    !sigmadecay is too fast if neighbors point more than 50% of expected pair distance  
    !call assert(RA/R0.LT.0.5,msg='group A nodes point more than 50% of expected pair distance.')
    !call assert(RB/R0.LT.0.5,msg='group B nodes point more than 50% of expected pair distance.')
    !call assert(RC/R0.LT.0.5,msg='group C nodes point more than 50% of expected pair distance.')
    call assert((RA+RB+RC)/R0.LT.1.5,msg='group nodes point more than 50% of expected pair distance.')

    !
    !check that expected gcell pointer distance is the same magnitude as expected dataset distance  
    !calculate expected dataset pair distance
    RD=0.0_double
    RDvar=0.0_double
    do k=1,nsubject
       do kprime=1,nsubject
          do i=1,nrecordpersubject(k)
             do j=1,nrecordpersubject(kprime)
                !if(.not.((i.EQ.j).and.(k.EQ.kprime)))then
                !accumulate distance
                x=sum((trainingbatch(nx-2:nx,i,k)-trainingbatch(nx-2:nx,j,kprime))**2)
                RD=RD+sqrt(x)
                RDvar=RDvar+x
             end do
          end do
       end do
    end do
    norm=real(sum(nrecordpersubject),double)
    norm=norm*norm
    RD=RD/norm
    RDvar=RDvar/norm-RD*RD
    write(*,*)'dataset average distance and stdev',RD,sqrt(RDvar)
    call assert((RD-sqrt(RDvar).lt.R0).and.(RD+sqrt(RDvar).gt.R0)&
         ,msg='average gcell pointer distance does not fit within one sigma of expected dataset distance&
         &, adjust learning rate')

    !cleanup training batches memory
    if(allocated(trainingbatch))deallocate(trainingbatch)
    !cleanup number of records per known subjects array memory
    if(allocated(nrecordpersubject))deallocate(nrecordpersubject)
    !cleanup som memory
    call kill(this)
    call kill(sourcelayer)
    !deallocate PCA eigenvector memory
    if(allocated(PCAeigenvector))deallocate(PCAeigenvector)
    if(allocated(PCAeigenvalue))deallocate(PCAeigenvalue)
    !cleanup dataset mean and variance memory
    if(allocated(datasetbar))deallocate(datasetbar)
    if(allocated(datasetvar))deallocate(datasetvar)
    !!cleanup dataset memory
    !if(allocated(dataset))deallocate(dataset)

!!!!!!!!!!!!!!!!!

    write(*,*)'test training som with all PCA perception data improves model fitness&
         & on set of 49 classifiers on 2D hcp map.'
    call make (sourcelayer,N=nx)
    call make(this,N=49,ND1=7,ND2=7,hcp=.true.,activation='cauchy')
    call link(this%ffn,sourcelayer)

    !get list of known subjects
    open(123,file='knownsubjects.dat')
    ierr=0
    nsubject=0
    subjectlist=0
    do while (ierr.EQ.0)
       nsubject=nsubject+1
       read(123,*,iostat=ierr) subjectlist(nsubject)
    end do
    close(123)
    !get list of known CIDs
    open(123,file='knownCIDs.dat')
    ierr=0
    nCID=0
    CIDlist=0
    do while (ierr.EQ.0)
       nCID=nCID+1
       read(123,*,iostat=ierr) CIDlist(nCID)
    end do
    close(123)
    !allocate number of records per known subjects array
    if(allocated(nrecordpersubject))deallocate(nrecordpersubject)
    allocate(nrecordpersubject(nsubject))
    nrecordpersubject=0
    !allocate dataset mean and variance 
    if(allocated(datasetbar))deallocate(datasetbar)
    allocate(datasetbar(nx))
    datasetbar=0.0_double
    if(allocated(datasetvar))deallocate(datasetvar)
    allocate(datasetvar(nx))
    datasetvar=0.0_double
    !allocate training batches
    nx=24
    if(allocated(trainingbatch))deallocate(trainingbatch)
    allocate(trainingbatch(nx,nmaxrecord,nsubject))
    if(allocated(trainingbatchlabel))deallocate(trainingbatchlabel)
    allocate(trainingbatchlabel(4,nmaxrecord,nsubject))
    !get training data as batches of known subjects
    open(111,file='data/allsubjectperception.dat')
    ierr=0
    nrecordpersubject=0
    do while(ierr.EQ.0)
       read(111,*,iostat=ierr)dummy1,dummy2,dummy3,dummy4,datasetbar
       do i=1,nsubject
          if (subjectlist(i).EQ.dummy1.and.any(CIDlist.EQ.int(dummy2)))then
             nrecordpersubject(i)=nrecordpersubject(i)+1
             call assert(nrecordpersubject(i).LE.nmaxrecord,msg='nrecord per subject is greater than nmaxrecord')
             trainingbatch(:,nrecordpersubject(i),i)=datasetbar
             trainingbatchlabel(1,nrecordpersubject(i),i)=dummy1
             trainingbatchlabel(2,nrecordpersubject(i),i)=dummy2
             trainingbatchlabel(3,nrecordpersubject(i),i)=dummy3
             trainingbatchlabel(4,nrecordpersubject(i),i)=dummy4
          end if
       end do
    end do
    close(111)

    !calculate dataset mean and variance per column
    datasetbar=0._double
    datasetbar=0._double
    do i=1,nsubject
       do j=1,nrecordpersubject(i)
          datasetbar=datasetbar+trainingbatch(:,j,i)
          datasetvar=datasetvar+trainingbatch(:,j,i)**2
       end do
    end do
    norm=real(sum(nrecordpersubject),double)
    datasetbar=datasetbar/norm
    datasetvar=sqrt(datasetvar/norm-datasetbar**2)

    !allocate PCA eigenvecotor and eignevalue memory
    if(allocated(PCAeigenvector))deallocate(PCAeigenvector)
    allocate(PCAeigenvector(nx,nx))
    if(allocated(PCAeigenvalue))deallocate(PCAeigenvalue)
    allocate(PCAeigenvalue(nx))
    if(allocated(PCAwork))deallocate(PCAwork)
    pcaworksizemax=3*nx
    allocate(PCAwork(pcaworksizemax))

    !calculate correlation matrix
    !write(*,*)'calculating perception correlation matrix'
    PCAeigenvector=0._double
    do k=1,nsubject
       do i=1,nx
          do j=i,nx
             PCAeigenvector(i,j)=PCAeigenvector(i,j)&
                  +sum((trainingbatch(i,1:nrecordpersubject(k),k)-datasetbar(i))&
                  *(trainingbatch(j,1:nrecordpersubject(k),k)-datasetbar(j)))&
                  /(datasetvar(i)*datasetvar(j))
             !PCAeigenvector(j,i)=PCAeigenvector(i,j)
          end do
       end do
    end do
    norm=real(sum(nrecordpersubject),double)
    PCAeigenvector=PCAeigenvector/norm

    !query the optimal workspace
    pcaworksize=-1
    call dsyev('V','U',nx,PCAeigenvector,nx,PCAeigenvalue,pcawork,pcaworksize,i)
    pcaworksize = min(pcaworksizemax,int(pcawork(1)))
    !compute eigenvalues and eigenvectors of a real symmetric matrix
    call dsyev('V','U',nx,PCAeigenvector,nx,PCAeigenvalue,pcawork,pcaworksize,i)
    !check for convergence
    call assert(i.EQ.0,msg='lapack routine dsyev failed to compute eigenvalues')

    !rotate perception data
    do k=1,nsubject
       do j=1,nrecordpersubject(k)
          trainingbatch(:,j,k)=matmul(transpose(pcaeigenvector),trainingbatch(:,j,k))
       end do
    end do

    !allocate activity array memory
    if(allocated(activity))deallocate(activity)
    allocate(activity(this%ffn%layer%N))

    !calculate initial model fitness
    fitness=0.0_double
    count=0
    !loop over subject pairs
    do i=1,nsubject
       do j=i,nsubject
          !loop over compound labels
          do k=1,nrecordpersubject(i)
             do kprime=1,nrecordpersubject(j)
                !match for same coupound ID and concentration
                if(trainingbatchlabel(2,k,i).EQ.trainingbatchlabel(2,kprime,j)&
                     .and.&
                     trainingbatchlabel(3,k,i).EQ.trainingbatchlabel(3,kprime,j))then
                   !increment sample counter
                   count=count+1
                   !write(*,*)count,fitness/real(count)

                   !load subject k data
                   sourcelayer%input=trainingbatch(:,k,i)
                   call update(sourcelayer)
                   !calculate gcell activation 
                   call update(this%ffn,weightcentered=.true.)
                   !save subject k activation density
                   activity=this%ffn%layer%node
                   !write(*,*)activity

                   !load subject kprime data
                   sourcelayer%input=trainingbatch(:,k,i)
                   call update(sourcelayer)
                   !calculate gcell activation 
                   call update(this%ffn,weightcentered=.true.)
                   !calculate overlap of subject kprime activation density with that of subject k
                   fitness=fitness+sum(activity*this%ffn%layer%node)

                end if
             end do
          end do
       end do
    end do
    !fitness=fitness*2.0_double/real(nsubject**2-nsubject,double)
    fitness=fitness/real(count,double)
    Finit=fitness
    write(*,*)'initial fitness',Finit

    !train som
    this%sigmastart=2.0
    this%sigmadecay=0.0002
    this%learningstart=1E-4
    this%learningdecay=0.00005
    nepoch=10000!0
!!$    this%sigmastart=2.0
!!$    this%sigmadecay=0.00001
!!$    this%learningstart=5E-1
!!$    this%learningdecay=0.0000001
!!$    nepoch=1000!0
    do epoch=0,nepoch
       call progress(epoch,nepoch)
       !k=ceiling(uran()*nsubject)
       k=mod(epoch,nsubject)+1
       do i=1,nrecordpersubject(k)
          sourcelayer%input=trainingbatch(:,i,k)
          call update(sourcelayer)
          call trainstep(this,epoch)
       end do
       if(mod(epoch,10).EQ.0)then !dump intermediate model fitness
          fitness=0.0_double
          count=0
          !loop over subject pairs
          do i=1,nsubject
             do j=i,nsubject
                !loop over compound labels
                do k=1,nrecordpersubject(i)
                   do kprime=1,nrecordpersubject(j)
                      !match for same coupound ID and concentration
                      if(trainingbatchlabel(2,k,i).EQ.trainingbatchlabel(2,kprime,j)&
                           .and.&
                           trainingbatchlabel(3,k,i).EQ.trainingbatchlabel(3,kprime,j))then
                         !increment sample counter
                         count=count+1
                         !write(*,*)count,fitness/real(count)

                         !load subject k data
                         sourcelayer%input=trainingbatch(:,k,i)
                         call update(sourcelayer)
                         !calculate gcell activation 
                         call update(this%ffn,weightcentered=.true.)
                         !save subject k activation density
                         activity=this%ffn%layer%node
                         !write(*,*)activity

                         !load subject kprime data
                         sourcelayer%input=trainingbatch(:,k,i)
                         call update(sourcelayer)
                         !calculate gcell activation 
                         call update(this%ffn,weightcentered=.true.)
                         !calculate overlap of subject kprime activation density with that of subject k
                         fitness=fitness+sum(activity*this%ffn%layer%node)

                      end if
                   end do
                end do
             end do
          end do
          !fitness=fitness*2.0_double/real(nsubject**2-nsubject,double)
          fitness=fitness/real(count,double)
          write(*,*)epoch,'fitness', fitness

       end if
    end do
    !calculate final model fitness
    fitness=0.0_double
    count=0
    !loop over subject pairs
    do i=1,nsubject
       do j=i,nsubject
          !loop over compound labels
          do k=1,nrecordpersubject(i)
             do kprime=1,nrecordpersubject(j)
                !match for same coupound ID and concentration
                if(trainingbatchlabel(2,k,i).EQ.trainingbatchlabel(2,kprime,j)&
                     .and.&
                     trainingbatchlabel(3,k,i).EQ.trainingbatchlabel(3,kprime,j))then
                   !increment sample counter
                   count=count+1
                   !write(*,*)count,fitness/real(count)

                   !load subject k data
                   sourcelayer%input=trainingbatch(:,k,i)
                   call update(sourcelayer)
                   !calculate gcell activation 
                   call update(this%ffn,weightcentered=.true.)
                   !save subject k activation density
                   activity=this%ffn%layer%node
                   !write(*,*)activity

                   !load subject kprime data
                   sourcelayer%input=trainingbatch(:,k,i)
                   call update(sourcelayer)
                   !calculate gcell activation 
                   call update(this%ffn,weightcentered=.true.)
                   !calculate overlap of subject kprime activation density with that of subject k
                   fitness=fitness+sum(activity*this%ffn%layer%node)

                end if
             end do
          end do
       end do
    end do
    !fitness=fitness*2.0_double/real(nsubject**2-nsubject,double)
    fitness=fitness/real(count,double)
    Ffinal=fitness
    write(*,*)'final fitness',Ffinal

    !cleanup activity array memory
    if(allocated(activity))deallocate(activity)

    stop

    !check neighboring gcells point to nearby coordinates in perception space
    !compared to average distance of corrdiantes by any pair of gcells
    !
    !calculate expected cell pair pointer distance
    R0=0.0_double
    R0var=0.0_double
    do i=1,this%ffn%layer%N
       do j=i,this%ffn%layer%N
          if(i.NE.j)then
             !accumulate distance
             x=sum((this%ffn%W(i,:)-this%ffn%W(j,:))**2)
             R0=R0+sqrt(x)
             R0var=R0var+x
          end if
       end do
    end do
    norm=.5_double*(this%ffn%layer%N**2-this%ffn%layer%N)
    R0=R0/norm
    R0var=R0var/norm-R0*R0
    write(*,*)'average gcell pointer distance and stdev',R0,sqrt(R0var)
    !Group A -compare cells around cell 09: 02,03,08,09,10,16,17 
    RA=0.0_double
    RA=RA+sqrt(sum((this%ffn%W(2,:)-this%ffn%W(9,:))**2))
    RA=RA+sqrt(sum((this%ffn%W(3,:)-this%ffn%W(9,:))**2))
    RA=RA+sqrt(sum((this%ffn%W(8,:)-this%ffn%W(9,:))**2))
    RA=RA+sqrt(sum((this%ffn%W(10,:)-this%ffn%W(9,:))**2))
    RA=RA+sqrt(sum((this%ffn%W(16,:)-this%ffn%W(9,:))**2))
    RA=RA+sqrt(sum((this%ffn%W(17,:)-this%ffn%W(9,:))**2))
    RA=RA/6.0_double
    write(*,*)'average group A gcell pointer distance',RA,RA/R0
    !Group B -compare cells around cell 20: 12,13,19,20,21,26,27
    RB=0.0_double
    RB=RB+sqrt(sum((this%ffn%W(12,:)-this%ffn%W(20,:))**2))
    RB=RB+sqrt(sum((this%ffn%W(13,:)-this%ffn%W(20,:))**2))
    RB=RB+sqrt(sum((this%ffn%W(19,:)-this%ffn%W(20,:))**2))
    RB=RB+sqrt(sum((this%ffn%W(21,:)-this%ffn%W(20,:))**2))
    RB=RB+sqrt(sum((this%ffn%W(26,:)-this%ffn%W(20,:))**2))
    RB=RB+sqrt(sum((this%ffn%W(27,:)-this%ffn%W(20,:))**2))
    RB=RB/6.0_double
    write(*,*)'average group B gcell pointer distance',RB,RB/R0
    !Group C -compare cells around cell 38: 31,32,37,38,39,45,46
    RC=0.0_double
    RC=RC+sqrt(sum((this%ffn%W(31,:)-this%ffn%W(38,:))**2))
    RC=RC+sqrt(sum((this%ffn%W(32,:)-this%ffn%W(38,:))**2))
    RC=RC+sqrt(sum((this%ffn%W(37,:)-this%ffn%W(38,:))**2))
    RC=RC+sqrt(sum((this%ffn%W(39,:)-this%ffn%W(38,:))**2))
    RC=RC+sqrt(sum((this%ffn%W(45,:)-this%ffn%W(38,:))**2))
    RC=RC+sqrt(sum((this%ffn%W(46,:)-this%ffn%W(38,:))**2))
    RC=RC/6.0_double
    write(*,*)'average group C gcell pointer distance',RC,RC/R0

    !sigmadecay is too fast if neighbors point more than 50% of expected pair distance  
    call assert(RA/R0.LT.0.5,msg='group A nodes point more than 50% of expected pair distance.')
    call assert(RB/R0.LT.0.5,msg='group B nodes point more than 50% of expected pair distance.')
    call assert(RC/R0.LT.0.5,msg='group C nodes point more than 50% of expected pair distance.')

    !
    !check that expected gcell pointer distance is the same magnitude as expected dataset distance  
    !calculate expected dataset pair distance
    RD=0.0_double
    RDvar=0.0_double
    do k=1,nsubject
       do kprime=1,nsubject
          do i=1,nrecordpersubject(k)
             do j=1,nrecordpersubject(kprime)
                !accumulate distance
                x=sum((trainingbatch(nx-2:nx,i,k)-trainingbatch(nx-2:nx,j,kprime))**2)
                RD=RD+sqrt(x)
                RDvar=RDvar+x
             end do
          end do
       end do
    end do
    norm=real(sum(nrecordpersubject),double)
    norm=norm*norm
    RD=RD/norm
    RDvar=RDvar/norm-RD*RD
    write(*,*)'dataset average distance and stdev',RD,sqrt(RDvar)

    call assert((RD-sqrt(RDvar).lt.R0).and.(RD+sqrt(RDvar).gt.R0)&
         ,msg='average gcell pointer distance does not fit within one sigma of expected dataset distance&
         &, adjust learning rate')

    !cleanup training batches memory
    if(allocated(trainingbatch))deallocate(trainingbatch)
    if(allocated(trainingbatchlabel))deallocate(trainingbatchlabel)
    !cleanup number of records per known subjects array memory
    if(allocated(nrecordpersubject))deallocate(nrecordpersubject)
    !cleanup som memory
    call kill(this)
    call kill(sourcelayer)
    !deallocate PCA eigenvector memory
    if(allocated(PCAeigenvector))deallocate(PCAeigenvector)
    if(allocated(PCAeigenvalue))deallocate(PCAeigenvalue)
    !cleanup dataset mean and variance memory
    if(allocated(datasetbar))deallocate(datasetbar)
    if(allocated(datasetvar))deallocate(datasetvar)
    !!cleanup dataset memory
    !if(allocated(dataset))deallocate(dataset)
    stop
!!!!!!!!!!!!!!!!!



    stop

    write(*,*)'test training som with all PCA all subject perception data&
         & set evolves 100 classifiers on 2D hcp map toward perception data.'
    !perception data column labels 
    !subjects #,CID,dilution,detection,intensity,pleasantness,familiarity,identification,edible,bakery,sweet,fruit,fish,garlic,spices,cold,sour,burnt,acid,warm,musky,sweaty,ammonia/urinous,decayed,wood,grass,flower,chemical
    !allocate dataset
    nx=24
    if(allocated(dataset))deallocate(dataset)
    allocate(dataset(nx,35090))

    call make (sourcelayer,N=nx)
    call make(this,N=100,ND1=10,ND2=10,hcp=.true.)
    call link(this%ffn,sourcelayer)

    !record training set
    dataset=0
    subjectlist=0
    CIDlist=0
    nsubject=0
    nCID=0
    open(111,file='data/allsubjectperception.dat')
    ierr=0
    i=0
    do while(ierr.EQ.0)
       i=i+1
       read(111,*,iostat=ierr)dummy1,dummy2,dummy3,dummy4,dataset(:,i)

       !count unique CIDs and number of subjects
       if(all(subjectlist.NE.int(dummy1)))then
          nsubject=nsubject+1
          subjectlist(nsubject)=int(dummy1)
       end if
       if(all(CIDlist.NE.int(dummy2)))then
          nCID=nCID+1
          CIDlist(nCID)=int(dummy2)
       end if

    end do
    close(111)

    write(*,*)nsubject,'subjects found'!,subjectlist(1:nsubject)
    write(*,*)nCID,'CIDs found'!,CIDlist(1:nCID)

    !randomly pick 70% of subjects 
    do k=1,1000 !swap values in subjectlist 1000 times
       i=ceiling(uran()*nsubject) !choose random array position to swap
       j=ceiling(uran()*nsubject) !choose random array position to swap
       !save value in position i
       ierr=subjectlist(i)
       !swap j value to position i
       subjectlist(i)=subjectlist(j)
       !return saved value to position j
       subjectlist(j)=ierr
    end do
    open(123,file='knownsubjects.dat')
    do i=1,int(nsubject*.7)
       write(123,*)subjectlist(i)
    end do
    close(123)

    !randomly pick 70% of CIDs 
    do k=1,1000 !swap values in CIDlist 1000 times
       i=ceiling(uran()*nCID) !choose random array position to swap
       j=ceiling(uran()*nCID) !choose random array position to swap
       !save value in position i
       ierr=CIDlist(i)
       !swap j value to position i
       CIDlist(i)=CIDlist(j)
       !return saved value to position j
       CIDlist(j)=ierr
    end do
    open(123,file='knownCIDs.dat')
    do i=1,int(nCID*.7)
       write(123,*)CIDlist(i)
    end do
    close(123)

    !stop
    !allocate dataset mean and variance 
    if(allocated(datasetbar))deallocate(datasetbar)
    allocate(datasetbar(nx))
    datasetbar=0.0_double
    if(allocated(datasetvar))deallocate(datasetvar)
    allocate(datasetvar(nx))
    datasetvar=0.0_double
    !calculate dataset mean and  variance per column
    do i=1,size(dataset,2)
       datasetbar=datasetbar+dataset(:,i)
       datasetvar=datasetvar+dataset(:,i)**2
    end do
    datasetbar=datasetbar/real(size(dataset,2),double)
    datasetvar=sqrt(datasetvar/real(size(dataset,2),double)-datasetbar**2)

    write(*,*)"Perception label      mean                       stdev                    ratio"    
    write(*,*)"intensity      ",datasetbar(1),datasetvar(1)  ,datasetbar(1)/datasetvar(1)  
    write(*,*)"pleasantness   ",datasetbar(2),datasetvar(2)  ,datasetbar(2)/datasetvar(2)  
    write(*,*)"familiarity    ",datasetbar(3),datasetvar(3)  ,datasetbar(3)/datasetvar(3)  
    write(*,*)"identification ",datasetbar(4),datasetvar(4)  ,datasetbar(4)/datasetvar(4)  
    write(*,*)"edible         ",datasetbar(5),datasetvar(5)  ,datasetbar(5)/datasetvar(5)  
    write(*,*)"bakery         ",datasetbar(6),datasetvar(6)  ,datasetbar(6)/datasetvar(6)  
    write(*,*)"sweet          ",datasetbar(7),datasetvar(7)  ,datasetbar(7)/datasetvar(7)  
    write(*,*)"fruit          ",datasetbar(8),datasetvar(8)  ,datasetbar(8)/datasetvar(8)  
    write(*,*)"fish           ",datasetbar(9),datasetvar(9)  ,datasetbar(9)/datasetvar(9)  
    write(*,*)"garlic         ",datasetbar(10),datasetvar(10),datasetbar(10)/datasetvar(10)
    write(*,*)"spices         ",datasetbar(11),datasetvar(11),datasetbar(11)/datasetvar(11)
    write(*,*)"cold           ",datasetbar(12),datasetvar(12),datasetbar(12)/datasetvar(12)
    write(*,*)"sour           ",datasetbar(13),datasetvar(13),datasetbar(13)/datasetvar(13)
    write(*,*)"burnt          ",datasetbar(14),datasetvar(14),datasetbar(14)/datasetvar(14)
    write(*,*)"acid           ",datasetbar(15),datasetvar(15),datasetbar(15)/datasetvar(15)
    write(*,*)"warm           ",datasetbar(16),datasetvar(16),datasetbar(16)/datasetvar(16)
    write(*,*)"musky          ",datasetbar(17),datasetvar(17),datasetbar(17)/datasetvar(17)
    write(*,*)"sweaty         ",datasetbar(18),datasetvar(18),datasetbar(18)/datasetvar(18)
    write(*,*)"ammonia/urinous",datasetbar(19),datasetvar(19),datasetbar(19)/datasetvar(19)
    write(*,*)"decayed        ",datasetbar(20),datasetvar(20),datasetbar(20)/datasetvar(20)
    write(*,*)"wood           ",datasetbar(21),datasetvar(21),datasetbar(21)/datasetvar(21)
    write(*,*)"grass          ",datasetbar(22),datasetvar(22),datasetbar(22)/datasetvar(22)
    write(*,*)"flower         ",datasetbar(23),datasetvar(23),datasetbar(23)/datasetvar(23)
    write(*,*)"chemical       ",datasetbar(24),datasetvar(24),datasetbar(24)/datasetvar(24)

    !allocate PCA eigenvecotor and eignevalue memory
    if(allocated(PCAeigenvector))deallocate(PCAeigenvector)
    allocate(PCAeigenvector(nx,nx))
    if(allocated(PCAeigenvalue))deallocate(PCAeigenvalue)
    allocate(PCAeigenvalue(nx))
    if(allocated(PCAwork))deallocate(PCAwork)
    pcaworksizemax=3*nx
    allocate(PCAwork(pcaworksizemax))

    !calculate PCA
    PCAeigenvector=0._double 
    do i=1,nx
       do j=i,nx
          PCAeigenvector(i,j)=sum((dataset(i,:)-datasetbar(i))*(dataset(j,:)-datasetbar(j)))&
               /(datasetvar(i)*datasetvar(j))
          !PCAeigenvector(j,i)=PCAeigenvector(i,j)
       end do
    end do
    PCAeigenvector=PCAeigenvector/real(size(dataset,2),double)

    !query the optimal workspace
    pcaworksize=-1
    call dsyev('V','U',nx,PCAeigenvector,nx,PCAeigenvalue,pcawork,pcaworksize,i)
    pcaworksize = min(pcaworksizemax,int(pcawork(1)))
    !compute eigenvalues and eigenvectors of a real symmetric matrix
    call dsyev('V','U',nx,PCAeigenvector,nx,PCAeigenvalue,pcawork,pcaworksize,i)
    !check for convergence
    call assert(i.EQ.0,msg='lapack routine dsyev failed to compute eigenvalues')

    write(*,*)'PCA eigen values, percentage, cumulative percentange'
    do i=1,nx
       write(*,*)pcaeigenvalue(i),pcaeigenvalue(i)/sum(pcaeigenvalue)&
            ,sum(pcaeigenvalue(1:i))/sum(pcaeigenvalue),1.0-sum(pcaeigenvalue(1:i))/sum(pcaeigenvalue)
    end do
    write(*,*)

    write(*,*)"Perception label  PCA1 vector,               PCA2 vector,              PCA3 vector"    
    write(*,*)"intensity      ",pcaeigenvector(nx,1),pcaeigenvector(nx-1,1),pcaeigenvector(nx-2,1)
    write(*,*)"pleasantness   ",pcaeigenvector(nx,2),pcaeigenvector(nx-1,2),pcaeigenvector(nx-2,2)
    write(*,*)"familiarity    ",pcaeigenvector(nx,3),pcaeigenvector(nx-1,3),pcaeigenvector(nx-2,3)
    write(*,*)"identification ",pcaeigenvector(nx,4),pcaeigenvector(nx-1,4),pcaeigenvector(nx-2,4)
    write(*,*)"edible         ",pcaeigenvector(nx,5),pcaeigenvector(nx-1,5),pcaeigenvector(nx-2,5)
    write(*,*)"bakery         ",pcaeigenvector(nx,6),pcaeigenvector(nx-1,6),pcaeigenvector(nx-2,6)
    write(*,*)"sweet          ",pcaeigenvector(nx,7),pcaeigenvector(nx-1,7),pcaeigenvector(nx-2,7)
    write(*,*)"fruit          ",pcaeigenvector(nx,8),pcaeigenvector(nx-1,8),pcaeigenvector(nx-2,8)
    write(*,*)"fish           ",pcaeigenvector(nx,9),pcaeigenvector(nx-1,9),pcaeigenvector(nx-2,9)
    write(*,*)"garlic         ",pcaeigenvector(nx,10),pcaeigenvector(nx-1,10),pcaeigenvector(nx-2,10)
    write(*,*)"spices         ",pcaeigenvector(nx,11),pcaeigenvector(nx-1,11),pcaeigenvector(nx-2,11)
    write(*,*)"cold           ",pcaeigenvector(nx,12),pcaeigenvector(nx-1,12),pcaeigenvector(nx-2,12)
    write(*,*)"sour           ",pcaeigenvector(nx,13),pcaeigenvector(nx-1,13),pcaeigenvector(nx-2,13)
    write(*,*)"burnt          ",pcaeigenvector(nx,14),pcaeigenvector(nx-1,14),pcaeigenvector(nx-2,14)
    write(*,*)"acid           ",pcaeigenvector(nx,15),pcaeigenvector(nx-1,15),pcaeigenvector(nx-2,15)
    write(*,*)"warm           ",pcaeigenvector(nx,16),pcaeigenvector(nx-1,16),pcaeigenvector(nx-2,16)
    write(*,*)"musky          ",pcaeigenvector(nx,17),pcaeigenvector(nx-1,17),pcaeigenvector(nx-2,17)
    write(*,*)"sweaty         ",pcaeigenvector(nx,18),pcaeigenvector(nx-1,18),pcaeigenvector(nx-2,18)
    write(*,*)"ammonia/urinous",pcaeigenvector(nx,19),pcaeigenvector(nx-1,19),pcaeigenvector(nx-2,19)
    write(*,*)"decayed        ",pcaeigenvector(nx,20),pcaeigenvector(nx-1,20),pcaeigenvector(nx-2,20)
    write(*,*)"wood           ",pcaeigenvector(nx,21),pcaeigenvector(nx-1,21),pcaeigenvector(nx-2,21)
    write(*,*)"grass          ",pcaeigenvector(nx,22),pcaeigenvector(nx-1,22),pcaeigenvector(nx-2,22)
    write(*,*)"flower         ",pcaeigenvector(nx,23),pcaeigenvector(nx-1,23),pcaeigenvector(nx-2,23)
    write(*,*)"chemical       ",pcaeigenvector(nx,24),pcaeigenvector(nx-1,24),pcaeigenvector(nx-2,24)
    write(*,*)

    !rotate perception data
    do i=1,size(dataset,2)
       dataset(:,i)=matmul(transpose(pcaeigenvector),dataset(:,i))
       !write(*,*) dataset(:,i)
    end do

    !dump som evolution in rgb coordinates
    open(222,file='SOMperceptionRGB-allPCAallSUB.dat')
    !get rgb domain of ffn weights
    RGBmin(1)=minval(this%ffn%W(:,nx-2))
    RGBL(1)=maxval(this%ffn%W(:,nx-2))
    RGBmin(2)=minval(this%ffn%W(:,nx-1))
    RGBL(2)=maxval(this%ffn%W(:,nx-1))
    RGBmin(3)=minval(this%ffn%W(:,nx))
    RGBL(3)=maxval(this%ffn%W(:,nx))
    RGBL=RGBL-RGBmin
    do j=1,50 !dump 50 frames of initial state
       do i=1,this%ffn%layer%N
          !transform weights into rgb coords from 0-255
          write(222,*)this%nodeCoord(i,1:2),int(256*(this%ffn%W(i,nx-2:nx)-RGBmin)/RGBL)
       end do
       write(222,*)!new data block for every epoch
    end do
    !train som
    !    this%sigmastart=3.0
    !    this%sigmadecay=0.002
    !    this%learningstart=0.05
    !    this%learningdecay=0.00001
    !    nepoch=10000
    !    !train som
    !    this%sigmadecay=0.0001
    !    this%learningdecay=0.001
    nepoch=1000
    do epoch=0,nepoch
       call progress(epoch,nepoch)
       do i=1,size(dataset,2)
          sourcelayer%input=dataset(:,i)
          !j=ceiling(uran()*size(dataset,2))
          sourcelayer%input=dataset(:,j)
          call update(sourcelayer)
          call trainstep(this,epoch)
       end do
       !get rgb domain of ffn weights
       RGBmin(1)=minval(this%ffn%W(:,nx-2))
       RGBL(1)=maxval(this%ffn%W(:,nx-2))
       RGBmin(2)=minval(this%ffn%W(:,nx-1))
       RGBL(2)=maxval(this%ffn%W(:,nx-1))
       RGBmin(3)=minval(this%ffn%W(:,nx))
       RGBL(3)=maxval(this%ffn%W(:,nx))
       RGBL=RGBL-RGBmin
       do i=1,this%ffn%layer%N
          !transform weights into rgb coords from 0-255
          write(222,*)this%nodeCoord(i,1:2),int(256*(this%ffn%W(i,nx-2:nx)-RGBmin)/RGBL)
       end do
       write(222,*)!new data block for every epoch
    end do
    !get rgb domain of ffn weights
    RGBmin(1)=minval(this%ffn%W(:,nx-2))
    RGBL(1)=maxval(this%ffn%W(:,nx-2))
    RGBmin(2)=minval(this%ffn%W(:,nx-1))
    RGBL(2)=maxval(this%ffn%W(:,nx-1))
    RGBmin(3)=minval(this%ffn%W(:,nx))
    RGBL(3)=maxval(this%ffn%W(:,nx))
    RGBL=RGBL-RGBmin
    do j=1,50 !dump 50 frames of final state
       do i=1,this%ffn%layer%N
          !transform weights into rgb coords from 0-255
          write(222,*)this%nodeCoord(i,1:2),int(256*(this%ffn%W(i,nx-2:nx)-RGBmin)/RGBL)
       end do
       write(222,*)!new data block for every epoch
    end do
    !close SOM evolution output file
    close(222)

    !cleanup som memory
    call kill(this)
    call kill(sourcelayer)
    !deallocate PCA eigenvector memory
    if(allocated(PCAeigenvector))deallocate(PCAeigenvector)
    if(allocated(PCAeigenvalue))deallocate(PCAeigenvalue)
    !cleanup dataset mean and variance memory
    if(allocated(datasetbar))deallocate(datasetbar)
    if(allocated(datasetvar))deallocate(datasetvar)
    !cleanup dataset memory
    if(allocated(dataset))deallocate(dataset)

  end subroutine som_test
  !-----------------------------------------

end module som_class

