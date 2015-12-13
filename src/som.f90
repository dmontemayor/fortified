!>\brief
!! Class som
!!\details
!! Self organizing map SOM is a type of feed foward network distinguished
!! by a map of interconnected nodes in the map space.
!! weights to external input nodes. 
!<------------------------------------------------------------------------
module SOM_class
  use type_kinds
  use ffn_class
  implicit none
  private

  public::SOM, SOM_test
  public::make, kill, display, store, update, reset, check, trainstep

  type SOM
     logical::initialized=.false.
     type(ffn)::ffn
     integer(long)::ND1,ND2,ND3
     real(double),dimension(:,:),pointer::nodeCoord
     logical::mexhat
     real(double)::sigmastart,sigmadecay
     real(double)::learningstart,learningdecay
  end type SOM

  !> Train SOM object by one step.
  interface trainstep
     module procedure SOM_trainstep
  end interface

  !> Creates the SOM object.
  interface make
     module procedure SOM_init
  end interface

  !> Destroys the SOM object.
  interface kill
     module procedure SOM_kill
  end interface

  !> Displays the current state of the SOM object.
  interface display
     module procedure SOM_display
  end interface

  !> Stores the current state of the SOM object.
  interface store
     module procedure SOM_store
  end interface

  !> Recaluclates the SOM object.
  interface update
     module procedure SOM_update
  end interface

  !> Reinitializes the SOM object.
  interface reset
     module procedure SOM_reset
  end interface

  !> Checks that the SOM object.
  interface check
     module procedure SOM_check
  end interface

contains
  !======================================================================
  !> \brief Evolves the SOM object by one training step.
  !> \param this is the SOM object to be evolved.
  !> \param[in] S is the current training step
  !=====================================================================
  subroutine SOM_trainstep(this,S)
    type(SOM),intent(inout)::this
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

  end subroutine SOM_trainstep

  !======================================================================
  !> \brief SOM learning function
  !> \param this is the SOM object to measure.
  !> \param[in] S is the current training step.
  !=====================================================================
  function Lfunc(this,S)
    use functions
    type(SOM),intent(inout)::this
    integer(long),intent(in)::S
    real(double)::Lfunc
    Lfunc=this%learningstart*exp(-this%learningdecay*S)
    return
  end function Lfunc

  !======================================================================
  !> \brief SOM neigborhood function
  !> \param this is the SOM object to measure.
  !> \param[in] A is the winning node.
  !> \param[in] S is the current training step.
  !=====================================================================
  function Nfunc(this,A,S)
    use functions
    type(SOM),intent(inout)::this
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
  !> \brief cartesian (L2norm) distance between two nodes in SOM.
  !> \param this is the SOM object to measure.
  !> \param[in] A is the starting node.
  !> \param[in] B is the ending node.
  !=====================================================================
  real(double) function cartdist(this,A,B)
    type(SOM),intent(inout)::this
    integer(long),intent(in)::A,B

    cartdist=sqrt(sum((this%nodecoord(A,:)-this%nodecoord(B,:))**2))

    return
  end function cartdist

  !======================================================================
  !> \brief Creates and initializes the SOM object.
  !> \param this is the SOM object to be initialized.
  !> \param[in] file is an optional string containing the name of a previously stored SOM file.
  !> \remark If no input file is provided the user must manually initialize THIS using stout.
  !=====================================================================
  subroutine SOM_init(this,N,ND1,ND2,ND3,hcp)
    type(SOM),intent(inout)::this
    integer(long),intent(in)::N
    integer(long),intent(in),optional::ND1,ND2,ND3

    logical,intent(in),optional::hcp
    real(double)::root3,third,twothirdsroot6

    integer(long)::i,j,k,M

    !set number of nodes on dimension 1
    this%ND1=N
    if(present(ND1))this%ND1=ND1

    !set number of nodes on dimension 2
    this%ND2=1
    if(present(ND2))this%ND2=ND2

    !set number of nodes on dimension 3
    this%ND3=1
    if(present(ND3))this%ND3=ND3

    call make(this%ffn,N=this%ND1*this%ND2*this%ND3)

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
  end subroutine SOM_init

  !======================================================================
  !> \brief Destroys the SOM object.
  !> \param this is the SOM object to be destroyed.
  !====================================================================
  subroutine SOM_kill(this)
    type(SOM),intent(inout)::this
 
    !kill the layer primitive
    call kill(this%ffn)

    !kill node coordinates
    if(associated(this%nodecoord))nullify(this%nodecoord)

    !un-initialize som object
    this%initialized=.false.

  end subroutine SOM_kill

  !======================================================================
  !> \brief Computes the current state of SOM object.
  !> \param this is the SOM  object to be updated.
  !======================================================================
  subroutine SOM_update(this)
    type(SOM),intent(inout)::this

  end subroutine SOM_update

  !======================================================================
  !> \brief Re-initiallizes the SOM object.
  !> \param this is the SOM  object to be re-initialized.
  !======================================================================
  subroutine SOM_reset(this)
    type(SOM),intent(inout)::this

  end subroutine SOM_reset

  !======================================================================
  !> \brief Stores the current state of the SOM object to file.
  !> \param[in] this is the SOM  object to be updated.
  !> \param[in] file is a string containing the location of the store file.
  !======================================================================
  subroutine SOM_store(this,file)
    type(SOM),intent(in)::this
    character*(*),intent(in)::file

  end subroutine SOM_store

  !======================================================================
  !> \brief Retrun the som object as a single line record entry.
  !> \param[in] this is the som object.
  !> \param[in] msg is an optional string message to annotate the displayed object.
  !======================================================================
  character(len=line) function som_display(this,msg)
    type(som),intent(in)::this
    character*(*),intent(in),optional::msg
    character(len=7)::FMT='(A10)'

    write(som_display,FMT)'helloworld'

   
  end function som_display

  !======================================================================
  !> \brief Checks the SOM object.
  !> \param[in] this is the SOM object to be checked.
  !> \return Nothing if all checks pass or 1 and a warn for the first failed check.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function SOM_check(this)
    use testing_class
    type(SOM),intent(in)::this

    !initiate with no problems found 
    SOM_check=0

    !check that object is initialized
    call assert(this%initialized,msg='SOM_check: SOM object not initialized.',iostat=SOM_check)
    if(SOM_check.NE.0)return

    !check the ffn
    call assert(check(this%ffn).EQ.0,msg='SOM_check: failed ffn check!',iostat=SOM_check)
    if(SOM_check.NE.0)return

    !check map dimensions
    call assert(this%ND1*this%ND2*this%ND3.EQ.this%ffn%layer%N,msg='SOM_check: all nodes must fit in rectangular 3D grid.'&
         ,iostat=SOM_check)
    if(SOM_check.NE.0)return

    !check node coordinates
    call assert(associated(this%nodeCoord),msg='SOM_check: node coordinates not associated',iostat=SOM_check)
    if(SOM_check.NE.0)return
    call assert(size(this%nodecoord).EQ.this%ffn%layer%N*3,msg='SOM_check: size of node coordinates not equal to number of nodes.',&
         iostat=SOM_check)
    if(SOM_check.NE.0)return
    call assert(all(this%nodecoord.EQ.this%nodecoord),msg='SOM_check: node coordinates have NaN values.',iostat=SOM_check)
    if(SOM_check.NE.0)return

    !check sigmastart
    call assert(this%sigmastart,this%sigmastart,msg='SOM_check: sigmastart has NaN or huge vaule.',iostat=SOM_check)
    if(SOM_check.NE.0)return
    call assert(this%sigmastart.GT.epsilon(1.0_double),msg='SOM_check: sigmastart has tiny or negative value.',iostat=SOM_check)
    if(SOM_check.NE.0)return

    !check sigmadecay
    call assert(this%sigmadecay,this%sigmadecay,msg='SOM_check: sigmadecay has NaN or huge vaule.',iostat=SOM_check)
    if(SOM_check.NE.0)return
    call assert(this%sigmadecay.GE.0.0_double,msg='SOM_check: sigmadecay has negative value.',iostat=SOM_check)
    if(SOM_check.NE.0)return

    !check learningstart
    call assert(this%learningstart,this%learningstart,msg='SOM_check: learningstart has NaN or huge vaule.',iostat=SOM_check)
    if(SOM_check.NE.0)return
    call assert(this%learningstart.GE.0.0_double,msg='SOM_check: learningstart has negative value.',iostat=SOM_check)
    if(SOM_check.NE.0)return

    !check learningdecay
    call assert(this%learningdecay,this%learningdecay,msg='SOM_check: learningdecay has NaN or huge vaule.',iostat=SOM_check)
    if(SOM_check.NE.0)return
    call assert(this%learningdecay.GE.0.0_double,msg='SOM_check: learningdecay has negative value.',iostat=SOM_check)
    if(SOM_check.NE.0)return

  end function SOM_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the SOM methods.
  !> \param[in] this is the SOM object whose methods will be excercised.
  !> \return Nothing if all tests pass or 1 and a stop for the first failed test.
  !> \remark Will stop after first failed check.
  !======================================================================
  subroutine SOM_test
    use testing_class
    use layer_class
    type(SOM)::this
    type(layer)::sourcelayer

    real(double),allocatable::vector(:),mat(:,:)
    real(double)::x,y,z
    integer(long)::i,j,k,ierr,epoch,nres,count
    character::residue,resname(20),label
    character(len=1000)::seq,labseq

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

    write(*,*)'test SOM training step evolves node toward source node at 1.0'
    call make (sourcelayer,N=1)
    sourcelayer%node=1.0_double
    call make(this,N=1)
    call link(this%ffn,sourcelayer)
    this%ffn%W=0.0_double
    x=0_double
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
    call assert(y.LT.x,msg='SOM node did not evolve node toward source node at 1.0.')
    call kill(this)
    call kill(sourcelayer)


    write(*,*)'test SOM training step evolves node toward source node vector at (-1.0,-1.0)'
    call make (sourcelayer,N=2)
    sourcelayer%node(1:2)=-1.0_double
    call make(this,N=1)
    call link(this%ffn,sourcelayer)
    x=0_double
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
    call assert(y.LT.x,msg='SOM node did not evolve node toward source node vector at (-1.0,-1.0).')
    call kill(this)
    call kill(sourcelayer)


    write(*,*)'test training SOM with doublewell data set evolves 2 labels to points (1,2) and (-1,-2).'
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
    call assert(x.LT.0.16.and.y.LT.0.05,msg='SOM did not evolve well 2 nodes toward the centers of the doublepotential.')
    call kill(this)
    call kill(sourcelayer)
    

    write(*,*)'test training SOM with zpotential data set evolves 10 labels on 1D map toward z potential.'
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
    call assert(z.LT.1.5E-3,msg='SOM zpotential training MSE is greater than 1.5E-3')
    call kill(this)
    call kill(sourcelayer)

!!$stop
!!$    write(*,*)'test training SOM with protstruct dataset to recover 3 labels&
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
!!$                   !load points into SOM and train
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
!!$       read(111,*,iostat=ierr)residue,label
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
!!$                !find nearest SOM label
!!$                do j=1,3
!!$                   !get distance of input coord with SOM node j
!!$                   x=sqrt(sum( (this%ffn%W(j,1:60)-sourcelayer%input(1:60))**2))
!!$                   if(j.EQ.1)then
!!$                      !autoamtic save mindistance z and label k if first label
!!$                      z=x
!!$                      k=j
!!$                   end if
!!$                   if(x.LE.z)then
!!$                      !save mindistance z and label k
!!$                      z=x
!!$                      k=j
!!$                   end if
!!$                end do
!!$                
!!$                !correlate SOM label k with training label j
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
!!$                   write(*,*)'unkown label program will stop. ***'//labseq(i:i)//'***'
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
!!$                labseq=trim(labseq)//label
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
!!$!    call assert(z.LT.1.07E-3,msg='SOM zpotential training MSE is greater than 1.07E-3')
!!$    call kill(this)
!!$    call kill(sourcelayer)


  end subroutine SOM_test
  !-----------------------------------------

end module SOM_class

