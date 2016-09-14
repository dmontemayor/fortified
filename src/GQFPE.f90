!>\brief
!! GQFPE class
!!\details
!! A Generalized quantum Fokker-Plank equation Engine appropriate for
!! photo-induced quantum dynamical processes. The Engine extends the path
!! integral method of Caldeira Legget (Physica A (1983) 121, 587) for a
!! non-equilibrium infulence functional. The Engine expands the paths in 
!! the exponents of the infulence functional to second order with respect
!! to time inorder to treat non-equilibrium and non-Markovian effects
!! consistently.
!<------------------------------------------------------------------------
module GQFPE_class
  use type_kinds
  use math
  use atomicunits
  implicit none
  private

  public::GQFPE, GQFPE_test
  public::describe, check, make, kill, backup, update, reset, status
  public::step

  type GQFPE
     logical::initialized=.false.
     character(len=label)::name='GQFPE'
     
     !***************      Enter GQFPE attributes here     ******************!
     !static parameters
     real(double)::dt         !simulation time increment
     real(double)::time       !current simulation time
     real(double)::omega_c    !cutoff frequency
     real(double)::beta       !inverse temperature
     integer(long)::nstate    !number of oscillator states

     !dynamic parameter pointer arrays
     complex(double),dimension(:,:),pointer::density

     !objects

     !calculated variables
     real(double)::R_pq
     real(double)::R_qq
     real(double)::R_pp
     real(double)::D

     !dynamic variable pointer arrays
     !complex(double),dimension(:,:,:),pointer::supermatrixL
     complex(double),dimension(:,:),pointer::supermatrixL
     complex(double),dimension(:),pointer::expancoef,expancoefdot

     !***********************************************************************!
     !=======================================================================!
     !***   Here are some example attributes your derived type may have   ***!
     ! type(object)::that                              !an other type        !
     ! integer(short)::ierr                            !a short integer      !
     ! integer(long)::ndim                             !a long integer       !
     ! real(double)::var                               !a real variable      !
     ! complex(double)::zed                            !a complex variable   !
     ! real(double),dimension(:,:),pointer::matrix     !a real matrix        !
     ! complex(double),dimension(:,:),pointer::Zmatrix !a complex matrix     !
     !***********************************************************************!
  end type GQFPE

  !> Runs engine for one integration step.
  interface step
     module procedure GQFPE_step
  end interface

  !> Creates the GQFPE object.
  interface make
     module procedure GQFPE_init
  end interface

  !> Destroys the GQFPE object.
  interface kill
     module procedure GQFPE_kill
  end interface

  !> Returns current state of the GQFPE object.
  interface status
     module procedure GQFPE_status
  end interface

  !> Returns a plain text description of the GQFPE object.
  interface describe
     module procedure GQFPE_describe
  end interface
  
  !> Returns the current state of the GQFPE object.
  interface backup
     module procedure GQFPE_backup
  end interface

  !> Recaluclates the GQFPE object.
  interface update
     module procedure GQFPE_update
  end interface

  !> Reinitializes the GQFPE object.
  interface reset
     module procedure GQFPE_reset
  end interface

  !> Checks that the GQFPE object.
  interface check
     module procedure GQFPE_check
  end interface

contains
  !======================================================================
  !> \brief Run GQFPE engine for one integration step.
  !! \param THIS is the GQFPE object.
  !=====================================================================
  subroutine GQFPE_step(this)
    type(GQFPE),intent(inout)::this
    integer(long)::N,INFO
    integer(long)::IPIV(size(this%supermatrixL,1))
    complex(double)::WORK(size(this%supermatrixL,1))
    complex(double)::prop(size(this%supermatrixL,1),size(this%supermatrixL,2))
    complex(double)::prevrate(size(this%supermatrixL,1))

!!$    !update to get current supermatrixL(t)
!!$    call update(this)
!!$    !square current supermatrixL(t)
!!$    prop=matmul(this%supermatrixL,this%supermatrixL)
!!$    !increment by half time step to time=t+dt/2
!!$    this%time=this%time+0.5_double*this%dt
!!$    !update to get supermatrixL(t+dt/2)
!!$    call update(this)
!!$    !accumulate propagator (1+dt*L(t+dt/2)+L(t)**2)
!!$    prop=iden(this%nstate**2)+prop+this%dt*this%supermatrixL
!!$    !increment by half time step to get currtime=t+dt
!!$    this%time=this%time+0.5_double*this%dt
!!$    !evolve expansion coeficients: f(t+dt)=(1+dt**L(t+dt/2)+L(t)**2) * f(t)
!!$    this%expancoef=matmul(prop,this%expancoef)
!!$

    !save previous expansion coeficient rate for prevtime=t-dt
    prevrate=this%expancoefdot
    !increment by time step to time=t+dt
    this%time=this%time+this%dt
    !update to get current supermatrixL(t+dt) and previous expansion coef rates(t)
    call update(this)
    !calculate propagator for time=t+dt
    prop=iden(this%nstate**2)-0.5_double*this%supermatrixL
    !use cgetrf followed by cgetri to invert propagator
    N=this%nstate**2
    call cgetrf(N,N,prop,N,IPIV,INFO)
    call cgetri(N,prop,N,IPIV,WORK,N,INFO)
    !!evolve expansion coeficients f(t+dt)=inverse_prop*(f(t)+frate(t-dt)/2)
    !this%expancoef=this%expancoef+0.5_double*prevrate
    !evolve expansion coeficients f(t+dt)=inverse_prop*(f(t)+frate(t)/2)
    this%expancoef=this%expancoef+0.5_double*this%expancoefdot
    this%expancoef=matmul(prop,this%expancoef)

    !write(*,*)'runtime=',this%time*this%omega_c

  end subroutine GQFPE_step
  !======================================================================
  !> \brief Retruns a description of GQFPE as a string.
  !> \param[in] THIS is the GQFPE object.
  !======================================================================
  character(len=comment) function GQFPE_describe(this)
    type(GQFPE),intent(in)::this
    character(len=5)::FMT='(A)'

    write(GQFPE_describe,FMT)'No description for GQFPE has been provided.'
   
  end function GQFPE_describe

  !======================================================================
  !> \brief Creates and initializes GQFPE.
  !! \param THIS is the GQFPE object.
  !! \param[in] FILE is an optional string containing the name of a
  !! previously backuped GQFPE file.
  !=====================================================================
  subroutine GQFPE_init(this,file)
    use filemanager
    use testing_class
    type(GQFPE),intent(inout)::this
    character*(*),intent(in),optional::file
    integer(long)::unit
    logical::fileisopen=.false.
    character(len=label)::header
    character(len=path)::infile

    integer(long)::i,j,k,l

    !initialize all sub-objects
    !***      example     ***
    !call make(this%object)
    !************************

    !check input file
    if(present(file))then 
       
       !check input file
       inquire(file=file,opened=fileisopen,number=unit)
       if(unit.LT.0)unit=newunit()
       if(.not.fileisopen)open(unit,file=file)
       
       !check if file is of type GQFPE
       read(unit,*)header
       call assert(trim(header).EQ.this%name,msg='GQFPE_init: bad input file header in file'//file)
       
       !read static parameters
       !***             example            ***
       !***read scalar parameters from file***
       !read(unit,*)this%XXX
       !**************************************
       read(unit,*)this%dt
       read(unit,*)this%time
       read(unit,*)this%omega_c
       read(unit,*)this%beta
       read(unit,*)this%nstate

       !use reset to manage dynamic memory, reset sub-objects, and set random parameters
       call reset(this,state=1)
       
       !READ dynamic parameter array values
       !***      example     ***
       !read(unit,*)(this%PPP(i),i=0,this%XXX-1)
       !************************
       read(unit,*)((this%density(i,j)&
            ,j=0,size(this%density,2)-1)&
            ,i=0,size(this%density,1)-1)
       
       !READ sub-objects
       !***      example     ***
       !read(unit,*)infile
       !call make(this%object,file=trim(infile))
       !************************
       
       !finished reading all attributes - now close backup file
       if(.not.fileisopen)close(unit)
       
       !declare initialization complete
       this%initialized=.true.
    else
       !Set static parameters to default settings
       !***       example      ***
       !***set scalar parameter***
       !this%XXX=123
       !**************************
       this%dt=fs
       this%time=0.0_double
       this%omega_c=200.0_double*invcm       
       this%beta=5.0_double/this%omega_c       
       this%nstate=2

       !Use reset to make a default object
       call reset(this,state=1)
    end if

    !finished making now update object
    call update(this)

  end subroutine GQFPE_init

  !======================================================================
  !> \brief Destroys the GQFPE object.
  !> \param THIS is the GQFPE object to be destroyed.
  !> \remarks kill is simply the reset method passed with a null flag 
  !====================================================================
  subroutine GQFPE_kill(this)
    type(GQFPE),intent(inout)::this
 
    call reset(this,0)

  end subroutine GQFPE_kill

  !======================================================================
  !> \brief Computes the current state of GQFPE object.
  !> \param THIS is the GQFPE  object to be updated.
  !======================================================================
  subroutine GQFPE_update(this)
    type(GQFPE),intent(inout)::this

    !parameters
    real(double),parameter::gamma_s=0.1_double
    real(double),parameter::omega_s=1.0_double
    integer(long),parameter::MatsubaraNmin=10 !minimum matsubara number

    !variables
    real(double)::beta_s
    real(double)::t_s
    real(double)::massratio,Gammaoft 
    real(double)::A1,A2
    real(double)::B1,B2
    real(double)::Lterm,LRe,LIm
    real(double)::twopioverbeta,cotbeta,expdts
    real(double)::D,D2,Dts,D2ts,term
    real(double)::Gzero,Gonetilde,Gtwo

    !time saving values
    real(double)::omegaA1,omegaB1
    real(double)::omegaA2,omegaB2
    real(double)::gammaGamma,gammaRpq

    !book keeping
    integer(long)::istate,jstate
    integer(long)::lstate,lpstate
    integer(long)::n,i,j,k,kp,statepair(2)



    !variables
    t_s=this%time*this%omega_c
    beta_s=this%beta*hbar*this%omega_c
    massratio=1.0_double-2.0_double*gamma_s*F2(t_s) !verified
    Gammaoft=F1tilde(t_s)/massratio   !verified
    A1=(1.0_double+massratio)/(2.0_double*massratio)&
         +gamma_s*exp(-t_s)/omega_s**2
    A2=(1.0_double-massratio)/(2.0_double*massratio)&
         -gamma_s*exp(-t_s)/omega_s**2

    twopioverbeta=twopi/beta_s
    cotbeta=1.0_double/tan(0.5_double*beta_s)

    !accumulate G functions over minimum Matsubara numbers
    Gzero=(0.0_double,0.0_double)
    Gonetilde=(0.0_double,0.0_double)
    Gtwo=(0.0_double,0.0_double)
    D=0.0_double
    do n=1,MatsubaraNmin
       D=D+twopioverbeta
       D2=D*D
       Dts=D*t_s
       D2ts=Dts*Dts
       expdts=exp(-Dts)

       !term=F2(D*t_s)/(D2*(D2-1))
       Gtwo=Gtwo+(1.0_double-(1.0_double+Dts+0.5_double*D2ts)*expdts)&
            /(D2*(D2-1))
       !term=F1tilde(D*t_s)/(D*(D2-1))
       Gonetilde=Gonetilde+(1.0_double-(1.0_double+Dts-0.5_double*D2ts)*expdts)&
            /(D*(D2-1))
       term=(1.0_double-expdts)/(D2-1)
       !term=F0(D*t_s)/(D2-1)
       Gzero=Gzero+term
         !if(t_s>1.0) write(444,*)Gzero
    end do
    !accumulate higher matsubara terms until G0 converges defined when new term
    ! contributes less than 0.1% of current G0.
    !Note, other G-terms converge faster than G0 so we need only focus
    ! on G0 convergence to assume global convergence. 
    do while(term/Gzero.GT.5E-4)
       D=D+twopioverbeta
       D2=D*D
       Dts=D*t_s
       D2ts=Dts*Dts
       expdts=exp(-Dts)
       
       !term=F2(D*t_s)/(D2*(D2-1))
       Gtwo=Gtwo+(1.0_double-(1.0_double+Dts+0.5_double*D2ts)*expdts)&
            /(D2*(D2-1))
       !term=F1tilde(D*t_s)/(D*(D2-1))
       Gonetilde=Gonetilde+(1.0_double-(1.0_double+Dts-0.5_double*D2ts)*expdts)&
            /(D*(D2-1))
       term=(1.0_double-expdts)/(D2-1)
       !term=F0(D*t_s)/(D2-1)
       Gzero=Gzero+term
         !if(t_s>1.0) write(444,*)Gzero
    end do
      !if(t_s>1.0)flush(444)
      !if(t_S>1.0)stop

    !multiply by units
    Gzero=Gzero*4.0_double/beta_s
    Gonetilde=Gonetilde*4.0_double/beta_s
    Gtwo=Gtwo*4.0_double/beta_s
    !add high temperature term
    Gzero=Gzero+F0(t_s)*cotbeta
    Gonetilde=Gonetilde+F1tilde(t_s)*cotbeta
    Gtwo=Gtwo+F2(t_s)*cotbeta
    
    !calculate R_pq
    this%R_pq=Gonetilde/massratio&
         -2.0_double*gamma_s*Gtwo&
         *F1tilde(t_s)/(massratio*massratio)

    !calculate R_qq
    this%R_qq=0.5_double*gamma_s*Gtwo/(massratio*massratio)

    !calculate R_pp
    this%R_pp=0.0_double !set default value
    this%R_pp=Gzero/gamma_s&
         -2.0_double*Gammaoft*Gonetilde&
         +2.0_double*gamma_s*Gammaoft*Gammaoft*Gtwo

!!$    !calculate R_pq
!!$    this%R_pq=0.0_double !set default value
!!$    this%R_pq=G1tilde(beta_s,t_s)/massratio&
!!$         -2.0_double*gamma_s*G2(beta_s,t_s)&
!!$         *F1tilde(t_s)/(massratio*massratio)
!!$
!!$    !calculate R_qq
!!$    this%R_qq=0.0_double !set default value
!!$    this%R_qq=0.5_double*gamma_s*G2(beta_s,t_s)/(massratio*massratio)
!!$
!!$    !calculate R_pp
!!$    this%R_pp=0.0_double !set default value
!!$    this%R_pp=G0(beta_s,t_s)/gamma_s&
!!$         -2.0_double*Gammaoft*G1tilde(beta_s,t_s)&
!!$         +2.0_double*gamma_s*Gammaoft*Gammaoft*G2(beta_s,t_s)

    !calculate D inequality term
    this%D=0.0_double !set default value
    this%D=4.0_double*this%R_pp*this%R_qq&
         -this%R_pq*this%R_pq&
         -Gammaoft*Gammaoft

    !calculate B1 and B2
    B2=this%R_pp*(gamma_s/omega_s)**2
    B1=this%R_qq+B2
    B2=this%R_qq-B2

    !calculate time saving values
    omegaA1=omega_s*A1
    omegaB1=omega_s*B1
    omegaA2=omega_s*A2
    omegaB2=omega_s*B2
    gammaGamma=gamma_s*Gammaoft
    gammaRpq=gamma_s*this%R_pq


    !Save previous expansion coeficient rate for old supermatrixL
    this%expancoefdot=matmul(this%supermatrixL,this%expancoef)

    !calculate supermatrixL
    this%supermatrixL=(0.0_double,0.0_double) !set default value

    !loop over Liouville states (lstate) defined as schrodinger state pair
    ! coordinates (istate,jstate) following a space filling curve. Space filling
    ! curve type options are 'rowmajor' and 'hilbert'
    do lstate=0,this%nstate*this%nstate-1
       !calculate schrodinger state pair
       statepair=spacefilling_index2coord(this%nstate,lstate,type='rowmajor')
       !statepair=spacefilling_index2coord(this%nstate,lstate,type='hilbert')
       
       istate=statepair(1)
       jstate=statepair(2)
       
       !---------------------accumulate all 9 terms---------------------
       !---------with i=(istate-ipstate) and j=(jstate-jpstate)---------
       !---------note: modulo function always returns positive----------
       !------note: these 9 terms are accessed by 3x3 index (k,kp)------ 
       ! i = -2   j = 0 
       ! i = -1   j = -1, +1 
       ! i = +0   j = -2, 0, +2
       ! i = +1   j = -1, +1
       ! i = +2   j = 0
       !----------------------------------------------------------------
       
       !acquire 9 terms
       do k=0,2
          do kp=0,2
             !calculate schrodinger prime-state pair
             i=k+kp-2
             j=k-kp
             statepair(1)=modulo(istate-i,this%nstate)
             statepair(2)=modulo(jstate-j,this%nstate)
             !get Loiuville prime-state (lpstate) for current prime-state pair
             lpstate=spacefilling_coord2index(this%nstate,statepair,type='rowmajor')
             !lpstate=spacefilling_coord2index(this%nstate,statepair,type='hilbert')

             !calculate L-term
             !default value
             Lterm=(0.0_double,0.0_double)
             
             !i=-2
             if(i.EQ.-2)then
                !j=0
                if(j.EQ.0)then
                   Lterm=sqrt(real((istate+1)*(istate+2),double))*0.5_double*&
                        (eye*(omegaA2-gammaRpq)&
                        +(omegaB2-gammaGamma))
                   !Lterm=sqrt(real((istate+1)*(istate+2),double))*0.5_double
                   !LRe=omegaB2-gammaGamma
                   !LIm=omegaA2-gammaRpq
                   !Lterm=Lterm*(LRe+eye*LIm) !inexplicably commented out
                end if
             end if
             
             !i=-1
             if(i.EQ.-1)then
                !j=-1,1
                if(j.EQ.-1)then
                   Lterm=sqrt(real((istate+1)*(jstate+1),double))&
                        *(omegaB1+gammaGamma)
                end if
                if(j.EQ.1)then
                   Lterm=sqrt(real((istate+1)*jstate,double))&
                        *(eye*gammaRpq-omegaB2)
                end if
             end if
             
             !i=0
             if(i.EQ.0)then
                !j=-2,0,2
                if(j.EQ.-2)then
                   Lterm=-sqrt(real((jstate+1)*(jstate+2),double))*0.5_double*&
                        (eye*(omegaA2-gammaRpq)&
                        -(omegaB2-gammaGamma))
                end if
                if(j.EQ.0)then
                   Lterm=(-eye*modulo(istate-jstate,this%nstate)*omegaA1&
                        -(istate+jstate+1)*omegaB1+gammaGamma)
                   !LRe=-(istate+jstate+1)*omegaB1+gammaGamma
                   !LIm=-modulo(istate-jstate,this%nstate)*omegaA1
                   !Lterm=LRe!+eye*LIm !inexplicably commented out
                end if
                if(j.EQ.2)then
                   Lterm=-sqrt(real(jstate*(jstate-1),double))*0.5_double*&
                        (eye*(omegaA2+gammaRpq)&
                        -(omegaB2+gammaGamma))
                end if
             end if
             
             !i=1
             if(i.EQ.1)then
                !j=-1,1
                if(j.EQ.-1)then
                   Lterm=-sqrt(real(istate*(jstate+1),double))&
                        *(eye*gammaRpq+omegaB2)
                end if
                if(j.EQ.1)then
                   Lterm=sqrt(real(istate*jstate,double))&
                        *(omegaB1-gammaGamma)
                end if
             end if
             
             !i=2
             if(i.EQ.2)then
                !j=0
                if(j.EQ.0)then
                   Lterm=sqrt(real(jstate*(jstate-1),double))*0.5_double*&
                        (eye*(omegaA2+gammaRpq)&
                        +(omegaB2+gammaGamma))
                end if
             end if
                          
             !accumulate supermatrixL term
             this%supermatrixL(lstate,lpstate)&
                  =this%supermatrixL(lstate,lpstate)+Lterm
             
          end do
       end do

    end do !lstate loop
    !multiply supermatrixL by units of omega_c
    this%supermatrixL=this%supermatrixL*this%omega_c


  contains
    real(double) function G0(beta,t)
      real(double),intent(in)::beta,t
      integer(long),parameter::MatsubaraNmin=10 !minimum matsubara number

      integer(long)::n
      real(double)::twopioverbeta,D,D2,term

      !constants
      twopioverbeta=twopi/beta

      !initial value
      G0=0.0_double

      !accumulate sum over minimum Matsubara numbers
      D=0.0_double
      do n=1,MatsubaraNmin
         D=D+twopioverbeta
         D2=D*D
         term=F0(D*t)/(D2-1)
         G0=G0+term
         !if(t>1.0) write(444,*)G0
      end do
      !accumulate additional matusbara terms while
      ! their relative contribution is greater than 0.05%  
      do while(term/G0.GT.5E-4)
         D=D+twopioverbeta
         D2=D*D
         term=F0(D*t)/(D2-1)
         G0=G0+term
         !if(t>1.0) write(444,*)G0
      end do
      !multiply by units
      G0=G0*4.0_double/beta
      !add high temperature term
      G0=G0+F0(t)/tan(0.5_double*beta)
      !if(t>1.0)flush(444)
      !if(t>1.0)stop
    end function G0 !verified
    !-----------------------
    real(double) function G1tilde(beta,t)
      real(double),intent(in)::beta,t
      integer(long),parameter::MatsubaraNmin=10 !minimum matsubara number 

      integer(long)::n
      real(double)::twopioverbeta,D,D2,term 

      !constants
      twopioverbeta=twopi/beta

      !initial value
      G1tilde=0.0_double

      !accumulate sum over minimum Matsubara numbers
      D=0.0_double
      do n=1,MatsubaraNmin
         D=D+twopioverbeta
         D2=D*D
         term=F1tilde(D*t)/(D*(D2-1))
         G1tilde=G1tilde+term
         !if(t>1.0) write(444,*)G1tilde
      end do
      !accumulate additional matusbara terms while
      ! their relative contribution is greater than 0.05%  
      do while(term/G1tilde.GT.5E-4)
         !n=n+1
         D=D+twopioverbeta
         D2=D*D
         term=F1tilde(D*t)/(D*(D2-1))
         G1tilde=G1tilde+term
         !if(t>1.0) write(444,*)G1tilde
      end do
      !multiply by units
      G1tilde=G1tilde*4.0_double/beta
      !add high temperature term
      G1tilde=G1tilde+F1tilde(t)/tan(0.5_double*beta)
      !if(t>1.0)flush(444)
      !if(t>1.0)stop

    end function G1tilde !verified
    !-----------------------
    real(double) function G2(beta,t)
      real(double),intent(in)::beta,t
      integer(long),parameter::MatsubaraNmin=10 !minimum matsubara number

      integer(long)::n
      real(double)::twopioverbeta,D,D2,term

      !constants
      twopioverbeta=twopi/beta

      !initial value
      G2=0.0_double

      !accumulate sum over minimum Matsubara numbers
      D=0.0_double
      do n=1,MatsubaraNmin
         D=D+twopioverbeta
         D2=D*D
         term=F2(D*t)/(D2*(D2-1))
         G2=G2+term
         !if(t>1.0) write(444,*)G2
      end do
      !accumulate additional matusbara terms while
      ! their relative contribution is greater than 0.05%  
      do while(term/G2.GT.5E-4)
         D=D+twopioverbeta
         D2=D*D
         term=F2(D*t)/(D2*(D2-1))
         G2=G2+term
         !if(t>1.0) write(444,*)G2
      end do
      !add units
      G2=G2*4.0_double/beta
      !add high temperature term
      G2=G2+F2(t)/tan(0.5_double*beta)
      !if(t>1.0)flush(444)
      !if(t>1.0)stop

    end function G2 !verified
    !-----------------------
    real(double) function F0(t)
      real(double),intent(in)::t
      !F0=1.0_double !initial value
      F0=1.0_double-exp(-t)
    end function F0 !verified
    !-----------------------
    real(double) function F1tilde(t)
      real(double),intent(in)::t
      !F1tilde=0.0_double !initial value
      F1tilde=1.0_double-(1.0_double+t-0.5_double*t*t)*exp(-t)
    end function F1tilde
    !-----------------------
    real(double) function F2(t)
      real(double),intent(in)::t
      !F2=0.0_double !initial value
      F2=1.0_double-(1.0_double+t+0.5_double*t*t)*exp(-t)
    end function F2
    !-----------------------
  end subroutine GQFPE_update

  !======================================================================
  !> \brief Re-initiallizes the GQFPE object.
  !> \param THIS is the GQFPE  object to be re-initialized.
  !> \param STATE is an optional integer:
  !>        when 0, will create a null state by deallocating all dynamic
  !>        memory and returning the object to an un-initiallized state;
  !>        when not 0, will return the object to the default settings;
  !>        when not present, object will reset based on current scalar
  !>        parameters.
  !======================================================================
  subroutine GQFPE_reset(this,state)
    type(GQFPE),intent(inout)::this
    integer(long),intent(in),optional::STATE
    integer(long)::statepair(2),i,j

    if(present(state))then
       if(state.EQ.0)then
          !nullify all pointers
          !******        Example - cleanup pointer attribute 'PPP'       ****
          !if(associated(this%PPP))nullify(this%PPP)
          !******************************************************************
          if(associated(this%supermatrixL))nullify(this%supermatrixL)
          if(associated(this%expancoef))nullify(this%expancoef)
          if(associated(this%expancoefdot))nullify(this%expancoefdot)
          if(associated(this%density))nullify(this%density)
          
          !kill all sub-objects
          !**** example **********
          !call kill(this%object)
          !***********************
          
          !set all scalar parameters to error values
          !**** example **********
          !this%nstate=-1
          !***********************
          this%dt=0.0_double
          this%time=Huge(this%time)
          this%R_pq=huge(this%R_pq)
          this%omega_c=0.0_double
          this%beta=-1.0_double
          this%R_qq=huge(this%R_qq)
          this%R_pp=huge(this%R_qq)
          this%D=huge(this%D)
          this%nstate=0

          !un-initialize object
          this%initialized=.false.
       else
          !allocate dynamic memory
          !***  Example - cleanup pointer attribute 'PPP'     ***
          !***            then reallocate memory              ***
          !if(associated(this%PPP))nullify(this%PPP)
          !allocate(this%PPP(0:this%npt-1))
          !******************************************************
          if(associated(this%supermatrixL))nullify(this%supermatrixL)
          allocate(this%supermatrixL(0:this%nstate**2-1,0:this%nstate**2-1))

          if(associated(this%expancoef))nullify(this%expancoef)
          allocate(this%expancoef(0:this%nstate**2-1))
          
          if(associated(this%expancoefdot))nullify(this%expancoefdot)
          allocate(this%expancoefdot(0:this%nstate**2-1))
          
          if(associated(this%density))nullify(this%density)
          allocate(this%density(0:this%nstate-1,0:this%nstate-1))
          
          !Set default dynamic memory values
          !***  Example - set values in pointer 'PPP' to zero ***
          !this%PPP(:)=0.0_double
          !******************************************************
          !default density is superposition of first two states
          this%density=(0.0_double,0.0_double)
          do i=0,1
             do j=0,1
                this%density(i,j)=(0.5_double,0.0_double)
             end do
          end do
          !overwrite sub-object default static parameters
          !***       example      ***
          !***set object static parameter***
          !this%object%XXX=123
          !**************************

          !reset all sub objects to correct any memory issues
          !***      example     ***
          !call reset(this%object,state=1)
          !************************

          !overwrite sub-object default dynamic pointer arrays
          !***       example      ***
          !***set object pointer array values***
          !this%object%PPP(:)=XXX
          !**************************

          !declare initialization complete
          this%initialized=.true.
          
       end if
    end if
    
    !reset object based on current static parameters
    if(this%initialized)then

       !Sample calculated variable initial values
       !***  Example - attribute 'var' samples a Gaussian random number
       ! this%var=gran()

       !update expansion coefs according to current density
       this%expancoef=rowmajorcurve(this%density)

       !set expansion coefs rate to zero
       this%expancoefdot=(0.0_double,0.0_double)

       !set supermatrixL to zero
       this%supermatrixL=(0.0_double,0.0_double)

       !Resample sub-objects
       !**** example **********
       !call reset(this%object)
       !***********************

       !update now that object is fully reset and not in null state
       call update(this)

    end if
    

  end subroutine GQFPE_reset

  !======================================================================
  !> \brief Backups the current state of the GQFPE object to file.
  !> \param[in] THIS is the GQFPE  object to be updated.
  !> \param[in] FILE is a string containing the location of the backup file.
  !======================================================================
  subroutine GQFPE_backup(this,file)
    use filemanager
    use string
    use testing_class
    type(GQFPE),intent(in)::this
    character*(*),intent(in)::file
    integer(short)::unit
    logical::fileisopen
    integer(long)::i,j,k,l
    
    !check input file
    inquire(file=file,opened=fileisopen,number=unit)
    if(unit.LT.0)unit=newunit()
    if(.not.fileisopen)open(unit,file=file)
    
    !check GQFPE object
    call assert(check(this).EQ.0,msg='GQFPE object does not pass check.')

    !always write the data type on the first line
    write(unit,*)'GQFPE'

    !******      Backup below all the derived type's attributes       ****
    !******         in the order the MAKE method reads them           ****


    !First, Scalar attributes
    !******          Example - Backup a scalar attribute            ******
    ! write(unit,*)this%var
    !*********************************************************************
    write(unit,*)this%dt
    write(unit,*)this%time
    write(unit,*)this%omega_c
    write(unit,*)this%beta
    write(unit,*)this%nstate
 

    !Second, Dynamic parameter arrays
    !***       Example - Backup an NxM matrix attribute                ***
    ! write(unit,*)((this%matrix(i,j),j=1,M),i=1,N)
    !*********************************************************************
    write(unit,*)((this%density(i,j)&
         ,j=0,size(this%density,2)-1)&
         ,i=0,size(this%density,1)-1)
    
    !Last,objects
    !******              Example - Backup an object            ***********
    ! call backup(this%object,file//'.object')
    ! write(unit,*)quote(file//'.object')!write the object location
    !*********************************************************************
    

    !finished writing all attributes - now close backup file
    close(unit)  
  end subroutine GQFPE_backup
  
  !======================================================================
  !> \brief Retrun the current state of GQFPE as a string.
  !> \param[in] THIS is the GQFPE object.
  !> \param[in] MSG is an optional string message to annotate the status.
  !======================================================================
  character(len=line) function GQFPE_status(this,msg)
    type(GQFPE),intent(in)::this
    character*(*),intent(in),optional::msg
    character(len=5)::FMT='(A)'
    
    !Edit the status prompt to suit your needs
    write(GQFPE_status,FMT)'GQFPE status is currently not available'
    
  end function GQFPE_status

 !======================================================================
  !> \brief Checks the GQFPE object.
  !> \param[in] THIS is the GQFPE object to be checked.
  !> \return 0 if all checks pass or exit at first failed check and returm non zero.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function GQFPE_check(this)
    use testing_class
    type(GQFPE),intent(in)::this

    !initiate with no problems found 
    GQFPE_check=0

    !check that object is initialized
    call assert(this%initialized&
         ,msg='GQFPE_check: GQFPE object not initialized.'&
         ,iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check that object has correct name
    call assert(this%name.EQ.'GQFPE'&
         ,msg='GQFPE_check: GQFPE name is not set.'&
         ,iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check dt is not tiny
    call assert(abs(this%dt).GT.epsilon(this%dt)&
         ,msg='GQFPE_check: dt is tiny',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check dt is well behaved
    call assert(check(this%dt).EQ.0&
         ,msg='GQFPE_check: dt failed check',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check time is well behaved
    call assert(check(this%time).EQ.0&
         ,msg='GQFPE_check: time failed check',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check calculated variable R_pq is well behaved
    call assert(check(this%R_pq).EQ.0&
         ,msg='GQFPE_check: R_pq failed check',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check calculated variable omega_c is well behaved
    call assert(check(this%omega_c).EQ.0&
         ,msg='GQFPE_check: omega_c failed check',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check omega_c is not negative
    call assert(this%omega_c.GT.epsilon(this%omega_c)&
         ,msg='GQFPE_check: omega_c is not positive',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check beta is well behaved
    call assert(check(this%beta).EQ.0&
         ,msg='GQFPE_check: beta failed check',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check beta is not negative
    call assert(this%beta.GT.epsilon(this%beta)&
         ,msg='GQFPE_check: beta is not positive',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check calculated variable R_qq is well behaved
    call assert(check(this%R_qq).EQ.0&
         ,msg='GQFPE_check: R_qq failed check',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check calculated variable R_pp is well behaved
    call assert(check(this%R_pp).EQ.0&
         ,msg='GQFPE_check: R_pp failed check',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check calculated variable D is well behaved
    call assert(check(this%D).EQ.0&
         ,msg='GQFPE_check: D failed check',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check nstate is well behaved
    call assert(check(this%nstate).EQ.0&
         ,msg='GQFPE_check: nstate failed check',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return
    
    !check nstate is positive
    call assert(this%nstate.GT.0&
         ,msg='GQFPE_check: nstate is not a positive integer'&
         ,iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return
    

    !check supermatrixL is of correct size
    call assert(size(this%supermatrixL,1).EQ.this%nstate**2&
         ,msg='GQFPE check: dynamic pointer array supermatrixL is not of size&
         & nstate**2 for dimension 1',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return
    call assert(size(this%supermatrixL,2).EQ.this%nstate**2&
         ,msg='GQFPE check: dynamic pointer array supermatrixL is not of size&
         & nstate**2 for dimension 2',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check supermatrixL is well behaved
    call assert(check(this%supermatrixL).EQ.0&
         ,msg='GQFPE_check: supermatrixL failed check',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check expancoef is of correct size
    call assert(size(this%expancoef,1).EQ.this%nstate**2&
         ,msg='GQFPE check: dynamic pointer array expancoef is not of size&
         & nstate**2 for dimension 1',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check expancoef is well behaved
    call assert(check(this%expancoef).EQ.0&
         ,msg='GQFPE_check: expancoef failed check',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check expancoefdot is of correct size
    call assert(size(this%expancoefdot,1).EQ.this%nstate**2&
         ,msg='GQFPE check: dynamic pointer array expancoefdot is not of size&
         & nstate**2 for dimension 1',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check expancoefdot is well behaved
    call assert(check(this%expancoefdot).EQ.0&
         ,msg='GQFPE_check: expancoefdot failed check',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check density is of correct size
    call assert(size(this%density,1).EQ.this%nstate&
         ,msg='GQFPE check: dynamic pointer array density is not of size&
         & nstate for dimension 1',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return
    call assert(size(this%density,2).EQ.this%nstate&
         ,msg='GQFPE check: dynamic pointer array density is not of size&
         & nstate for dimension 2',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return

    !check density is well behaved
    call assert(check(this%density).EQ.0&
         ,msg='GQFPE_check: density failed check',iostat=GQFPE_check)
    if(GQFPE_check.NE.0)return


    !Check all attributes are within acceptable values
    !**********   Example - check an object attribute 'that'  *********
    !call assert(check(this%that).EQ.0&
    !     ,msg='GQFPE_check: that sub-object failed check'&
    !     ,iostat=GQFPE_check)
    !if(GQFPE_check.NE.0)return
    !***********************************************************************
    
    !***   Example - check an integer attribute 'ndim' is well behaved   ***
    !call assert(check(this%ndim).EQ.0&
    !     ,msg='GQFPE_check: ndim failed check',iostat=GQFPE_check)
    !if(GQFPE_check.NE.0)return
    !***********************************************************************
 
    !*** Example - add a constrain that says 'ndim' can only be positive ***
    !call assert(this%ndim.GT.0&
    !     ,msg='GQFPE_check: ndim is not positive',iostat=GQFPE_check)
    !if(GQFPE_check.NE.0)return
    !***********************************************************************

    !***  Example - check a real valued attribute 'var' is well behaved  ***
    !call assert(check(this%var).EQ.0&
    !     ,msg='GQFPE_check: var failed check',iostat=GQFPE_check)
    !if(GQFPE_check.NE.0)return
    !***********************************************************************

    !***  Example - add a constrain that says 'var' can not be zero     ***
    !call assert(abs(this%var).GT.epsilon(this%var)&
    !     ,msg='GQFPE_check: var is tiny',iostat=GQFPE_check)
    !if(GQFPE_check.NE.0)return
    !***********************************************************************

    !***  Example - check a real valued pointer attribute 'matrix'       ***
    !***            is well behaved                                      ***
    !call assert(check(this%matrix).EQ.0&
    !     ,msg='GQFPE_check: matrix failed check',iostat=GQFPE_check)
    !if(GQFPE_check.NE.0)return
    !***********************************************************************

    !********* Example - check an NxM matrix has right dimensions **********
    !call assert(size(this%matrix).EQ.N*M&
    !     ,msg='GQFPE_check: number of matrix elements not = N*M.'&
    !     ,iostat=GQFPE_check)
    !if(GQFPE_check.NE.0)return
    !***********************************************************************

  end function GQFPE_check
  !-----------------------------------------
  !======================================================================
  !> \brief Tests the GQFPE methods.
  !> \param[in] this is the GQFPE object whose methods will be excercised.
  !> \remark Will stop after first failed test.
  !======================================================================
  subroutine GQFPE_test
    use testing_class
    use filemanager
    type(GQFPE)::this
    character(len=label)::string
    integer(long)::unit
    integer(short)::ierr

    !diagnostic variables
    integer(long)::nstep
    real(double)::prevtime
    !arbitrary gussian density
    integer(long),parameter::nstate=100,ngauss=3,nparam=2,npt=1023,ndump=4 
    real(double),parameter::qmax=20.0_double,qmin=-qmax
    real(double)::GaussParam(ngauss,nparam) !mu=1 and sigma=2 for nparam
    real(double)::GaussWeight(ngauss)
    real(double)::psi(0:npt-1),phi(0:npt-1),grid(0:npt-1)
    real(double)::overlap(0:nstate-1),chi(0:nstate-1,0:npt-1)
    integer(long)::basismap(0:nstate-1)
    real(double)::q,dq
    integer(long)::ipt,n,i,j,statepair(2),k,kprime
    integer(long)::istate,jstate,lstate

    !verify GQFPE is compatible with current version
    include 'verification'
   
    write(*,*)'test static parameter dt is stored properly in backup file'
    call make(this) !make GQFPE
    this%dt=ps    !manually set static parameter to non-default value
    call reset(this,state=1)
    call system('rm -f GQFPE.tmpfile*') !remove any previous backup file(s)
    call backup(this,file='GQFPE.tmpfile') !create backup file
    call kill(this) !destroy GQFPE
    call make(this,file='GQFPE.tmpfile') !make new GQFPE from backup file
    !assert non default parameter is conserved
    call assert(this%dt.EQ.ps,&
         msg='GQFPE static parameter dt is not stored properly')
    call kill(this)    !destroy GQFPE to clean up memory
    call system('rm -f GQFPE.tmpfile*') !remove backup file(s)

    write(*,*)'test make sets correct default value for static parameter dt'
    call make(this) !make GQFPE
    call assert(this%dt.EQ.fs,&
         msg='GQFPE default static parameter dt is not fs')
    call kill(this)

    write(*,*)'test edge case for static parameter dt breaks GQFPE'
    call make(this) !make GQFPE
    this%dt=0.0_double
    call assert(check(this).NE.0,&
         msg='edge case value 0.0 for static parameter dt &
         &does not break GQFPE')
    call kill(this)

    write(*,*)'test edge case for static parameter dt breaks GQFPE'
    call make(this) !make GQFPE
    this%dt=huge(this%dt)
    call assert(check(this).NE.0,&
         msg='edge case huge value for static parameter dt &
         &does not break GQFPE')
    call kill(this)

    write(*,*)'test static parameter time is stored properly in backup file'
    call make(this) !make GQFPE
    this%time=ps    !manually set static parameter to non-default value
    call reset(this,state=1)
    call system('rm -f GQFPE.tmpfile*') !remove any previous backup file(s)
    call backup(this,file='GQFPE.tmpfile') !create backup file
    call kill(this) !destroy GQFPE
    call make(this,file='GQFPE.tmpfile') !make new GQFPE from backup file
    !assert non default parameter is conserved
    call assert(this%time.EQ.ps,&
         msg='GQFPE static parameter time is not stored properly')
    call kill(this)    !destroy GQFPE to clean up memory
    call system('rm -f GQFPE.tmpfile*') !remove backup file(s)

    write(*,*)'test make sets correct default value for static parameter time'
    call make(this) !make GQFPE
    call assert(this%time.EQ.0.0_double,&
         msg='GQFPE default static parameter time is not 0.0')
    call kill(this)

    write(*,*)'test edge case for static parameter time breaks GQFPE'
    call make(this) !make GQFPE
    this%time=huge(this%time)
    call assert(check(this).NE.0,&
         msg='edge case huge value for static parameter time &
         &does not break GQFPE')
    call kill(this)

    write(*,*)'test make sets default value for calculated variable R_pq'
    call make(this) !make GQFPE
    call assert(this%R_pq.EQ.0.0_double,&
         msg='GQFPE default calculated variable R_pq is not 0.0_double')
    call kill(this)

    write(*,*)'test edge case for calculated variable R_pq breaks GQFPE'
    call make(this) !make GQFPE
    this%R_pq=huge(this%R_pq)
    call assert(check(this).NE.0,&
         msg='edge case huge value for calculated variable R_pq &
         &does not break GQFPE')
    call kill(this)

    write(*,*)'check step method increments time'
    call make(this)
    prevtime=this%time
    call step(this)
    call assert(this%time.GT.prevtime&
         ,msg='step method does not increment time')
    call kill(this)

    write(*,*)'test static parameter omega_c is stored properly in backup file'
    call make(this) !make GQFPE
    this%omega_c=invcm !manually set static parameter to non-default value
    call reset(this,state=1)
    call system('rm -f GQFPE.tmpfile*') !remove any previous backup file(s)
    call backup(this,file='GQFPE.tmpfile') !create backup file
    call kill(this) !destroy GQFPE
    call make(this,file='GQFPE.tmpfile') !make new GQFPE from backup file
    !assert non default parameter is conserved
    call assert(this%omega_c.EQ.invcm,&
         msg='GQFPE static parameter omega_c is not stored properly')
    call kill(this)    !destroy GQFPE to clean up memory
    call system('rm -f GQFPE.tmpfile*') !remove backup file(s)

    write(*,*)'test make sets correct default value for static parameter omega_c'
    call make(this) !make GQFPE
    call assert(this%omega_c.EQ.200._double*invcm,&
         msg='GQFPE default static parameter omega_c is not 200 wavenumbers')
    call kill(this)

    write(*,*)'test edge case for static parameter omega_c breaks GQFPE'
    call make(this) !make GQFPE
    this%omega_c=0.0_double
    call assert(check(this).NE.0,&
         msg='edge case value 0.0 for static parameter omega_c &
         &does not break GQFPE')
    call kill(this)

    write(*,*)'test edge case for static parameter omega_c breaks GQFPE'
    call make(this) !make GQFPE
    this%omega_c=huge(this%omega_c)
    call assert(check(this).NE.0,&
         msg='edge case huge value for static parameter omega_c &
         &does not break GQFPE')
    call kill(this)

    write(*,*)'test edge case for static parameter omega_c breaks GQFPE'
    call make(this) !make GQFPE
    this%omega_c=-1.0_double
    call assert(check(this).NE.0,&
         msg='edge case negative value for static parameter omega_c &
         &does not break GQFPE')
    call kill(this)

    write(*,*)'test static parameter beta is stored properly in backup file'
    call make(this) !make GQFPE
    !manually set static parameter to non-default value
    this%beta=1.0_double/(500*Kel*kB)
    call reset(this,state=1)
    call system('rm -f GQFPE.tmpfile*') !remove any previous backup file(s)
    call backup(this,file='GQFPE.tmpfile') !create backup file
    call kill(this) !destroy GQFPE
    call make(this,file='GQFPE.tmpfile') !make new GQFPE from backup file
    !assert non default parameter is conserved
    call assert(this%beta.EQ.1.0_double/(500*Kel*kB),&
         msg='GQFPE static parameter beta is not stored properly')
    call kill(this)    !destroy GQFPE to clean up memory
    call system('rm -f GQFPE.tmpfile*') !remove backup file(s)

    write(*,*)'test make sets correct default value for static parameter beta'
    call make(this) !make GQFPE
    call assert(this%beta.EQ.5.0_double/(this%omega_c),&
         msg='GQFPE default static parameter beta is not 5/(hbar*omega_c)')
    call kill(this)

    write(*,*)'test edge case for static parameter beta breaks GQFPE'
    call make(this) !make GQFPE
    this%beta=huge(this%beta)
    call assert(check(this).NE.0,&
         msg='edge case huge value for static parameter beta &
         &does not break GQFPE')
    call kill(this)

    write(*,*)'test edge case for static parameter beta breaks GQFPE'
    call make(this) !make GQFPE
    this%beta=-1.0_double
    call assert(check(this).NE.0,&
         msg='edge case negative value for static parameter beta &
         &does not break GQFPE')
    call kill(this)

    write(*,*)'test 1.5<R_pq<2.0 at t_s=10 for beta_s=1.0'
    call make(this)
    this%beta=1.0_double/(hbar*this%omega_c)
    nstep=0
    do while(this%time*this%omega_c<10.0_double.and.nstep.lt.10000)
       call step(this) !perform one evolution step
       nstep=nstep+1
       !write(123,*)this%time*this%omega_c,this%R_pq
    end do
    call assert(this%R_pq,1.75_double,0.25_double&
         ,msg='Value of R_pq is not between 1.5 and 2.0 at t_s=10')
    call kill(this)

    write(*,*)'test 3.5<R_pq<4.0 at t_s=10 for beta_s=0.5'
    call make(this)
    this%beta=0.5_double/(hbar*this%omega_c)
    nstep=0
    do while(this%time*this%omega_c<10.0_double.and.nstep.lt.10000)
       call step(this) !perform one evolution step
       nstep=nstep+1
       !write(123,*)this%time*this%omega_c,this%R_pq
    end do
    call assert(this%R_pq,3.75_double,0.25_double&
         ,msg='Value of R_pq is not between 3.5 and 4.0 at t_s=10')
    call kill(this)

    write(*,*)'test R_pq is negative at t_s=10 when beta_s=5'
    call make(this)
    this%beta=5.0_double/(hbar*this%omega_c)
    nstep=0
    do while(this%time*this%omega_c<10.0_double.and.nstep.lt.10000)
       call step(this) !perform one evolution step
       nstep=nstep+1
       !write(123,*)this%time*this%omega_c,this%R_pq
    end do
    call assert(this%R_pq.LT.0.0&
         ,msg='Value of R_pq is negative at t_s=10')
    call kill(this)

    !----- Rqq
    write(*,*)'test make sets default value for calculated variable R_qq'
    call make(this) !make GQFPE
    call assert(this%R_qq.EQ.0.0_double,&
         msg='GQFPE default calculated variable R_qq is not 0.0_double')
    call kill(this)

    write(*,*)'test edge case for calculated variable R_qq breaks GQFPE'
    call make(this) !make GQFPE
    this%R_qq=huge(this%R_qq)
    call assert(check(this).NE.0,&
         msg='edge case huge value for calculated variable R_qq &
         &does not break GQFPE')
    call kill(this)

    write(*,*)'test 0.125<R_qq<0.150 at t_s=10 for beta_s=1.0'
    call make(this)
    this%beta=1.0_double/(hbar*this%omega_c)
    nstep=0
    do while(this%time*this%omega_c<10.0_double.and.nstep.lt.10000)
       call step(this) !perform one evolution step
       nstep=nstep+1
       !write(123,*)this%time*this%omega_c,this%R_qq
    end do
    call assert(this%R_qq,0.1375_double,0.0125_double&
         ,msg='Value of R_qq is not between 0.125 and 0.150 at t_s=10')
    call kill(this)

    write(*,*)'test R_qq=0.3+/-0.01 at t_s=10 for beta_s=0.5'
    call make(this)
    this%beta=0.5_double/(hbar*this%omega_c)
    nstep=0
    do while(this%time*this%omega_c<10.0_double.and.nstep.lt.10000)
       call step(this) !perform one evolution step
       nstep=nstep+1
       !write(123,*)this%time*this%omega_c,this%R_qq
    end do
    call assert(this%R_qq,0.3_double,0.01_double&
         ,msg='Value of R_qq is not 3.0+/-0.125 at t_s=10')
    call kill(this)

    write(*,*)'test R_qq is negative at t_s=10 for beta_s=5'
    call make(this)
    this%beta=5_double/(hbar*this%omega_c)
    nstep=0
    do while(this%time*this%omega_c<10.0_double.and.nstep.lt.10000)
       call step(this) !perform one evolution step
       nstep=nstep+1
       !write(123,*)this%time*this%omega_c,this%R_qq
    end do
    call assert(this%R_qq.LT.0.0_double&
         ,msg='Value of R_qq is not negative at t_s=10')
    call kill(this)


    !----- Rpp
    write(*,*)'test make sets default value for calculated variable R_pp'
    call make(this) !make GQFPE
    call assert(this%R_pp.EQ.0.0_double,&
         msg='GQFPE default calculated variable R_pp is not 0.0_double')
    call kill(this)

    write(*,*)'test edge case for calculated variable R_pp breaks GQFPE'
    call make(this) !make GQFPE
    this%R_pp=huge(this%R_pp)
    call assert(check(this).NE.0,&
         msg='edge case huge value for calculated variable R_pp &
         &does not break GQFPE')
    call kill(this)

    write(*,*)'test R_pp=15.0+/-2.5 at t_s=10 for beta_s=1.0'
    call make(this)
    this%beta=1.0_double/(hbar*this%omega_c)
    nstep=0
    do while(this%time*this%omega_c<10.0_double.and.nstep.lt.10000)
       call step(this) !perform one evolution step
       nstep=nstep+1
       !write(123,*)this%time*this%omega_c,this%R_pp
    end do
    call assert(this%R_pp,15.0_double,2.5_double&
         ,msg='Value of R_pp is not 15.0 +/- 2.5 at t_s=10')
    call kill(this)

    write(*,*)'test 30.0<R_pp<32.5 t_s=10 for beta_s=0.5'
    call make(this)
    this%beta=0.5_double/(hbar*this%omega_c)
    nstep=0
    do while(this%time*this%omega_c<10.0_double.and.nstep.lt.10000)
       call step(this) !perform one evolution step
       nstep=nstep+1
       !write(123,*)this%time*this%omega_c,this%R_pp
    end do
    call assert(this%R_pp,31.25_double,1.25_double&
         ,msg='Value of R_pp is not between 30 and 32.5 at t_s=10')
    call kill(this)

    write(*,*)'test R_pp is in range [4:5] t_s=10 for beta_s=5'
    call make(this)
    this%beta=5_double/(hbar*this%omega_c)
    nstep=0
    do while(this%time*this%omega_c<10.0_double.and.nstep.lt.10000)
       call step(this) !perform one evolution step
       nstep=nstep+1
       !write(123,*)this%time*this%omega_c,this%R_pp
    end do
    call assert(this%R_pp,4.5_double,0.5_double&
         ,msg='Value of R_pp is not between [4:5] at t_s=10')
    call kill(this)

    !----- D
    write(*,*)'test make sets default value for calculated variable D'
    call make(this) !make GQFPE
    call assert(this%D.EQ.0.0_double,&
         msg='GQFPE default calculated variable D is not 0.0_double')
    call kill(this)

    write(*,*)'test edge case for calculated variable D breaks GQFPE'
    call make(this) !make GQFPE
    this%D=huge(this%D)
    call assert(check(this).NE.0,&
         msg='edge case huge value for calculated variable D &
         &does not break GQFPE')
    call kill(this)

    write(*,*)'test D=5.0+/-2.5 at t_s=10 for beta_s=1.0'
    call make(this)
    this%beta=1.0_double/(hbar*this%omega_c)
    nstep=0
    do while(this%time*this%omega_c<10.0_double.and.nstep.lt.10000)
       call step(this) !perform one evolution step
       nstep=nstep+1
       write(123,*)this%time*this%omega_c,this%D
    end do
    call assert(this%D,5.0_double,2.5_double&
         ,msg='Value of D is not 5.0 +/- 2.5 at t_s=10')
    call kill(this)

    write(*,*)'test D=22.5 +/- 2.5 t_s=10 for beta_s=0.5'
    call make(this)
    this%beta=0.5_double/(hbar*this%omega_c)
    nstep=0
    do while(this%time*this%omega_c<10.0_double.and.nstep.lt.10000)
       call step(this) !perform one evolution step
       nstep=nstep+1
       write(123,*)this%time*this%omega_c,this%D
    end do
    call assert(this%D,22.5_double,2.5_double&
         ,msg='Value of D is not between 22.5 +/- 2.5 at t_s=10')
    call kill(this)

    write(*,*)'test D=0.0 +/- 2.5 at t_s=10 for beta_s=5'
    call make(this)
    this%beta=5_double/(hbar*this%omega_c)
    nstep=0
    do while(this%time*this%omega_c<10.0_double.and.nstep.lt.10000)
       call step(this) !perform one evolution step
       nstep=nstep+1
       write(123,*)this%time*this%omega_c,this%D
    end do
    call assert(this%D,0.0_double,2.5_double&
         ,msg='Value of D is not 0.0 +/- 2.5 at t_s=10')
    call kill(this)

    !write(*,*)'Diagnostic: plot Beta_s parameterizd by temperature and omega_c'
    !open(123,file='beta_omega.dat')
    !do i=1,400
    !   do j=1,2000
    !      q=hbar*j*invcm/(kb*i)
    !      write(123,*)i,j,q
    !   end do
    !   write(123,*)
    !end do
    !close(123)

    write(*,*)'Diagnostic: writing data for Fig1 in Jang Manuscript.'
    open(123,file='GQFPE_Fig1_beta_s0p5.dat')
    write(123,*)'# Jang Manuscript Figure 1 data for beta_s=0.5'
    write(123,*)'# time_s, R_pq, R_qq, R_pp, D'
    call make(this)
    this%beta=0.5_double/(hbar*this%omega_c)
    nstep=0
    do while(this%time*this%omega_c<10.0_double.and.nstep.lt.10000)
       call step(this) !perform one evolution step
       nstep=nstep+1
       write(123,*)this%time*this%omega_c,this%R_pq,this%R_qq,this%R_pp,this%D
    end do
    close(123)
    write(*,*)'beta_s=',this%beta*hbar*this%omega_c
    write(*,*)'gamma_s = 0.1 (defined)'
    write(*,*)'t_s=',this%time*this%omega_c
    write(*,*)'dt_s=',this%dt*this%omega_c
    write(*,*)'omega_c=',this%omega_c
    write(*,*)'temperature=',1.0_double/(kB*this%beta)
    call kill(this)
    open(123,file='GQFPE_Fig1_beta_s1p0.dat')
    write(123,*)'# Jang Manuscript Figure 1 data for beta_s=1.0'
    write(123,*)'# time_s, R_pq, R_qq, R_pp, D'
    call make(this)
    this%beta=1.0_double/(hbar*this%omega_c)
    nstep=0
    do while(this%time*this%omega_c<10.0_double.and.nstep.lt.10000)
       call step(this) !perform one evolution step
       nstep=nstep+1
       write(123,*)this%time*this%omega_c,this%R_pq,this%R_qq,this%R_pp,this%D
    end do
    close(123)
    call kill(this)
    open(123,file='GQFPE_Fig1_beta_s5p0.dat')
    write(123,*)'# Jang Manuscript Figure 1 data for beta_s=5.0'
    write(123,*)'# time_s, R_pq, R_qq, R_pp, D'
    call make(this)
    this%beta=5.0_double/(hbar*this%omega_c)
    nstep=0
    do while(this%time*this%omega_c<10.0_double.and.nstep.lt.10000)
       call step(this) !perform one evolution step
       nstep=nstep+1
       write(123,*)this%time*this%omega_c,this%R_pq,this%R_qq,this%R_pp,this%D
    end do
    close(123)
    call kill(this)


    write(*,*)'test static parameter nstate is stored properly in backup file'
    call make(this) !make GQFPE
    this%nstate=10    !manually set static parameter to non-default value
    call reset(this,state=1)
    call system('rm -f GQFPE.tmpfile*') !remove any previous backup file(s)
    call backup(this,file='GQFPE.tmpfile') !create backup file
    call kill(this) !destroy GQFPE
    call make(this,file='GQFPE.tmpfile') !make new GQFPE from backup file
    !assert non default parameter is conserved
    call assert(this%nstate.EQ.10,&
         msg='GQFPE static parameter nstate is not stored properly')
    call kill(this)    !destroy GQFPE to clean up memory
    call system('rm -f GQFPE.tmpfile*') !remove backup file(s)

    write(*,*)'test make sets correct default value for static parameter nstate'
    call make(this) !make GQFPE
    call assert(this%nstate.EQ.2,&
         msg='GQFPE default static parameter nstate is not 2')
    call kill(this)

    write(*,*)'test edge case for static parameter nstate breaks GQFPE'
    call make(this) !make GQFPE
    this%nstate=huge(this%nstate)
    call assert(check(this).NE.0,&
         msg='edge case huge value for static parameter nstate &
         &does not break GQFPE')
    call kill(this)

    write(*,*)'test edge case for static parameter nstate breaks GQFPE'
    call make(this) !make GQFPE
    this%nstate=0
    call assert(check(this).NE.0,&
         msg='edge case value 0 for static parameter nstate &
         &does not break GQFPE')
    call kill(this)

    write(*,*)'test basis can reproduce arbitrary 3 gaussian mixed state&
         & density with overlap>95%.'
    !define xgrid
    dq=(qmax-qmin)/real(npt,double)
    do ipt=0,npt-1
       q=qmin+(ipt*dq)
       grid(ipt)=q
    end do
    !first gaussian
    gaussweight(1)=1.0_double
    gaussparam(1,1)=1.1_double !mu
    gaussparam(1,2)=0.5_double !sigma
    !second gaussian
    gaussweight(2)=0.7_double
    gaussparam(2,1)=-0.7_double !mu
    gaussparam(2,2)=1.0_double !sigma
    !third gaussian
    gaussweight(3)=0.5_double
    gaussparam(3,1)=0.3_double !mu
    gaussparam(3,2)=0.2_double !sigma
    !accumulate and normalize density
    do ipt=0,npt-1
       psi(ipt)=sum(gaussweight(:)&
            *exp(-.5_double*((grid(ipt)-gaussparam(:,1))/gaussparam(:,2))**2))
    end do
    psi=psi/sqrt(sum(psi*psi))
    !define basis map (use consecutive HO states i.e. 1:1 mapping)
    do istate=0,nstate-1
       basismap(istate)=istate
    end do
    !calcualte basis functions
    do istate=0,nstate-1
       call HOwf(basismap(istate),grid,1.0_double,1.0_double,chi(istate,:))
       !check basis function
       do ipt=0,npt-1
          if(check(chi(istate,ipt)).NE.0)then
             write(*,*)'bad basis function at state',istate,ipt
             write(*,*)'HOstate',basismap(istate),chi(istate,ipt)
             stop
          end if
       end do
    end do
    !calculate expansion coefs
    do istate=0,nstate-1
       overlap(istate)=sum(psi(:)*chi(istate,:))
    end do
    !calculate mixed state
    do ipt=0,npt-1
       phi(ipt)=sum(overlap(:)*chi(:,ipt))
    end do
    call assert(sqrt(sum(psi*phi)).GT.0.98&
         ,msg='overlap of arb density and calculated density is less than 98%'&
         ,iostat=ierr)
    if(ierr.NE.0)then
       write(*,*)'Diagnostic: dumping arbitrary density'
       open(123,file='arbitrarydensity.dat')
       write(123,*)'# Compare density reproducibility'
       write(123,*)'# col1: xgrid'
       write(123,*)'# col2: target density (sum of 3 arb gaussians)'
       write(123,*)'# col3: linear combo N of Harmonic osc. states'
       write(123,*)'# col4: Nth Harmonic osc. wavefunction'
       do ipt=0,npt-1
          write(123,*)grid(ipt),psi(ipt),phi(ipt),chi(nstate-1,ipt)
       end do
       close(123)
    end if

    write(*,*)'test make allocates memory for dyanamic pointer array supermatrixL'
    call make(this)
    call assert(associated(this%supermatrixL),msg='GQFPE dynamic pointer array &
         &supermatrixL is not associated after make.')
    call kill(this)

    write(*,*)'test kill deallocates memory for dynamic pointer array supermatrixL!'
    call make(this)
    call kill(this)
    call assert(.not.associated(this%supermatrixL),msg='GQFPE dynamic pointer array &
         &supermatrixL remains associated after killed.')

    write(*,*)'test supermatrixL with dimensions not (nstate**2)x(nstate**2) breaks GQFPE.'
    call make(this)
    this%nstate=3 !change nstate without reseting GQFPE so supermatrixL is wrong size
    call assert(check(this).NE.0,msg='GQFPE does not break when supermatrixL dimensions&
         & are not (nstate**2)x(nstate**2)')
         !& are not (nstate**2)x3x3')
    call kill(this)
    
    write(*,*)'test dynamic pointer array supermatrixL can be resized by adjusting &
         & static parameter nstate then reseting with state=1'
    call make(this) !make GQFPE
    this%nstate=5 !adjust static parameter nstate to non-default value
    call reset(this,state=1) !reset GQFPE to reallocate dynamic memory
    !assert dynamic pointer array size has changed properly
    call assert(size(this%supermatrixL,1).EQ.5*5,&
         msg='GQFPE dynamic pointer array supermatrixL did not change size upon &
         &reset with state=1')
    call kill(this)    !destroy GQFPE to clean up memory

    write(*,*)'test edge case for dynamic pointer array supermatrixL breaks GQFPE'
    call make(this) !make GQFPE
    this%supermatrixL=Huge(1.0_double) !set edge case value
    !assert edge case value breaks GQFPE
    call assert(check(this).NE.0,&
         msg='edge case Huge values for dynamic pointer array supermatrixL &
         &does not break GQFPE')
    call kill(this)    !destroy GQFPE to clean up memory

    write(*,*)'Diagnostic: dumping supermatrixL at ts=5?'
    call make(this)
    this%nstate=8
    call reset(this,state=1)
    nstep=0
    do while(this%time*this%omega_c<5.0_double.and.nstep.lt.10000)
       call step(this) !perform one evolution step
       nstep=nstep+1
    end do
    open(123,file='supermatrixL.dat')
    !write(*,'(4X,25(I2,1X))')(jstate,jstate=0,size(this%supermatrixL,2)-1)
    do istate=0,size(this%supermatrixL,1)-1
       do jstate=0,size(this%supermatrixL,2)-1
          write(123,*)istate,jstate&
               ,real(this%supermatrixL(istate,jstate))!&
               !,aimag(this%supermatrixL(istate,jstate))
       end do
       write(123,*)
       !write(*,'(I2,3X,25(I1,2X))')istate,(int(this%supermatrixL(istate,jstate))&
       !     ,jstate=0,size(this%supermatrixL,2)-1)
    end do
    close(123)
    call kill(this)
    call system('gnuplot supermatrixL.plt')

    write(*,*)'test make allocates memory for dyanamic pointer array expancoef'
    call make(this)
    call assert(associated(this%expancoef),msg='GQFPE dynamic pointer array &
         &expancoef is not associated after make.')
    call kill(this)

    write(*,*)'test kill deallocates memory for dynamic pointer array expancoef!'
    call make(this)
    call kill(this)
    call assert(.not.associated(this%expancoef),msg='GQFPE dynamic pointer array &
         &expancoef remains associated after killed.')

    write(*,*)'test expancoef with dimensions not nstatexnstate breaks GQFPE.'
    call make(this)
    this%nstate=3 !change nstate without reseting GQFPE so expancoef is wrong size
    call assert(check(this).NE.0,msg='GQFPE does not break when expancoef dimensions&
         & are not nstatexnstate')
    call kill(this)
    
    write(*,*)'test dynamic pointer array expancoef can be resized by adjusting &
         & static parameter nstate then reseting with state=1'
    call make(this) !make GQFPE
    this%nstate=5 !adjust static parameter nstate to non-default value
    call reset(this,state=1) !reset GQFPE to reallocate dynamic memory
    !assert dynamic pointer array size has changed properly
    call assert(size(this%expancoef,1).EQ.5*5,&
         msg='GQFPE dynamic pointer array expancoef did not change size upon &
         &reset with state=1')
    call kill(this)    !destroy GQFPE to clean up memory

    write(*,*)'test edge case for dynamic pointer array expancoef breaks GQFPE'
    call make(this) !make GQFPE
    this%expancoef=Huge(1.0_double) !set edge case value
    !assert edge case value breaks GQFPE
    call assert(check(this).NE.0,&
         msg='edge case Huge values for dynamic pointer array expancoef &
         &does not break GQFPE')
    call kill(this)    !destroy GQFPE to clean up memory

    write(*,*)'test make allocates memory for dyanamic pointer array expancoefdot'
    call make(this)
    call assert(associated(this%expancoefdot),msg='GQFPE dynamic pointer array &
         &expancoefdot is not associated after make.')
    call kill(this)

    write(*,*)'test kill deallocates memory for dynamic pointer array expancoefdot!'
    call make(this)
    call kill(this)
    call assert(.not.associated(this%expancoefdot),msg='GQFPE dynamic pointer array &
         &expancoefdot remains associated after killed.')

    write(*,*)'test expancoefdot with dimensions not nstatexnstate breaks GQFPE.'
    call make(this)
    this%nstate=3 !change nstate without reseting GQFPE so expancoefdot is wrong size
    call assert(check(this).NE.0,msg='GQFPE does not break when expancoefdot dimensions&
         & are not nstatexnstate')
    call kill(this)
    
    write(*,*)'test dynamic pointer array expancoefdot can be resized by adjusting &
         & static parameter nstate then reseting with state=1'
    call make(this) !make GQFPE
    this%nstate=5 !adjust static parameter nstate to non-default value
    call reset(this,state=1) !reset GQFPE to reallocate dynamic memory
    !assert dynamic pointer array size has changed properly
    call assert(size(this%expancoefdot,1).EQ.5*5,&
         msg='GQFPE dynamic pointer array expancoefdot did not change size upon &
         &reset with state=1')
    call kill(this)    !destroy GQFPE to clean up memory

    write(*,*)'test edge case for dynamic pointer array expancoefdot breaks GQFPE'
    call make(this) !make GQFPE
    this%expancoefdot=Huge(1.0_double) !set edge case value
    !assert edge case value breaks GQFPE
    call assert(check(this).NE.0,&
         msg='edge case Huge values for dynamic pointer array expancoefdot &
         &does not break GQFPE')
    call kill(this)    !destroy GQFPE to clean up memory


    !density tests
    write(*,*)'test make allocates memory for dyanamic pointer array density'
    call make(this)
    call assert(associated(this%density),msg='GQFPE dynamic pointer array &
         &density is not associated after make.')
    call kill(this)

    write(*,*)'test kill deallocates memory for dynamic pointer array density'
    call make(this)
    call kill(this)
    call assert(.not.associated(this%density),msg='GQFPE dynamic pointer array &
         &density remains associated after killed.')

    write(*,*)'test dynamic pointer array density is saved properly in backup file'
    call make(this) !make GQFPE
    this%density=eye    !manually set dynamic pointer array to non-default value
    call system('rm -f GQFPE.tmpfile*') !remove any previous backup file(s)
    call backup(this,file='GQFPE.tmpfile') !create backup file
    call kill(this) !destroy GQFPE
    call make(this,file='GQFPE.tmpfile') !make new GQFPE from backup file
    !assert non default parameter is conserved
    call assert(all(this%density.EQ.eye),&
         msg='GQFPE dynamic pointer array density is not stored properly')
    call kill(this)    !destroy GQFPE to clean up memory
    call system('rm -f GQFPE.tmpfile*') !remove backup file(s)

    write(*,*)'test make dynamic pointer array density with non default nstate&
         & will break GQFPE when not reset with state=1'
    call make(this) !make GQFPE
    this%nstate=5
    call assert(check(this).NE.0,msg='changing nstate without call to reset did not break&
         & GQFPE')
    call kill(this)    !destroy GQFPE to clean up memory
    
    write(*,*)'test dynamic pointer array density can be resized by adjusting &
         & static parameter nstate then reseting with state=1'
    call make(this) !make GQFPE
    this%nstate=5 !adjust static parameter nstate to non-default value
    call reset(this,state=1) !reset GQFPE to reallocate dynamic memory
    !assert dynamic pointer array size has changed properly
    call assert(size(this%density,1).EQ.5,&
         msg='GQFPE dynamic pointer array density did not change size upon &
         &reset with state=1')
    call kill(this)    !destroy GQFPE to clean up memory

    write(*,*)'test default value for dynamic pointer array density'
    call make(this) !make GQFPE
    !assert defualt value is superposition of first two states
    call assert(all(this%density.EQ.(0.5_double,0.0_double)),&
         msg='default initial state is not state1+state2 superposition')
    call kill(this)    !destroy GQFPE to clean up memory

    write(*,*)'test edge case for dynamic pointer array density breaks GQFPE'
    call make(this) !make GQFPE
    this%density=huge(1.0_double) !set edge case value
    !assert edge case value breaks GQFPE
    call assert(check(this).NE.0,&
         msg='edge case huge value for dynamic pointer array density &
         &does not break GQFPE')
    call kill(this)    !destroy GQFPE to clean up memory

    write(*,*)'Diagnostic: dumping density simlulation for runtime ts=[0:50]?'
    call make(this)
    N=5
    this%omega_c=200.0_double*invcm
    !this%omega_c=1.0_double/Eh 
    this%beta=5.0_double/this%omega_c
    this%nstate=50!2**N
    !nyquist time
    this%dt=1.0_double/(2.0_double*hbar*this%omega_c*(this%nstate+0.5_double))
    !this%dt=0.005_double/this%omega_c
    call reset(this,state=1)
    write(*,*)'cutoff frequency in hartrees =',this%omega_c
    write(*,*)'cutoff frequency in wavenumbers =',this%omega_c/invcm
    write(*,*)'temperature in Kelvin=',1.0_double/(this%beta*kB)
    write(*,*)'kBT in hartrees=',1.0_double/this%beta
    write(*,*)'kBT in wavenumbers=',1.0_double/this%beta/invcm
    write(*,*)'bath friction in hartrees =',0.1*this%omega_c
    write(*,*)'bath friction in wavenumbers =',0.1*this%omega_c/invcm
    write(*,*)'simulation time step in au =',this%dt
    write(*,*)'simulation time step in fs =',this%dt/fs
    write(*,*)'simulation step frequency in wavenumbers =',1.0_double/this%dt/invcm
    write(*,*)'simulation step frequency in hartrees/hbar =',1.0_double/this%dt
    write(*,*)'Nyquist frequency in wavenumbers =',0.5_double/this%dt/invcm
    write(*,*)'Nyquist frequency in hartrees/hbar =',0.5_double/this%dt
    write(*,*)'highest HO state energy in wavenumbers ='&
         ,hbar*this%omega_c*(this%nstate+.5)/invcm
    write(*,*)'highest HO state energy in hartrees ='&
         ,hbar*this%omega_c*(this%nstate+.5)
    write(*,*)'-------scaled units-------'
    write(*,*)'scaled inverse temperature=',this%beta*hbar*this%omega_c
    write(*,*)'step of scaled time=',this%dt*this%omega_c
    write(*,*)'size of HO oscillator basis=',this%nstate
 !stop
   open(123,file='density.dat')
    !retrieve density from expansion coefs
!    do lstate=0,this%nstate**2-1
!       statepair=spacefilling_index2coord(this%nstate,lstate,type='rowmajor')
!       this%density(statepair(1),statepair(2))=this%expancoef(lstate)
!    end do
!    write(123,*)this%time*this%omega_c,(&
!         real(this%density(istate,istate)),istate=0,this%nstate-1)
    do istate=0,ndump-1
       statepair=istate
       lstate=spacefilling_coord2index(this%nstate,statepair,type='rowmajor')
       this%density(istate,istate)=this%expancoef(lstate)
    end do
    write(123,*)this%time*this%omega_c,(&
         real(this%density(istate,istate)),istate=0,ndump-1)
    nstep=0
    do while(this%time*this%omega_c<50.0_double.and.nstep.lt.10000)
       call step(this) !perform one evolution step
       nstep=nstep+1

       write(*,*)nstep,'runtime=',this%time*this%omega_c
       
       if(mod(nstep,10).EQ.0)then !dump density
          !retrieve density from expansion coefs
          do lstate=0,this%nstate**2-1
             statepair=spacefilling_index2coord(this%nstate,lstate,type='rowmajor')
             this%density(statepair(1),statepair(2))=this%expancoef(lstate)
          end do
          write(123,*)this%time*this%omega_c,(&
               real(this%density(istate,istate)),istate=0,7)!this%nstate-1)
          flush(123)
       end if
    end do
    !write(*,'(4X,25(I2,1X))')(jstate,jstate=0,size(this%density,2)-1)
    !do istate=0,size(this%density,1)-1
    !   do jstate=0,size(this%density,2)-1
    !      write(123,*)istate,jstate&
    !           ,real(this%density(istate,jstate))!&
    !           !,aimag(this%density(istate,jstate))
    !   end do
    !   write(123,*)
    !   !write(*,'(I2,3X,25(I1,2X))')istate,(int(this%density(istate,jstate))&
    !   !     ,jstate=0,size(this%density,2)-1)
    !end do
    close(123)
    call kill(this)
    call system('gnuplot density.plt')




    !================== consider the following tests for =====================
    !=========================================================================

    !==================       static parameters             ==================

    !write(*,*)'test static parameter NNN is stored properly in backup file'
    !call make(this) !make GQFPE
    !this%NNN=MMM    !manually set static parameter to non-default value
    !call reset(this,state=1)
    !call system('rm -f GQFPE.tmpfile*') !remove any previous backup file(s)
    !call backup(this,file='GQFPE.tmpfile') !create backup file
    !call kill(this) !destroy GQFPE
    !call make(this,file='GQFPE.tmpfile') !make new GQFPE from backup file
    !!assert non default parameter is conserved
    !call assert(this%NNN.EQ.MMM,&
    !     msg='GQFPE static parameter NNN is not stored properly')
    !call kill(this)    !destroy GQFPE to clean up memory
    !call system('rm -f GQFPE.tmpfile*') !remove backup file(s)

    !write(*,*)'test make sets correct default value for static parameter NNN'
    !call make(this) !make GQFPE
    !call assert(this%NNN.EQ.MMM,&
    !     msg='GQFPE default static parameter NNN is not MMM')
    !call kill(this)

    !write(*,*)'test edge case for static parameter NNN breaks GQFPE'
    !call make(this) !make GQFPE
    !this%NNN=MMMM
    !call assert(check(this).NE.0,&
    !     msg='edge case value MMM for static parameter NNN &
    !     &does not break GQFPE')
    !call kill(this)

    !==================    dynamic arrays and pointers      ==================

    !write(*,*)'test make allocates memory for dyanamic pointer array PPP'
    !call make(this)
    !call assert(associated(this%PPP),msg='GQFPE dynamic pointer array &
    !     &PPP is not associated after make.')
    !call kill(this)

    !write(*,*)'test kill deallocates memory for dynamic pointer array PPP'
    !call make(this)
    !call kill(this)
    !call assert(.not.associated(this%PPP),msg='GQFPE dynamic pointer array &
    !     &PPP remains associated after killed.')

    !write(*,*)'test dynamic pointer array PPP is saved properly in backup file'
    !call make(this) !make GQFPE
    !this%PPP=YYY    !manually set dynamic pointer array to non-default value
    !call system('rm -f GQFPE.tmpfile*') !remove any previous backup file(s)
    !call backup(this,file='GQFPE.tmpfile') !create backup file
    !call kill(this) !destroy GQFPE
    !call make(this,file='GQFPE.tmpfile') !make new GQFPE from backup file
    !!assert non default parameter is conserved
    !call assert(all(this%PPP.EQ.YYY),&
    !     msg='GQFPE dynamic pointer array PPP is not stored properly')
    !call kill(this)    !destroy GQFPE to clean up memory
    !call system('rm -f GQFPE.tmpfile*') !remove backup file(s)

    !write(*,*)'test make allocates dynamic pointer array PPP of default size'
    !call make(this) !make GQFPE
    !!assert dynamic pointer array has proper size for all dimensions
    !call assert(size(this%PPP,J).EQ.N,&
    !     msg='GQFPE dynamic pointer array PPP is not of size N for &
    !     &dimension J')
    !call assert(size(this%PPP,I).EQ.N,&
    !     msg='GQFPE dynamic pointer array PPP is not of size N for &
    !     &dimension I')
    !call kill(this)    !destroy GQFPE to clean up memory
    
    !write(*,*)'test dynamic pointer array PPP can be resized by adjusting &
    !     & static parameter NNN then reseting with state=1'
    !call make(this) !make GQFPE
    !this%NNN=MMM !adjust static parameter NNN to non-default value
    !call reset(this,state=1) !reset GQFPE to reallocate dynamic memory
    !!assert dynamic pointer array size has changed properly
    !call assert(size(this%PPP,J).EQ.MMM,&
    !     msg='GQFPE dynamic pointer array PPP did not change size upon &
    !     &reset with state=1')
    !call kill(this)    !destroy GQFPE to clean up memory

    !write(*,*)'test edge case for dynamic pointer array PPP breaks GQFPE'
    !call make(this) !make GQFPE
    !this%PPP=YYY !set edge case value
    !assert edge case value breaks GQFPE
    !call assert(check(this).NE.0,&
    !     msg='edge case value YYY for dynamic pointer array PPP &
    !     &does not break GQFPE')
    !call kill(this)    !destroy GQFPE to clean up memory

    !==================        calculated variables         ==================

    !write(*,*)'test make sets default value for calculated variable XXX'
    !call make(this) !make GQFPE
    !call assert(this%XXX.EQ.YYY,&
    !     msg='GQFPE default calculated variable XXX is not YYY')
    !call kill(this)

    !write(*,*)'test edge case for calculated variable XXX breaks GQFPE'
    !call make(this) !make GQFPE
    !this%XXX=YYY
    !call assert(check(this).NE.0,&
    !     msg='edge case value YYY for calculated variable XXX &
    !     &does not break GQFPE')
    !call kill(this)
    !========================================================================


    write(*,*)'ALL GQFPE TESTS PASSED!'
  end subroutine GQFPE_test
end module GQFPE_class
