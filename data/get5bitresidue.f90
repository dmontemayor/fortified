! program will read residuemap.txt and use it to convert protstruct file
! into 5bit by 2*nwindow+1 residue window code
program main
  implicit none

  integer,parameter::nres=20,nstruct=3,nwindow=4
  character(len=50)::filename
  integer::len,ierr

  !residue character symbols
  character(len=1)::reschar(nres)

  !residue 5bit code and structure 2bit code 
  integer::resbit(5,0:nres),structbit(2,0:nstruct)

  !residue and structure window
  integer::reswindow(-nwindow:nwindow)
  integer::structwindow(-nwindow:nwindow)

  !book keeping
  integer::input5bit(5), nrecord, ires, istruct, jres, nseq, seqlen, maxseqlen
  character(len=1)::char1
  character(len=2)::char2
  logical::resok

  !load residue map
  open(111,file='residuemap.txt', action='read', status='old')
  read(111,*,iostat=ierr)!skip header line
  !read null resiude
  if(ierr.EQ.0)then
     read(111,FMT='(5(I1),1X,1A)',iostat=ierr) input5bit, char1
     nrecord=0
     resbit(:,nrecord)=input5bit
     reschar(nrecord)=char1
  end if
  do while(ierr.EQ.0)
     read(111,FMT='(5(I1),1X,1A)',iostat=ierr) input5bit, char1
     if(ierr.EQ.0.and.char1.NE.'X')then
        nrecord=nrecord+1
        resbit(:,nrecord)=input5bit
        reschar(nrecord)=char1
     end if
  end do
  close(111)

  !stop if all N residues were not assigned
  if(nrecord.NE.nres)then
     write(*,*)'ERROR: did not assign all residues.'
     write(*,*)'nrecord=',nrecord
     write(*,*)'program will stop.'
     stop
  end if

  !laod 2bit struct code
  !struct 0: null struct
  structbit(1,0)=0
  structbit(2,0)=0
  !struct 1: random coil '_'
  structbit(1,1)=1
  structbit(2,1)=1
  !struct 2: beta sheet 'e'
  structbit(1,2)=0
  structbit(2,2)=1
  !struct 3: alpha-helix 'h'
  structbit(1,3)=1
  structbit(2,3)=0

  !get input file
  call get_command_argument(1,filename,len,ierr)
  if(ierr.GT.0)write(*,*)'file name cannot be read properly! Program will stop.'
  if(ierr.LT.0)write(*,*)trim(filename)//&
       'file name was truncated!Program will stop.'
  if(ierr.NE.0)stop
  
  !open output file
  open(222,file=trim(filename)//'.dat')
  write(*,*)'program will output to file: '//trim(filename)//'.dat'

  !read input file and translate to output
  open(111,file=trim(filename),iostat=ierr)
  if(ierr.EQ.0) write(*,*)'reading input file ',trim(filename)
  nrecord=0
  nseq=0
  seqlen=0
  maxseqlen=0
  !initialize window with null residues and null structures
  reswindow=0
  structwindow=0
  do while(ierr.EQ.0)
     read(111,'(1A,2A)',iostat=ierr) char1,char2
     if(ierr.EQ.0)then
        nrecord=nrecord+1
        
        if(char1.EQ.'<')then
           !sequence breaks
           nseq=nseq+1
           if(seqlen.GT.maxseqlen)maxseqlen=seqlen
           seqlen=0

           !dump end of previous sequence if present
           do while(reswindow(0).NE.0)
              
              !push reisude and structure to window
              call push(0,reswindow)
              call push(0,structwindow)
              
              if(reswindow(0).NE.0)write(222,FMT='(47(I1,1X))')&
                   (resbit(:,reswindow(jres)),jres=-nwindow,nwindow)&
                   ,structbit(:,structwindow(0))
           end do

           !output sequence break
           write(222,*)'<>'

           !reinitialize window with null residues and null structures
           reswindow=0
           structwindow=0

        else

           !translate residue
           resok=.false.
           ires=0
           do jres=1,nres
              if(char1.EQ.reschar(jres))ires=jres
           end do

           if(ires.NE.0)resok=.true.           
           !check residue was assigned correctly
           !if(ires.EQ.0)then
           !   write(*,*)'unknown residue symbol encountered in record',nrecord
           !   write(*,*)'valid residue symbols are:'
           !   write(*,*)reschar(1:nres)
           !   write(*,*)'program will stop'
           !   stop
           !end if

           if(resok)then

              !translate structure
              select case (adjustl(char2))
              case ('_') !random coil
                 istruct=1
              case ('e') !beta sheet
                 istruct=2
              case ('h') !alpha-helix
                 istruct=3
              case default !null struct
                 istruct=0
              end select
                         
              !check structure was assigned correctly
              if(istruct.EQ.0)then
                 write(*,*)'unknown structure label encountered in record'&
                      ,nrecord
                 write(*,*)'unkown label is "'//adjustl(char2)//'"'
                 write(*,*)'valid labels are: "_", "e", and "h"'
                 write(*,*)'program will stop'
                 stop
              end if
           
              seqlen=seqlen+1

              !push reisude and structure to window
              call push(ires,reswindow)
              call push(istruct,structwindow)
              
              !output after residue window has been primed
              !(element 0 is not null)
              if(reswindow(0).NE.0)write(222,FMT='(47(I1,1X))')&
                   (resbit(:,reswindow(jres)),jres=-nwindow,nwindow)&
                   ,structbit(:,structwindow(0))
           end if
     
        end if
     end if
  end do
  !close output file
  close(222)
  write(*,*)'program found',nseq-1,'sequences'
  write(*,*)'longest sequence residue length',maxseqlen
  !close input file
  close(111)

contains
  subroutine push(I,A)
    !discards first element of A then
    !shifts all elements of A down by 1 then
    !assigns I to last element of A
    implicit none
    
    integer,intent(in)::I
    integer,intent(inout)::A(:)
    integer::N,j
    
    N=size(A)
    
    !shifts all elements of A down by 1 then
    do j=1,N-1
       A(j)=A(j+1)
    end do
    
    !assigns I to last element of A
    A(N)=I
    
  end subroutine push
  
end program


