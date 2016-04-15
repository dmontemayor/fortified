!> \brief
!! Get strings from numbers and other string operations
module string
  use type_kinds
!  use MPIFRAMEWORK
!  use filemanager
  private
  
  Public::quote
  Public::int2str,float2str,int2hex
  Public::str2int,str2float,hex2int

  interface int2str
     module procedure long2str
     module procedure short2str
  end interface
  interface float2str
     module procedure double2str
     module procedure single2str
  end interface

  interface str2int
     module procedure str2long
  end interface
  interface str2float
     module procedure str2double
  end interface
  interface hex2int
     module procedure hex2long
  end interface
  interface int2hex
     module procedure long2hex
  end interface

contains
  !----------------------------
  function quote(str)
    character*(*),intent(in)::str
    character(len=len(str)+2)::quote

    quote="'"//str//"'"
  end function quote
  !----------------------------
  function long2str(I)

    integer(long),intent(in)::I
    character(len=label)::long2str
    integer(long)::ierr

    if((I.NE.I).or.(abs(I).GT.huge(I)))then
       long2str='int2str_ERROR'
    else
       write(long2str,'(I10)',iostat=ierr)I
       if(ierr.NE.0)long2str='int2str_ERROR'
    end if

    long2str=adjustl(long2str)

  end function long2str
  !----------------------------
  function short2str(I)

    integer(short),intent(in)::I
    integer(long)::longI
    character(len=label)::short2str

    longI=I
    short2str=trim(int2str(longI))
    
  end function short2str
  !----------------------------
  function str2long(string)

    character*(*),intent(in)::string
    integer(long)::str2long,order,i
    integer(long)::ierr
    logical::tryhex

    !make sure all characters are appropriate
    order=len(string)
    tryhex=.false.
    do i=1,order
       select case(string(i:i))
       case(' ')
       case('0')
       case('1')
       case('2')
       case('3')
       case('4')
       case('5')
       case('6')
       case('7')
       case('8')
       case('9')
       case default
          !write(*,*)'str2int error: '//string&
          !     //' contains non integer characters!'
          tryhex=.true.
       end select
    end do
    if(tryhex)then
       !attempt hexadecimal
       str2long=hex2long(string)
       !write(*,*)'in hex: '//string//' =',str2long
    else
       read(string,*,iostat=ierr)str2long
       if(ierr.NE.0)then
          write(*,*)'str2int error: '//string//' is not an integer! Program will stop!'
          stop
       end if
    end if
  end function str2long
  !----------------------------
  function double2str(X)

    real(double),intent(in)::X
    character(len=label)::double2str
    integer(long)::ierr

    if((X.NE.X).or.(abs(X).GT.huge(X)))then
       double2str='float2str_ERROR'
    else
       write(double2str,'(E14.6E3)',iostat=ierr)X
       if(ierr.NE.0)double2str='float2str_ERROR'
    end if

    double2str=adjustl(double2str)
  end function double2str
  !----------------------------
  function str2double(string)

    character*(*),intent(in)::string
    real(double)::str2double
    integer(long)::ierr

    read(string,*,iostat=ierr)str2double
    if(ierr.NE.0)then
       write(*,*)'str2float error: '//string//' is not a float! Program will stop!'
       stop
    end if

  end function str2double
  !----------------------------
  function single2str(X)
    
    real(single),intent(in)::X
    real(double)::doubleX
    character(len=label)::single2str
    
    doubleX=X
    single2str=trim(double2str(doubleX))
    
  end function single2str
  !----------------------------
  function hex2long(string)

    character*(*),intent(in)::string
    integer(long)::hex2long,order,i,j

    !get length of hex
    order=len(string)

    !make sure all characters are appropriate
    do i=1,order
       select case(string(i:i))
       case('0')
       case('1')
       case('2')
       case('3')
       case('4')
       case('5')
       case('6')
       case('7')
       case('8')
       case('9')
       case('a')
       case('b')
       case('c')
       case('d')
       case('e')
       case('f')
       case('A')
       case('B')
       case('C')
       case('D')
       case('E')
       case('F')
       case default
          write(*,*)'hex2int error: '//string//' is not a hexadecimal number! Program will stop.'
          stop
       end select
    end do

    hex2long=0
    do i=1, order
       select case(string(i:i))
       case('0')
          j=0
       case('1')
          j=1
       case('2')
          j=2
       case('3')
          j=3
       case('4')
          j=4
       case('5')
          j=5
       case('6')
          j=6
       case('7')
          j=7
       case('8')
          j=8
       case('9')
          j=9
       case('a')
          j=10
       case('b')
          j=11
       case('c')
          j=12
       case('d')
          j=13
       case('e')
          j=14
       case('f')
          j=15
       case('A')
          j=10
       case('B')
          j=11
       case('C')
          j=12
       case('D')
          j=13
       case('E')
          j=14
       case('F')
          j=15
       end select
       hex2long=hex2long+j*16**(order-i)
    end do
  end function hex2long
  !----------------------------
  function long2hex(I)

    integer(long),intent(in)::I
    character(len=label)::long2hex
    character::char
    integer(long)::order,j,remainder,digit

    !get length of hex
    order=0
    do while(I.GE.16**order)
       order=order+1
    end do

    remainder=I
    long2hex=''
    do j=1, order
       digit=0
       do while(remainder.GE.(digit+1)*16**(order-j))
          digit=digit+1
       end do
       remainder=remainder-digit*16**(order-j)
       if(digit.LT.10)then
          char=trim(int2str(digit))
       else
          select case(digit)
          case(10)
             char='a'
          case(11)
             char='b'
          case(12)
             char='c'
          case(13)
             char='d'
          case(14)
             char='e'
          case(15)
             char='f'
          case default
             write(*,*)'error hex2int: cannot convert integer to hexadecimal! Program will stop!'
             stop
          end select
       end if 
          long2hex=trim(long2hex)//char
    end do
  end function long2hex
end module string
