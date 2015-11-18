module MDutils
  implicit none
  private

  public::secondarystruct,residuetable

  type residue
     logical::bitcode(5)
     character(len=1)::symbol
     character(len=4)::resname
  end type residue

  type(residue)::residuetable(0:2**5-1)

  data residuetable(00)%bitcode/.false.,.false.,.false.,.false.,.false./
  data residuetable(00)%symbol/ '_'/
  data residuetable(00)%resname/'NUL '/!null

  data residuetable(01)%bitcode/.false.,.false.,.false.,.false.,.true./
  data residuetable(01)%symbol/ 'X'/
  data residuetable(01)%resname/'XXX '/!unassigned

  data residuetable(02)%bitcode/.false.,.false.,.false.,.true.,.false./
  data residuetable(02)%symbol/ 'X'/
  data residuetable(02)%resname/'XXX '/!unassigned

  data residuetable(03)%bitcode/.false.,.false.,.false.,.true.,.true./
  data residuetable(03)%symbol/ 'A'/
  data residuetable(03)%resname/'ALA '/!alanine

  data residuetable(04)%bitcode/.false.,.false.,.true.,.false.,.false./
  data residuetable(04)%symbol/ 'X'/
  data residuetable(04)%resname/'XXX '/!unassigned

  data residuetable(05)%bitcode/.false.,.false.,.true.,.false.,.true./
  data residuetable(05)%symbol/ 'R'/
  data residuetable(05)%resname/'ARG '/!arginine

  data residuetable(06)%bitcode/.false.,.false.,.true.,.true.,.false./
  data residuetable(06)%symbol/ 'N'/
  data residuetable(06)%resname/'ASN '/!asparagine

  data residuetable(07)%bitcode/.false.,.false.,.true.,.true.,.true./
  data residuetable(07)%symbol/ 'D'/
  data residuetable(07)%resname/'ASP '/!apartic_acid

  data residuetable(08)%bitcode/.false.,.true.,.false.,.false.,.false./
  data residuetable(08)%symbol/ 'X'/
  data residuetable(08)%resname/'XXX '/!unassigned

  data residuetable(09)%bitcode/.false.,.true.,.false.,.false.,.true./
  data residuetable(09)%symbol/ 'C'/
  data residuetable(09)%resname/'CYS '/!cysteine

  data residuetable(10)%bitcode/.false.,.true.,.false.,.true.,.false./
  data residuetable(10)%symbol/ 'Q'/
  data residuetable(10)%resname/'GLN '/!glutamine

  data residuetable(11)%bitcode/.false.,.true.,.false.,.true.,.true./
  data residuetable(11)%symbol/ 'E'/
  data residuetable(11)%resname/'GLU '/!glutamic_acid

  data residuetable(12)%bitcode/.false.,.true.,.true.,.false.,.false./
  data residuetable(12)%symbol/ 'G'/
  data residuetable(12)%resname/'GLY '/!glycine

  data residuetable(13)%bitcode/.false.,.true.,.true.,.false.,.true./
  data residuetable(13)%symbol/ 'H'/
  data residuetable(13)%resname/'HIS '/!histidine

  data residuetable(14)%bitcode/.false.,.true.,.true.,.true.,.false./
  data residuetable(14)%symbol/ 'I'/
  data residuetable(14)%resname/'ILE '/!isoleucine

  data residuetable(15)%bitcode/.false.,.true.,.true.,.true.,.true./
  data residuetable(15)%symbol/ 'X'/
  data residuetable(15)%resname/'XXX '/!unassigned

  data residuetable(16)%bitcode/.true.,.false.,.false.,.false.,.false./
  data residuetable(16)%symbol/ 'X'/
  data residuetable(16)%resname/'XXX '/!unassigned

  data residuetable(17)%bitcode/.true.,.false.,.false.,.false.,.true./
  data residuetable(17)%symbol/ 'L'/
  data residuetable(17)%resname/'LEU '/!leucine

  data residuetable(18)%bitcode/.true.,.false.,.false.,.true.,.false./
  data residuetable(18)%symbol/ 'K'/
  data residuetable(18)%resname/'LYS '/!lysine

  data residuetable(19)%bitcode/.true.,.false.,.false.,.true.,.true./
  data residuetable(19)%symbol/ 'M'/
  data residuetable(19)%resname/'MET '/!methionine

  data residuetable(20)%bitcode/.true.,.false.,.true.,.false.,.false./
  data residuetable(20)%symbol/ 'F'/
  data residuetable(20)%resname/'PHE '/!phenylalanine

  data residuetable(21)%bitcode/.true.,.false.,.true.,.false.,.true./
  data residuetable(21)%symbol/ 'P'/
  data residuetable(21)%resname/'PRO '/!proline

  data residuetable(22)%bitcode/.true.,.false.,.true.,.true.,.false./
  data residuetable(22)%symbol/ 'S'/
  data residuetable(22)%resname/'SER '/!serine

  data residuetable(23)%bitcode/.true.,.false.,.true.,.true.,.true./
  data residuetable(23)%symbol/ 'X'/
  data residuetable(23)%resname/'XXX '/!unassigned

  data residuetable(24)%bitcode/.true.,.true.,.false.,.false.,.false./
  data residuetable(24)%symbol/ 'T'/
  data residuetable(24)%resname/'THR '/!threonine

  data residuetable(25)%bitcode/.true.,.true.,.false.,.false.,.true./
  data residuetable(25)%symbol/ 'W'/
  data residuetable(25)%resname/'TRP '/!tryptophan

  data residuetable(26)%bitcode/.true.,.true.,.false.,.true.,.false./
  data residuetable(26)%symbol/ 'Y'/
  data residuetable(26)%resname/'TYR '/!tyrosine

  data residuetable(27)%bitcode/.true.,.true.,.false.,.true.,.true./
  data residuetable(27)%symbol/ 'X'/
  data residuetable(27)%resname/'XXX '/!unassigned

  data residuetable(28)%bitcode/.true.,.true.,.true.,.false.,.false./
  data residuetable(28)%symbol/ 'V'/
  data residuetable(28)%resname/'VAL '/!valine

  data residuetable(29)%bitcode/.true.,.true.,.true.,.false.,.true./
  data residuetable(29)%symbol/ 'X'/
  data residuetable(29)%resname/'XXX '/!unassigned

  data residuetable(30)%bitcode/.true.,.true.,.true.,.true.,.false./
  data residuetable(30)%symbol/ 'X'/
  data residuetable(30)%resname/'XXX '/!unassigned

  data residuetable(31)%bitcode/.true.,.true.,.true.,.true.,.true./
  data residuetable(31)%symbol/ 'X'/
  data residuetable(31)%resname/'XXX '/!unassigned


  interface secondarystruct
     module procedure bool2struct
     module procedure int2struct
  end interface secondarystruct


contains
  !--------------------------------------
  function int2struct(bit1,bit2)
    ! wrapper for bool2struct to accept integer input
    !integer values: 0 for false, true otherwise 
    integer,intent(in)::bit1,bit2
    character(len=1)::int2struct
    logical::lbit1,lbit2
    lbit1=.false.
    lbit2=.false.
    if(bit1.NE.0)lbit1=.true.
    if(bit2.NE.0)lbit2=.true.
    int2struct=bool2struct(lbit1,lbit2)
  end function int2struct
  !--------------------------------------
  function bool2struct(bit1,bit2)
    ! returns character symbol for protein secondary structure
    ! contained in two input logical values
    ! return values are:
    ! 'a' = alpha-helix (1,0)
    ! 'b' = beta sheet  (0,1)
    ! 'c' = random coil (1,1)
    ! 'x' = null        (0,0)
    logical,intent(in)::bit1,bit2
    character(len=1)::bool2struct
    if(bit1)then!a or c
       if(bit2)then!c
          bool2struct='c'
       else!a
          bool2struct='a'
       end if
    else!b or x
       if(bit2)then!b
          bool2struct='b'
       else!x
          bool2struct='x'
       end if
    end if
  end function bool2struct
end module MDutils
