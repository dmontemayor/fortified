module MDutils
  implicit none
  private

  public::secondarystruct


contains
  function secondarystruct (bit1,bit2)
    ! returns character symbol for protein secondary structure
    ! contained in two input bit values
    ! return values are:
    ! 'a' = alpha-helix (1,0)
    ! 'b' = beta sheet  (0,1)
    ! 'c' = random coil (1,1)
    ! 'x' = null        (0,0)
    logical,intent(in)::bit1,bit2
    character(len=1)::secondarystruct
    if(bit1)then!a or c
       if(bit2)then!c
          secondarystruct='c'
       else!a
          secondarystruct='a'
       end if
    else!b or x
       if(bit2)then!b
          secondarystruct='b'
       else!x
          secondarystruct='x'
       end if
    end if
  end function secondarystruct
end module MDutils
