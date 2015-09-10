!============================================================
!> \brief
!! Definition of integer, floating point, and string kinds.
!! \details
!! Cross-Platform portability is the goal such that type
!! precsision is maintained from build to build of the software.
!! \authors
!! Author: Daniel Montemayor 2010
!<===========================================================
module type_kinds
  implicit none

  Private

  !Integer Kinds
  Public :: byte,short,long

  !Floating point Kinds
  Public :: single,double

  !String Lengths
  Public :: line,extline,label,path,comment 

  !> string length for a line in a file
  integer, parameter:: line=72
  !> string length for an extended line in a file
  integer, parameter:: extline=132
  !> string length for a name/title/label
  integer, parameter:: label=31 
  !> string length for a directory path
  integer, parameter:: path=255
  !> string length for a message/comment
  integer, parameter:: comment=255

  !> byte integer size 
  integer, parameter:: byte=selected_int_kind(1)
  !> short integer size 
  integer, parameter:: short=selected_int_kind(4)
  !> long integer size 
  integer, parameter:: long=selected_int_kind(8)

  !> single precision float size 32-bit
  integer, parameter:: single=selected_real_kind(6,37)
  !> double precision float size 64-bit
  integer, parameter:: double=selected_real_kind(15,307)
  !!> quad precision float size 128-bit
  !integer, parameter:: quad=selected_real_kind(33,4931)

end module type_kinds
