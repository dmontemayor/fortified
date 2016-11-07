!===========================================================
!> \brief
!! Definition of atomic units.
!! \details
!! Atomic units are the default used throughout nonadmd.
!! \author
!! Daniel Montemayor
!<===========================================================
module atomicunits
  use type_kinds
  use math
  implicit none

  !Fundamental Units
  real(double),parameter::a0=1.0_double   !<(Fundamental length) Bohr radius
  real(double),parameter::me=1.0_double   !<(Fundamental mass) electron mass
  real(double),parameter::hbar=1.0_double !<(Fundamental action) reduced Planks const
  real(double),parameter::e=1.0_double    !<(Fundamental charge) electron charge
  real(double),parameter::Kel=1.0_double  !<(Fundamental temperature) Kelvin
  real(double),parameter::mol=1.0_double  !<(Fundamental amount) Avogadro's number
  real(double),parameter::rad=1.0_double  !<(Fundamental angle) radian

  !Derived Units
  real(double),parameter::Eh=(hbar/a0)**2/me                      !<Hartree (unit value) 
  real(double),parameter::kC=Eh*a0/(e**2)                         !<Coulomb force (unit value)
  real(double),parameter::mp=1836.152663_double*me                !<proton rest mass
  real(double),parameter::mn=1838.685239_double*me                !<neutron rest mass
  real(double),parameter::kb=3.166815203_double*1E-6*Eh/Kel       !<boltzman constant
  real(double),parameter::eV=3.674932587_double*1E-2*Eh           !<Electon volt
  real(double),parameter::c0=137.0359996287515_double*hbar/(a0*me)!<speed of light
  real(double),parameter::eps0=1.0_double/(4*pi*kC)             !<vacuum permittivity
  real(double),parameter::Debye=.3934302014076827_double*e*a0     !<Dipole moment
  real(double),parameter::Deg=twopi/360_double*rad                !<Degrees

  !Converstion Factors
  real(double),parameter::angstrom=1.889726133921252_double*a0    !<angstrom
  real(double),parameter::fs=41.34137337414122_double*hbar/Eh     !<femtosecond
  real(double),parameter::ps=fs*1_double*1E3                      !<picosecond
  real(double),parameter::invcm=4.5663352672_double*1E-6*Eh/hbar  !<wavenumber(1/cm)

  
end module atomicunits
