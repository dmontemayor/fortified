ACLOCAL_AMFLAGS = -I m4 --install
EXTRA_DIST = m4 data

corefiles = src/type_kinds.f90 src/math.f90 src/rand.f90 src/testing.f90 src/functions.f90 src/progressbar.f90 src/filemanager.f90 src/string.f90 src/template.f90

testfiles = src/coretests.f90

auxfiles = src/MDutils.f90 src/atomicunits.f90

# All fortified executable modules
bin_PROGRAMS = nonadmd cann test

###### CANN chemistry artificial neural networks ###### 	
cann_FCFLAGS =  -fopenmp -J src/
cannfiles = src/layer.f90 src/bpg.f90 src/ffn.f90 src/som.f90 src/device.f90
cann_SOURCES = $(corefiles) $(auxfiles) $(cannfiles) src/cann_driver.f90

###### NONADMD Non-Adiabatic Molecular Dynamics ###### 	
nonadmd_FCFLAGS =  -fopenmp -J src/
nonadmdfiles = src/quasiparticle.f90 src/ham.f90 src/wavelet.f90 src/metiu.f90
nonadmd_SOURCES = $(corefiles) $(auxfiles) $(nonadmdfiles) src/nonadmd_driver.f90

###### TEST module ###### 	
test_FCFLAGS = -g -fbounds-check -fbacktrace -fopenmp -J src/
test_SOURCES = $(corefiles) $(testfiles) $(auxfiles) $(cannfiles) $(nonadmdfiles) src/test_driver.f90


