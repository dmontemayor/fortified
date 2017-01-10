ACLOCAL_AMFLAGS = -I m4 --install
EXTRA_DIST = m4 data testing

corefiles = src/type_kinds.f90 src/math.f90 src/rand.f90 src/testing.f90 src/functions.f90 src/progressbar.f90 src/filemanager.f90 src/string.f90 src/template.f90
testfiles = src/coretests.f90

auxfiles = src/MDutils.f90 src/atomicunits.f90 src/csv_io.f90 src/molreader.f90

lapack =  lapack-3.6.1/libtmglib.a lapack-3.6.1/liblapack.a lapack-3.6.1/librefblas.a 

# All fortified executable modules
bin_PROGRAMS = nonadmd cann test ontologyparser

############################
###### ONTOLOGYPARSER ###### 	
ontologyparser_FCFLAGS = -g -fbounds-check -fbacktrace -fopenmp -J src/
ontologyfiles = src/ontologyparser.f90 testing/test_ontologyparser.f90
ontologyparser_SOURCES = $(corefiles) $(ontologyfiles)

#######################################################
###### CANN chemistry artificial neural networks ###### 	
cann_FCFLAGS =  -fopenmp -J src/
cannfiles = src/layer.f90 src/bpg.f90 src/ffn.f90 src/som.f90 src/device.f90
cann_SOURCES = $(corefiles) $(auxfiles) $(cannfiles) src/cann_driver.f90
cann_LDADD = $(lapack)
######################################################
###### NONADMD Non-Adiabatic Molecular Dynamics ###### 	
nonadmd_FCFLAGS =  -fopenmp -J src/
elements = src/wf.f90 src/quasiparticle.f90 src/ham.f90 src/wavelet.f90 src/fauxnon.f90 #src/lapack.f
models= #src/metiu.f90 src/cukier.f90 
engines= #src/GQFPE.f90
nonadmd_SOURCES = $(corefiles) $(auxfiles) $(elements) $(models) $(engines) src/nonadmd_driver.f90
nonadmd_LDADD = $(lapack)

#########################
###### TEST module ######
test_FCFLAGS = -g -fbounds-check -fbacktrace -fopenmp -J src/
test_SOURCES = $(corefiles) $(testfiles) $(auxfiles) $(cannfiles) $(elements) $(models) $(engines) src/test_driver.f90
test_LDADD = $(lapack)

