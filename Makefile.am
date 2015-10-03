ACLOCAL_AMFLAGS = -I m4 --install
EXTRA_DIST = m4/NOTES data/NOTES data/gencubicdataset.f90 data/cubicdataset.dat

#SUBDIRS = src
#EXTRA_DIST = etc libexec share Doxyfile DoxygenLayout.xml
#EXTRA_DIST = subsystem1.in ... subsystem9.in subsystemcap.in coupling3.in 

corefiles = src/type_kinds.f90 src/math.f90 src/rand.f90 src/testing.f90
classfiles = src/TEMPLATE.f90 src/layer.f90
testfiles = src/coretests.f90

bin_PROGRAMS = App test

App_FCFLAGS = -O3 -J src/
test_FCFLAGS = -O3 -J src/

App_SOURCES = $(corefiles) $(classfiles) src/driver.f90
test_SOURCES = $(corefiles) $(testfiles) $(classfiles) src/test_driver.f90

