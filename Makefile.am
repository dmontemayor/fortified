ACLOCAL_AMFLAGS = -I m4 --install
EXTRA_DIST = mr/NOTES

#SUBDIRS = src
#EXTRA_DIST = etc libexec share Doxyfile DoxygenLayout.xml
#EXTRA_DIST = subsystem1.in ... subsystem9.in subsystemcap.in coupling3.in 

corefiles = src/type_kinds.f90 src/math.f90 src/rand.f90
testfiles = src/testing.f90 src/coretests.f90
classfiles = src/TEMPLATE.f90

bin_PROGRAMS = App test

App_FCFLAGS = -O3 -I $(EXEDIR)
test_FCFLAGS = -O3 -I $(EXEDIR)

App_SOURCES = $(corefiles) $(classfiles) src/driver.f90
test_SOURCES = $(corefiles) $(testfiles) $(classfiles) src/test_driver.f90

