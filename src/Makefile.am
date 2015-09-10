bin_PROGRAMS = App test

App_FCFLAGS = -O3 -I $(EXEDIR)
test_FCFLACS = -03 -I $(EXEDIR)

corefiles = type_kinds.f90 math.f90 rand.f90
testfiles = testing.f90 coretests.f90
classfiles = TEMPLATE.f90

App_SOURCES = $(corefiles) $(classfiles) driver.f90
test_SOURCES = $(corefiles) $(testfiles) $(classfiles) test_driver.f90

#EXTRA_DIST = subsystem1.in subsystem2.in subsystem3.in subsystem4.in subsystem5.in subsystem6.in subsystem7.in subsystem8.in subsystem9.in subsystemcap.in coupling3.in TEMPLATE.f90

