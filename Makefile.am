ACLOCAL_AMFLAGS = -I m4 --install
EXTRA_DIST = m4 data

#EXTRA_DIST = m4/NOTES data/NOTES data/gendoublewelldataset.f90 data/doublewell.dat data/genzpotdataset.f90 data/zpot.dat data/zpot.plt data/gencubicdataset.f90 data/cubic.dat data/cubic.plt data/genxordataset.f90 data/xor.dat data/xor.plt data/genxorsequentialdataset.f90 data/xor.seqdat data/xorseq.plt data/protstruct.dat data/protstruct_test.dat

corefiles = src/type_kinds.f90 src/math.f90 src/rand.f90 src/testing.f90 src/functions.f90 src/progressbar.f90

classfiles = src/TEMPLATE.f90 src/layer.f90 src/ffn.f90 src/som.f90 src/device.f90 src/rnn.f90

testfiles = src/coretests.f90

bin_PROGRAMS = App test

App_FCFLAGS = -J src/
test_FCFLAGS = -fbounds-check -fbacktrace -J src/

App_SOURCES = $(corefiles) $(classfiles) src/driver.f90
test_SOURCES = $(corefiles) $(testfiles) $(classfiles) src/test_driver.f90

