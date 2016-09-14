#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# minimum autoconf version number  
AC_PREREQ([2.63])

# Initialize autoconf
AC_INIT([fortified], [0.0.0], [dmontemayor@qc.cuny.edu], [fortified-gnu])

# Check for source files
AC_CONFIG_SRCDIR([src/string.f90])
AC_CONFIG_SRCDIR([src/verification])
AC_CONFIG_SRCDIR([src/filemanager.f90])
AC_CONFIG_SRCDIR([src/testing.f90])
AC_CONFIG_SRCDIR([src/type_kinds.f90])
AC_CONFIG_SRCDIR([src/math.f90])
AC_CONFIG_SRCDIR([src/rand.f90])
AC_CONFIG_SRCDIR([src/template.f90])
AC_CONFIG_SRCDIR([src/functions.f90])
AC_CONFIG_SRCDIR([src/progressbar.f90])
AC_CONFIG_SRCDIR([src/testdriver1.in])
AC_CONFIG_SRCDIR([src/testdriver2.in])
AC_CONFIG_SRCDIR([src/testdriver3.in])
AC_CONFIG_SRCDIR([src/testdriver4.in])
AC_CONFIG_SRCDIR([src/testdrivercap.in])

# auxilary files
AC_CONFIG_SRCDIR([src/MDutils.f90])
AC_CONFIG_SRCDIR([src/atomicunits.f90])
AC_CONFIG_SRCDIR([src/lapack.f])
AC_CONFIG_SRCDIR([src/csv_io.f90])

# Setup where configuration results go so we can #include them easily 
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

# Setup automake
AM_INIT_AUTOMAKE([subdir-objects -Wall -Werror])
AC_CONFIG_FILES(Makefile)

# Checks for programs.
AC_PROG_F77
AC_PROG_FC

# make sure src and data directories are present in build directory
AC_CONFIG_COMMANDS([create_user_srcdir],[mkdir ./src/])
AC_CONFIG_COMMANDS([create_user_datadir],[mkdir ./data/])

################### F O R T I F I E D   M O D U L E S ###################

###### C A N N   m o d u l e ######
# cann classes
AC_CONFIG_SRCDIR([src/layer.f90])
AC_CONFIG_SRCDIR([src/ffn.f90])
AC_CONFIG_SRCDIR([src/som.f90])
AC_CONFIG_SRCDIR([src/device.f90])
AC_CONFIG_SRCDIR([src/rnn.f90])
AC_CONFIG_SRCDIR([src/bpg.f90])
 
# cann data files
AC_CONFIG_SRCDIR([data/gendoublewelldataset.f90])
AC_CONFIG_SRCDIR([data/doublewell.dat])
AC_CONFIG_SRCDIR([data/genzpotdataset.f90])
AC_CONFIG_SRCDIR([data/zpot.dat])
AC_CONFIG_SRCDIR([data/zpot.plt])
AC_CONFIG_SRCDIR([data/gencubicdataset.f90])
AC_CONFIG_SRCDIR([data/cubic.dat])
AC_CONFIG_SRCDIR([data/cubic.plt])
AC_CONFIG_SRCDIR([data/genxordataset.f90])
AC_CONFIG_SRCDIR([data/xor.dat])
AC_CONFIG_SRCDIR([data/xor.plt])
AC_CONFIG_SRCDIR([data/xorsinusoid.plt])
AC_CONFIG_SRCDIR([data/genxorsequentialdataset.f90])
AC_CONFIG_SRCDIR([data/xor.seqdat])
AC_CONFIG_SRCDIR([data/xorseq.plt])
AC_CONFIG_SRCDIR([data/protstruct.seqdat])
AC_CONFIG_SRCDIR([data/protstruct_test.seqdat])

# make local copy of datasets and scripts 
AC_CONFIG_COMMANDS([cp_datasets],[cp $ac_top_srcdir/data/*.dat $ac_top_srcdir/data/*.seqdat $ac_top_srcdir/data/*.plt ./data/])

# create unique drivers in build directory
AC_CONFIG_COMMANDS([create_cann_driver],[cp $ac_top_srcdir/src/cann_driver.f90 src/cann_driver.f90])

###### N O N A D M D   m o d u l e ######
# nonadmd elements
AC_CONFIG_SRCDIR([src/quasiparticle.f90])
AC_CONFIG_SRCDIR([src/ham.f90])
AC_CONFIG_SRCDIR([src/wavelet.f90])
AC_CONFIG_SRCDIR([src/fauxnon.f90])
 
# nonadmd models
AC_CONFIG_SRCDIR([src/metiu.f90])
AC_CONFIG_SRCDIR([src/cukier.f90])

# nonadmd engines
AC_CONFIG_SRCDIR([src/GQFPE.f90])

# nonadmd data files
#AC_CONFIG_SRCDIR([data/XXXX.f90])

# make local copy of datasets and scripts 
#AC_CONFIG_COMMANDS([cp_datasets],[cp $ac_top_srcdir/data/*.dat $ac_top_srcdir/data/*.seqdat $ac_top_srcdir/data/*.plt ./data/])

# create unique drivers in build directory
AC_CONFIG_COMMANDS([create_nonadmd_driver],[cp $ac_top_srcdir/src/nonadmd_driver.f90 src/nonadmd_driver.f90])


######  T E S T   D R I V E R   ######

# create test driver
AC_CONFIG_COMMANDS([create_test_driver],[cp $ac_top_srcdir/src/testdriver1.in $ac_top_srcdir/src/utest.in
# include mods in test driver
$ac_top_srcdir/libexec/sandwich '  use ' $ac_top_srcdir/cann.mods '_class' >> $ac_top_srcdir/src/utest.in
$ac_top_srcdir/libexec/sandwich '  use ' $ac_top_srcdir/nonadmd.mods '_class' >> $ac_top_srcdir/src/utest.in
# add help argument dialogue for test driver
cat $ac_top_srcdir/src/testdriver2.in >> $ac_top_srcdir/src/utest.in
$ac_top_srcdir/libexec/sandwich "     write(*,*)'" $ac_top_srcdir/cann.mods "'" >> $ac_top_srcdir/src/utest.in
$ac_top_srcdir/libexec/sandwich "     write(*,*)'" $ac_top_srcdir/nonadmd.mods "'" >> $ac_top_srcdir/src/utest.in
# add all tests option to test driver
cat $ac_top_srcdir/src/testdriver3.in >> $ac_top_srcdir/src/utest.in
$ac_top_srcdir/libexec/sandwich '     call ' $ac_top_srcdir/cann.mods '_test' >> $ac_top_srcdir/src/utest.in
$ac_top_srcdir/libexec/sandwich '     call ' $ac_top_srcdir/nonadmd.mods '_test' >> $ac_top_srcdir/src/utest.in
# add cann module test option to test driver
echo "  case('cann')" >> $ac_top_srcdir/src/utest.in
$ac_top_srcdir/libexec/sandwich '     call ' $ac_top_srcdir/cann.mods '_test' >> $ac_top_srcdir/src/utest.in
# add nonadmd module test option to test driver
echo "  case('nonadmd')" >> $ac_top_srcdir/src/utest.in
$ac_top_srcdir/libexec/sandwich '     call ' $ac_top_srcdir/nonadmd.mods '_test' >> $ac_top_srcdir/src/utest.in
# add individual class test options to test driver
cat $ac_top_srcdir/src/testdriver4.in >> $ac_top_srcdir/src/utest.in
$ac_top_srcdir/libexec/sandwich2 "  case('" $ac_top_srcdir/cann.mods "'); call " '_test' >> $ac_top_srcdir/src/utest.in
$ac_top_srcdir/libexec/sandwich2 "  case('" $ac_top_srcdir/nonadmd.mods "'); call " '_test' >> $ac_top_srcdir/src/utest.in
#complete test driver
cat $ac_top_srcdir/src/testdrivercap.in >> $ac_top_srcdir/src/utest.in
mv $ac_top_srcdir/src/utest.in ./src/test_driver.f90])


AC_OUTPUT
