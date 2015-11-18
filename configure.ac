#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# minimum autoconf version number  
AC_PREREQ([2.63])

# Initialize autoconf
AC_INIT([fortified], [0.0.0], [dmontemayor@qc.cuny.edu], [fortified-gnu])

# Check for source files
AC_CONFIG_SRCDIR([src/testing.f90])
AC_CONFIG_SRCDIR([src/type_kinds.f90])
AC_CONFIG_SRCDIR([src/math.f90])
AC_CONFIG_SRCDIR([src/rand.f90])
AC_CONFIG_SRCDIR([src/TEMPLATE.f90])
AC_CONFIG_SRCDIR([src/driver.f90])
AC_CONFIG_SRCDIR([src/functions.f90])
AC_CONFIG_SRCDIR([src/progressbar.f90])
AC_CONFIG_SRCDIR([src/testdriver1.in])
AC_CONFIG_SRCDIR([src/testdriver2.in])
AC_CONFIG_SRCDIR([src/testdriver3.in])
AC_CONFIG_SRCDIR([src/testdriver4.in])
AC_CONFIG_SRCDIR([src/testdrivercap.in])

# derived classes
AC_CONFIG_SRCDIR([src/layer.f90])
AC_CONFIG_SRCDIR([src/ffn.f90])
AC_CONFIG_SRCDIR([src/som.f90])
AC_CONFIG_SRCDIR([src/device.f90])
AC_CONFIG_SRCDIR([src/rnn.f90])
AC_CONFIG_SRCDIR([src/bpg.f90])

# auxilary modules
AC_CONFIG_SRCDIR([src/MDutils.f90])

# auxilary data files
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


# Setup where configuration results go so we can #include them easily 
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

# Setup automake
AM_INIT_AUTOMAKE([subdir-objects -Wall -Werror])
AC_CONFIG_FILES(Makefile)

# Checks for programs.
AC_PROG_FC
#AC_CHECK_PROGS([PWD],[pwd],[:])
#if test "$PWD" = :; then
#   AC_MSG_ERROR([This package needs pwd.])
#fi
#AC_DEFINE_UNQUOTED([EXEDIR],["$PWD"],[Build directory.])
#AC_SUBST([EXEDIR],[$PWD])

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

# make sure src/ is present in build directory
AC_CONFIG_COMMANDS([create_user_srcdir],[mkdir ./src/])

# make sure src/ is present in build directory
AC_CONFIG_COMMANDS([create_user_datadir],[mkdir ./data/])
AC_CONFIG_COMMANDS([cp_datasets],[cp $ac_top_srcdir/data/*.dat $ac_top_srcdir/data/*.seqdat $ac_top_srcdir/data/*.plt ./data/])

# create unique driver in build directory
AC_CONFIG_COMMANDS([create_system_driver],[cp $ac_top_srcdir/src/driver.f90 src/driver.f90])

# find all unit tests and create a test driver
AC_CONFIG_COMMANDS([create_test_driver],[cp $ac_top_srcdir/src/testdriver1.in $ac_top_srcdir/src/utest.in
$ac_top_srcdir/libexec/sandwich '  use ' $ac_top_srcdir/MODS '_class' >> $ac_top_srcdir/src/utest.in
cat $ac_top_srcdir/src/testdriver2.in >> $ac_top_srcdir/src/utest.in
$ac_top_srcdir/libexec/sandwich "     write(*,*)'" $ac_top_srcdir/MODS "'" >> $ac_top_srcdir/src/utest.in
cat $ac_top_srcdir/src/testdriver3.in >> $ac_top_srcdir/src/utest.in
$ac_top_srcdir/libexec/sandwich '     call ' $ac_top_srcdir/MODS '_test' >> $ac_top_srcdir/src/utest.in
cat $ac_top_srcdir/src/testdriver4.in >> $ac_top_srcdir/src/utest.in
$ac_top_srcdir/libexec/sandwich2 "  case('" $ac_top_srcdir/MODS "'); call " '_test' >> $ac_top_srcdir/src/utest.in
cat $ac_top_srcdir/src/testdrivercap.in >> $ac_top_srcdir/src/utest.in
mv $ac_top_srcdir/src/utest.in ./src/test_driver.f90])

AC_OUTPUT
