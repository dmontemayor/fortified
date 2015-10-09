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
#AC_CONFIG_SRCDIR([src/test_driver.f90])
AC_CONFIG_SRCDIR([src/testdriver1.in])
AC_CONFIG_SRCDIR([src/testdriver2.in])
AC_CONFIG_SRCDIR([src/testdriver3.in])
AC_CONFIG_SRCDIR([src/testdriver4.in])
AC_CONFIG_SRCDIR([src/testdrivercap.in])

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
AC_CONFIG_COMMANDS([cp_datasets],[cp $ac_top_srcdir/data/doublewelldataset.dat $ac_top_srcdir/data/zpotdataset.dat $ac_top_srcdir/data/protstructdataset.dat ./data/])

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
