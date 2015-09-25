#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# minimum autoconf version number  
AC_PREREQ([2.68])

# Initialize autoconf
AC_INIT([fortified], [0.0.0], [dmontemayor@qc.cuny.edu], [fortified-gnu])

# Check for source files
AC_CONFIG_SRCDIR([src/testing.f90])
AC_CONFIG_SRCDIR([src/type_kinds.f90])
AC_CONFIG_SRCDIR([src/math.f90])
AC_CONFIG_SRCDIR([src/rand.f90])
AC_CONFIG_SRCDIR([src/TEMPLATE.f90])
AC_CONFIG_SRCDIR([src/driver.f90])
AC_CONFIG_SRCDIR([src/test_driver.f90])

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

AC_OUTPUT
