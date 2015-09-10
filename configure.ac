#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([App], [0.0.0], [dmontemayor@qc.cuny.edu])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_CHECK_PROGS([PWD],[pwd],[:])
if test "$PWD" = :; then
   AC_MSG_ERROR([This package needs pwd.])
fi
AC_DEFINE_UNQUOTED([EXEDIR],["$PWD"],[Build directory.])
AC_SUBST([EXEDIR],[$PWD])

# Check for source files
AC_CONFIG_SRCDIR([src/testing.f90])
AC_CONFIG_SRCDIR([src/type_kinds.f90])
AC_CONFIG_SRCDIR([src/math.f90])
AC_CONFIG_SRCDIR([src/rand.f90])
AC_CONFIG_SRCDIR([src/TEMPLATE.f90])
AC_CONFIG_SRCDIR([src/driver.f90])
AC_CONFIG_SRCDIR([src/test_driver.f90])


# Checks for libraries.

# Checks for header files.


# Checks for typedefs, structures, and compiler characteristics.
AC_PROG_FC

# Checks for library functions.
AC_CONFIG_FILES([Makefile src/Makefile])

AC_OUTPUT
