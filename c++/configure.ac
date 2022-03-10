#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT(pr-calc, 1.0, malvarez@lbl.gov)
AM_INIT_AUTOMAKE(pr-calc, 1.0)
AC_CONFIG_SRCDIR([src]) 

# Checks for programs.
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h mach/mach.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

AC_OUTPUT(Makefile src/Makefile)
