#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.64])
AC_INIT([fastsg], [0.0], [murarasu@in.tum.de])

AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE(fastsg, 0.0)
LT_INIT

AC_CONFIG_SRCDIR([src/Helper.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset])

AC_CONFIG_FILES([Makefile src/Makefile examples/Makefile])

AC_OUTPUT
