#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS([stdlib.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_HEADER_TIME

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_FORK
AC_FUNC_STRTOD
AC_HAVE_LIBRARY([gsl],[],[echo "couldn't link -lgsl"])
AC_HAVE_LIBRARY([gslcblas],[],[echo "couldn't link -lgsl"])
AC_HAVE_LIBRARY([grace_np],[],[echo "couldn't link -lgsl"])
AC_HAVE_LIBRARY([mysqlpp],[],[echo "couldn't link -lmysqlpp"])
AC_CHECK_FUNCS([gettimeofday])

AC_CONFIG_FILES([Makefile
                 cgi-bin/Makefile
                 htdocs/Makefile])
AC_OUTPUT
