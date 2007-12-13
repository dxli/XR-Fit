#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)
MY_SUBDIRS="cgi-bin"
AC_SUBST(MY_SUBDIRS)
AC_PROG_LIBTOOL
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL

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
AC_SEARCH_LIBS([pow, sqrt], [m])
AC_SEARCH_LIBS([gsl_integration_qags], [gsl])
AC_SEARCH_LIBS([GracePrintf], [grace_np],[],[echo "You have to install the grace_np devel files from the grace package"])
AC_CHECK_FUNCS([gettimeofday])
AC_PROG_CXX
AC_HAVE_LIBRARY(mysqlpp,[],[echo "\tmysql++ package not found"])
AC_CONFIG_FILES([Makefile
                 cgi-bin/Makefile])
AC_OUTPUT
