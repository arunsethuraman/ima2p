#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.63])
AC_INIT(IMa2p, 1.0,arun@temple.edu)
#PKG_PROG_PKG_CONFIG
#AC_CONFIG_MACRO_DIR([m4])
m4_include([m4/ax_prog_cxx_mpi.m4])
AM_INIT_AUTOMAKE(IMa2p, 1.0)
AC_OUTPUT(Makefile src/Makefile)
AC_CONFIG_SRCDIR([src/ima_main_mpi.cpp])
AC_CONFIG_HEADERS([src/config.h])
AM_MAINTAINER_MODE


# Checks for programs.
#AC_PROG_CXX
#AC_OPENMP
AC_PROG_AWK
#AC_PROG_CC
#AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB

#AC_CONFIG_MACRO_DIR([m4])
#m4_include([m4/ax_prog_cxx_mpi.m4])

AC_TYPE_SIZE_T
AC_C_VOLATILE
AC_C_INLINE
AC_C_CHAR_UNSIGNED
AC_FUNC_VPRINTF
#AC_DEFINE([MPI_ENABLED],[1],[MPI enabling])



AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
         [compile with MPI (parallelization) support. If none is found,
         MPI is not used. Default: auto])
     ],,[with_mpi=auto])

AX_PROG_CXX_MPI([test x"$with_mpi" != xno],[use_mpi=yes AC_DEFINE([MPI_ENABLED],[],[MPI enabling])],[
       use_mpi=no
       if test x"$with_mpi" = xyes; then
         AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
       else
         AC_MSG_WARN([No MPI compiler found, won't use MPI.])
       fi
	])






#MPI_CXXLIBS='mpicxx --showme:link'
#MPI_CXXFLAGS='mpicxx --showme:compile'
#AC_SUBST(MPI_CXXLIBS)
#AC_SUBST(MPI_CXXFLAGS)

# Checks for libraries.
AC_CHECK_LIB([m], [main])
# Checks for header files.
AC_CHECK_HEADERS([arpa/inet.h float.h limits.h netdb.h netinet/in.h stdlib.h string.h sys/socket.h unistd.h mpi.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_HEADER_STDBOOL

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor memset pow socket sqrt strchr strdup strrchr strstr strtol strtoul])


#AC_CONFIG_FILES([src/Makefile Makefile])
AC_OUTPUT
