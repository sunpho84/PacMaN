AC_INIT([template],[1.0],[fr.sanfilippo@gmail.com])

AC_CONFIG_SRCDIR([include])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.0 no-define foreign subdir-objects])
AM_CONFIG_HEADER(include/config.hpp)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

AC_PROG_RANLIB

AC_PROG_CC
AC_PROG_CXX
AC_LANG([C++])

#silent automake
AM_SILENT_RULES([yes])

LX_FIND_MPI
CPPFLAGS="$MPI_CXXFLAGS $CPPFLAGS"
LIBS="$MPI_CXXLIBS $LIBS"
LDFLAGS="$MPI_CXXLDFLAGS $LDFLAGS"

AX_CXX_COMPILE_STDCXX_11(noext,mandatory)
AX_CXXFLAGS_WARN_ALL

dnl AC_CHECK_LIB(gsl, gsl_sf_airy_Ai_e,
dnl              [LIBS="${LIBS} -lgsl -lgslcblas"], 
dnl              [AC_MSG_ERROR([Cannot find Gnu Scientific Library.])],
dnl              [-lgslcblas])
dnl AC_CHECK_HEADER([gsl/gsl_version.h], ,
dnl              [AC_MSG_ERROR([Cannot find GSL headers.])])

#openMP
AC_OPENMP
CFLAGS="$CXXFLAGS $OPENMP_CFLAGS"
CPPFLAGS="$CPPFLAGS $OPENMP_CPPFLAGS"
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"

CXXFLAGS="-O3 $CXXFLAGS"

AC_CONFIG_FILES(Makefile bin/Makefile)

AC_OUTPUT
