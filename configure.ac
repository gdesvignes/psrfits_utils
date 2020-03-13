AC_INIT([psrfits_utils], [1.0])
AC_CONFIG_SRCDIR([downsample.c])
AC_CONFIG_MACRO_DIR([config])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

m4_pattern_allow([AM_PROG_AR])
AM_PROG_AR
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CC_STDC
AC_PROG_LIBTOOL

PKG_CHECK_MODULES([GLIB], glib-2.0 >= 2.16)
#
AC_MSG_NOTICE([Adding glib-2.0 CFLAGS ${GLIB_CFLAGS}])
AM_CFLAGS="${GLIB_CFLAGS} ${AM_CFLAGS}"
AC_MSG_NOTICE([Adding glib-2.0 LIBS ${GLIB_LIBS}])
LIBS="${GLIB_LIBS} ${LIBS}"

AX_PTHREAD

SWIN_LIB_PSRCHIVE
SWIN_LIB_PGPLOT
SWIN_LIB_CFITSIO
SWIN_LIB_FFTW
SWIN_LIB_PRESTO


AC_CHECK_LIB([m],[sin])

AC_ARG_ENABLE([sse],
	      [AS_HELP_STRING([--enable-sse], 
			      [enable SSE folding code])],
	      [],
	      [enable_sse=no])
if test "x$enable_sse" = "xyes"; then
	AC_DEFINE([FOLD_USE_INTRINSICS], [1], [Use SSE intrinsics for folding])
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
