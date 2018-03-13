dnl @synopsis SWIN_LIB_PRESTO
dnl 
AC_DEFUN([SWIN_LIB_PRESTO],
[
  AC_PROVIDE([SWIN_LIB_PRESTO])

  SWIN_PACKAGE_OPTIONS([presto])

  AC_MSG_CHECKING([for PRESTO functions])

  PRESTO_CFLAGS="-I$PRESTO/include"
  PRESTO_LIBS="-L$PRESTO/lib -lpresto"

  ac_save_CFLAGS="$CFLAGS"
  ac_save_LIBS="$LIBS"

  LIBS="$ac_save_LIBS $PRESTO_LIBS"
  CFLAGS="$ac_save_CFLAGS $PRESTO_CFLAGS"

  SWIN_PACKAGE_FIND([presto],[presto.h])
  SWIN_PACKAGE_TRY_COMPILE([presto],[#include <presto.h>])

  if test $have_presto = yes; then
    SWIN_PACKAGE_FIND([presto],[libpresto.*])
    SWIN_PACKAGE_TRY_LINK([presto],[#include <presto.h>],
                            [delay_from_dm(0,1400.0);],
                            [-lpresto])
  fi


  AC_MSG_RESULT([$have_presto])

  LIBS="$ac_save_LIBS"
  CFLAGS="$ac_save_CFLAGS"


  if test x$have_presto = xyes; then
    AC_DEFINE([HAVE_PRESTO], [1],[Define to 1 if the PRESTO library is installed])
    [$1]
  fi

  AC_SUBST(PRESTO_LIBS)
  AC_SUBST(PRESTO_CFLAGS)

  AM_CONDITIONAL(HAVE_PRESTO,[test x"$have_presto" = xyes])

])

