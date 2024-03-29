
# require at least autoconf 2.69
AC_PREREQ([2.69])

# Process this file with autoconf to produce a configure script.
AC_INIT([corels],[0.0.4])

# Ensure C++ is set up as R expects
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])
fi
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`
AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CC
AC_PROG_CXX

# Default to no GNU GMP
GMPFLAG=

# If gmp.h is found, define GMP
AC_CHECK_HEADER([gmp.h],
   [GMPFLAG=-DGMP],
   [AC_MSG_WARN([gmp.h not found, so skipping GMP resulting in slower code.])])

# If libgmp is found with __gmpz_init symbol, link to GMP
AC_CHECK_LIB(gmp, __gmpz_init, ,
   [AC_MSG_WARN([gmp library not found, so skipping GMP.])])

# Now substitute these variables in src/Makevars.in to create src/Makevars
AC_SUBST(GMP_CFLAGS, ["${GMPFLAG} ${CFLAGS}"])
AC_SUBST(GMP_LIBS, ${LIBS})
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
