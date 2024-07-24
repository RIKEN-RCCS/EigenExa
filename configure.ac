#--------*---------*---------*---------*---------*---------*---------*-*
#
#  File     : configure.ac
#
# Copyright(C) 2012-2024 RIKEN.
# Copyright(C) 2011-2012 Toshiyuki Imamura
#                        Graduate School of Informatics and Engineering,
#                        The University of Electro-Communications.
# Copyright (C) 2011- 2015 Japan Atomic Energy Agency.
# 
# Redistribution  and  use  in  source and binary forms, with or without
# modification,  are  permitted  provided  that the following conditions
# are met:
# 
# * Redistributions  of  source  code  must  retain  the above copyright
#   notice,  this  list  of  conditions  and  the  following  disclaimer.
# * Redistributions  in  binary  form must reproduce the above copyright
#   notice,  this list of conditions and the following disclaimer in the
#   documentation  and/or other materials provided with the distribution.
# 
# THIS  SOFTWARE  IS  PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# ``AS IS''  AND  ANY  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A  PARTICULAR  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL  DAMAGES  (INCLUDING,  BUT NOT
# LIMITED  TO,  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA,  OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY  OF  LIABILITY,  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF  THIS  SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#--------*---------*---------*---------*---------*---------*---------*-*

AC_PREREQ(2.59)
AC_INIT(EigenExa, [2.13], EigenExa@ml.riken.jp, EigenExa, https://https://github.com/RIKEN-RCCS/EigenExa)
AC_CONFIG_SRCDIR([src/eigen_libs.F])
AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE([foreign])

AC_PREFIX_DEFAULT([$PWD])


## make set
AC_PROG_MAKE_SET


## ranlib
AC_PROG_RANLIB

## check MacOS
if test `uname` = Darwin; then 
   OS_MACOS=yes
   SUFFIX_DYLIB=dylib
else
   SUFFIX_DYLIB=so
fi
AC_SUBST(SUFFIX_DYLIB)
AC_SUBST(OS_MACOS)

## shared lib
AC_ARG_ENABLE([shared-library],
    [AS_HELP_STRING(
      [--enable-shared-library],
      [enable Shared library (default).])],
    [
      if test x"${enableval}" = x"yes"; then
	    USE_SHARED_LIB=1
      else
	    USE_SHARED_LIB=0
      fi
    ],
	[
      USE_SHARED_LIB=1
      ]
  )
AC_SUBST(USE_SHARED_LIB)


## Timer
#

AC_ARG_ENABLE([timer],
    [AS_HELP_STRING(
      [--enable-timer@<:@=ARG@:>@],
      [set Timer Print Level (from 0 to 2, default = 0)])],
    [
      if test x"${enableval}" = x"yes"; then
	    OPT_TIMER_PRINT=0
      elif test x"${enableval}" = x"0" -o \
                x"${enableval}" = x"1" -o \
                x"${enableval}" = x"2"; then
        OPT_TIMER_PRINT=${enableval}
      else
        AC_MSG_ERROR([Not supported Timer Print Level ${enableval}.])

      fi
    ],
	[
      OPT_TIMER_PRINT=0
    ])

AC_SUBST(OPT_TIMER_PRINT)


## Cache alignment
#
AC_ARG_ENABLE([cache-alignment],
    [AS_HELP_STRING(
      [--enable-cache-alignment],
      [Boost by Cache alignment (default= no)])],
    [
      if test x"${enableval}" = x"yes"; then
	    OPT_CACHE_ALIGN=1
      else
	    OPT_CACHE_ALIGN=0
      fi
    ],
	[
      OPT_CACHE_ALIGN=0
    ])

AC_SUBST(OPT_CACHE_ALIGN)



AC_ARG_ENABLE([mpi],
    [AS_HELP_STRING(
      [--enable-mpi={intel,mpich,openmpi}],
      [MPI library.])],
    [
      if test x"${enableval}" = x"intel"; then
        if test ${CC_ACT} != icc; then
	    CC=mpiicc
	fi
        if test ${CC_ACT} != icx; then
	    CC=mpiicx
	fi
	    F77=mpiifort
	    FC=mpiifort
	    echo 'Intel MPI and MKL are forced to use'
	    LAPACK_LIBS=
	    LAPACK_PATH=
      else
	    CC=mpicc
	    F77=mpif90
	    FC=mpif90
      fi
    ],
    []
  )

AC_ARG_WITH([mkl],
    [AS_HELP_STRING(
      [--with-mkl],
      [MKL (BLAS and LAPACK) library (only x86).])],
    []
    [with_mkl=yes]
  )

AC_ARG_WITH([mkl-scalapack],
    [AS_HELP_STRING(
      [--with-mkl-scalapack],
      [MKL Scalapack library (only x86, link with Intel MPI)])],
    []
    [with_mkl=check; with_mkl_scalapack=check]
  )

## Intel CPU specific optimization flags
#
AC_ARG_ENABLE([avx],
    [AS_HELP_STRING(
      [--enable-avx],
      [enable AVX features if available])],
    [
      if test x"${enableval}" = x"yes"; then
	    OPT_AVX=1
      else
	    OPT_AVX=0
      fi
    ],
	[
      OPT_AVX=0
    ])
AC_SUBST(OPT_AVX)

AC_ARG_ENABLE([avx2],
    [AS_HELP_STRING(
      [--enable-avx2],
      [enable AVX2 features if available])],
    [
      if test x"${enableval}" = x"yes"; then
	    OPT_AVX2=1
      else
	    OPT_AVX2=0
      fi
    ],
	[
      OPT_AVX2=0
    ])
AC_SUBST(OPT_AVX2)

AC_ARG_ENABLE([mic-avx512],
    [AS_HELP_STRING(
      [--enable-mic-avx512],
      [enable MIC-AVX512 features if available])],
    [
      if test x"${enableval}" = x"yes"; then
	    OPT_MIC_AVX512=1
      else
	    OPT_MIC_AVX512=0
      fi
    ],
	[
      OPT_MIC_AVX512=0
    ])
AC_SUBST(OPT_MIC_AVX512)

AC_ARG_ENABLE([skylake-avx512],
    [AS_HELP_STRING(
      [--enable-skylake-avx512],
      [enable SKYLAKE-AVX512 features if available])],
    [
      if test x"${enableval}" = x"yes"; then
	    OPT_SKYLAKE_AVX512=1
      else
	    OPT_SKYLAKE_AVX512=0
      fi
    ],
	[
      OPT_SKYLAKE_AVX512=0
    ])
AC_SUBST(OPT_SKYLAKE_AVX512)


AC_ARG_ENABLE([debug],
    [AS_HELP_STRING(
      [--enable-debug],
      [enable DEBUG macro(DEBUG=1)])],
    [
      if test x"${enableval}" = x"yes"; then
	    DEBUG=1
      else
	    DEBUG=0
      fi
    ],
	[
      DEBUG=0
    ])
AC_SUBST(DEBUG)


## fortran compiler
#

eigenexa_default_fc="mpifrtpx mpiifort mpixlf90_r sxmpif90 mpif90"

# check usablility
AC_PROG_F77($eigenexa_default_fc)
AC_PROG_FC($eigenexa_default_fc)

if test x"${F77}" != x -a x"${FC}" = x; then
  FC=${F77}
fi
if test x"${F77}" = x -a x"${FC}" != x; then
  F77=${FC}
fi
if test x"${F77}" != x -a x"${FC}" != x; then
  FC=${F77}
fi

# check actual compiler
AC_MSG_CHECKING([for actual compiler])

fc_ver_info=$(${F77} -E --version)
for vtok in ${fc_ver_info} ; do

  if test x"${vtok}" = x"ifort"; then
    FC_ACT="ifort"
    break

  elif test x"${vtok}" = x"GNU"; then
    FC_ACT="gfortran"
    break

  fi

done

if test -z "${FC_ACT}"; then
  FC_ACT=${F77}
fi

if test x"${FC_ACT}" = x"ifort"; then
  fc_ver_info=$(${F77} -E --version)
  ICC_VERSION=$(echo ${fc_ver_info} | awk '{ print $3}')
  ICC_MAJOR_VERSION=${ICC_VERSION%%.*}
  if test ${ICC_MAJOR_VERSION} -ge 15 ; then
    OMP="-qopenmp"
  else
    OMP="-openmp"
  fi

  AC_MSG_RESULT([$FC_ACT $ICC_VERSION])

else

  AC_MSG_RESULT([$FC_ACT])

fi


## C compiler
#

eigenexa_default_cc="mpifccpx mpiicx mpiicc mpixlc_r sxmpic++ mpicc"

# check usablility
AC_PROG_CC($eigenexa_default_cc)

# check actual compiler
AC_MSG_CHECKING([for actual compiler])

cc_ver_info=$(${CC} -E --version)
for vtok in ${cc_ver_info} ; do

  if test x"${vtok}" = x"oneAPI"; then
    CC_ACT="icx"
    break

  elif test x"${vtok}" = x"icc"; then
    CC_ACT="icc"
    break

  elif test x"${vtok}" = x"GNU"; then
    CC_ACT="gcc"
    break

  fi

done

if test -z "${CC_ACT}"; then
  CC_ACT=${CC}
fi

if test ${CC_ACT} = icc -o ${CC_ACT} = icx; then
if test ${FC_ACT} != ifort; then
      AC_MSG_ERROR([C and Fortran compilers must be the same seriese])
fi
fi
if test ${CC_ACT} != icc -a ${CC_ACT} != icx -a ${FC_ACT} = ifort; then
      AC_MSG_ERROR([C and Fortran compilers must be the same seriese])
fi
if test ${CC} = mpiicc -o ${CC} = mpiicx; then
if test ${F77} != mpiifort; then
      AC_MSG_ERROR([C and Fortran compilers must be the same seriese])
fi
fi
if test ${CC} != mpiicc -a ${CC} != mpiicx -a ${F77} = mpiifort; then
      AC_MSG_ERROR([C and Fortran compilers must be the same seriese])
fi

AC_MSG_RESULT([$CC_ACT])


## lapack
#

AC_ARG_VAR([LAPACK_PATH], 
    [LAPACK library path (linker flags, e.g. -L<lapack dir>)])

AC_ARG_VAR([LAPACK_LIBS], 
    [LAPACK library files (libraries, e.g. -l<lapack library>)])


#
if test -n "${LAPACK_PATH}"; then
  for lib in ${LAPACK_PATH}
  do
    case ${lib} in
    -L*) OPT_LD_LAPACK+="${lib} " ;;
    *) OPT_LD_LAPACK+="-L${lib} " ;;
    esac
  done
fi

#
if test x"${with_mkl_scalapack}" != x; then
  if test x"${with_mkl}" = x; then
      AC_MSG_ERROR([MKL Scalapack needs Intel MKL])
  fi
  if test x"${F77}" != x"mpiifort"; then
      AC_MSG_ERROR([MKL Scalapack should be linked with Intel MPI/Intel compiler])
  fi
fi

USE_MKL=0
if test -n "${LAPACK_LIBS}"; then
  OPT_LD_LAPACK+=" ${LAPACK_LIBS} "

 if test x"${FC_ACT}" = x"ifort" -o x"${FC_ACT}" = x"gfortran"; then
  if test x"${with_mkl}" != x -o x"${with_mkl_scalapack}" != x; then
    if test x"${FC_ACT}" = x"ifort"; then
      OPT_LD_LAPACK+=" "
    fi
    if test x"${OS_MACOS}" = x"yes"; then
      OPT_LD_LAPACK+=" -I${MKLROOT}/include -L${MKLROOT}/lib "
    else
      OPT_LD_LAPACK+=" -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 "
    fi
  fi
  if test x${F77} = xmpiifort -a  x"${with_mkl_scalapack}" != x; then
    OPT_LD_LAPACK+=" -lmkl_scalapack_lp64 "
  fi
  if test x"${with_mkl}" != x -o x"${with_mkl_scalapack}" != x; then
    if test x"${OS_MACOS}" = x"yes"; then
      OPT_LD_LAPACK+=" -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_mpich_lp64 -liomp5 -lpthread -lm -ldl "
    else
      OPT_LD_LAPACK+=" -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl "
    fi
    USE_MKL=1
  fi
 fi

elif test x"${FC_ACT}" = x"ifort"; then
  if test x"${OS_MACOS}" = x"yes"; then
    OPT_LD_LAPACK+=" -I${MKLROOT}/include -L${MKLROOT}/lib "
  else
    OPT_LD_LAPACK+=" -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 "
  fi
  if test x${F77} = xmpiifort -o  x"${with_mkl_scalapack}" = x; then
    OPT_LD_LAPACK+=" -lmkl_scalapack_lp64 "
  fi
  if test x"${OS_MACOS}" = x"yes"; then
    OPT_LD_LAPACK+=" -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_mpich_lp64 -liomp5 -lpthread -lm -ldl "
  else
    OPT_LD_LAPACK+=" -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl "
  fi
  USE_MKL=1

elif test x"${FC_ACT}" = x"gfortran"; then
  if test x"${with_mkl}" != x -o x"${with_mkl_scalapack}" != x; then
    OPT_LD_LAPACK+=" -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 "
  fi
  if test x"${with_mkl_scalapack}" = x; then
    OPT_LD_LAPACK+=" -lscalapack"
  fi
  if test x"${with_mkl}" != x -o x"${with_mkl_scalapack}" != x; then
    OPT_LD_LAPACK+=" -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl "
    USE_MKL=1
  else
    OPT_LD_LAPACK+=" -llapack -lblas "
  fi

elif test x"${FC_ACT}" = x"mpifrtpx"; then
  OPT_LD_LAPACK+=" -SCALAPACK -SSL2BLAMP "

elif test x"${FC_ACT}" = x"mpixlf90_r"; then
  OPT_LD_LAPACK+=" -L${SCALAPACK_LIB} -L${LAPACK_LIB} -L/bgsys/local/lib -lscalapack -llapack -lesslsmpbg "

elif test x"${FC_ACT}" = x"sxmpif90"; then
  OPT_LD_LAPACK+=" -lscalapack -lblacsF90init -lblacs -lblacsF90init -llapack -lparblas "

fi

AC_SUBST(USE_MKL)
AC_SUBST(OPT_LD_LAPACK)


## For intel compiler architecture-oriented options
if test x"${FC_ACT}" = x"ifort"; then
  I_ARCH=""
  if test ${OPT_AVX} -eq 1 ; then
    I_ARCH=" -axAVX "
  elif test ${OPT_AVX2} -eq 1 ; then
    I_ARCH=" -axAVX2 "
  elif test ${OPT_MIC_AVX512} -eq 1 ; then
    I_ARCH=" -axMIC-AVX512 "
  elif test ${OPT_SKYLAKE_AVX512} -eq 1 ; then
    I_ARCH=" -axSKYLAKE-AVX512 "
  else
#    I_ARCH=" -xHOST "
    I_ARCH=" "
  fi
  I_BASE=${OMP}" -g -fp-model strict -inline-level=1 "
fi

## For GNU compiler architecture-oriented options
if test x"${FC_ACT}" = x"gfortran"; then
  G_ARCH=""
  if test ${OPT_AVX} -eq 1 ; then
    G_ARCH=" -march=avx "
  elif test ${OPT_AVX2} -eq 1 ; then
    G_ARCH=" -march=avx2 "
  elif test ${OPT_MIC_AVX512} -eq 1 ; then
#    G_ARCH=" -march=mic-avx512 "
    G_ARCH=" -march=knl "
  elif test ${OPT_SKYLAKE_AVX512} -eq 1 ; then
    G_ARCH=" -march=skylake-avx512 "
  fi
  G_BASE=" -g -fopenmp "
fi

## For Fujitsu compiler architecture-oriented options
if test x"${FC_ACT}" = x"mpifrtpx"; then
  F_ARCH=""
  F_BASE=" -g -Kopenmp,openmp_ordered_reduction "
fi


## compile flags
#

# FFLAGS
if test x"${ac_test_FFLAGS}" = x"set"; then
  FFLAGS=${ac_save_FFLAGS}

else

  if   test x"${FC_ACT}" = x"ifort"; then
    FFLAGS=${I_BASE}" -Ofast "${I_ARCH}" -W1 -Winline"
    
  elif test x"${FC_ACT}" = x"gfortran"; then
    FFLAGS=${G_BASE}" -O3 "${G_ARCH}
    
  elif test x"${FC_ACT}" = x"mpifrtpx"; then
    FFLAGS=${F_BASE}" -Kfast,ocl,simd,swp "${F_ARCH}
    
  elif test x"${FC_ACT}" = x"mpixlf90_r"; then
    FFLAGS="-qfixed -qsmp=omp -qthreaded -O3 -qessl -qstrict -qarch=qp -qtune=qp -qsimd=auto -qassert=contig -qdirective"
  
  elif test x"${FC_ACT}" = x"sxmpif90"; then
    FFLAGS="-Popenmp -Cvsafe -Ep"
  
  fi

fi

# FCFLAGS
if test x"${ac_test_FCFLAGS}" = x"set"; then
  FCFLAGS=${ac_save_FCFLAGS}

else

  if   test x"${FC_ACT}" = x"ifort"; then
    FCFLAGS=${I_BASE}" -Ofast "${I_ARCH}" -W1 -Winline"
    
  elif test x"${FC_ACT}" = x"gfortran"; then
    FCFLAGS=${G_BASE}" -O3 "${G_ARCH}
    
  elif test x"${FC_ACT}" = x"mpifrtpx"; then
    FCFLAGS=${F_BASE}" -Kfast,ocl,simd,swp "${F_ARCH}
    
  elif test x"${FC_ACT}" = x"mpixlf90_r"; then
    FCFLAGS="-qfixed -qsmp=omp -qthreaded -O3 -qessl -qstrict -qarch=qp -qtune=qp -qsimd=auto -qassert=contig -qdirective"
  
  elif test x"${FC_ACT}" = x"sxmpif90"; then
    FCFLAGS="-Popenmp -Cvsafe -Ep"
  
  fi

fi

# FFLAGS0 (special flags)

if   test x"${FC_ACT}" = x"ifort"; then
  FFLAGS0=${I_BASE}" -O0 -W1 -Winline"
  
elif test x"${FC_ACT}" = x"gfortran"; then
  FFLAGS0=${G_BASE}" -O0"
  
elif test x"${FC_ACT}" = x"mpifrtpx"; then
  FFLAGS0="-Kopenmp,openmp_ordered_reduction -Cpp"
  
elif test x"${FC_ACT}" = x"mpixlf90_r"; then
  FFLAGS0="-qfixed -qsmp=omp -qthreaded -O0 -qessl -qstrict"

elif test x"${FC_ACT}" = x"sxmpif90"; then
  FFLAGS0="-Popenmp -Ep"

fi

AC_SUBST(FFLAGS0)


# CFLAGS
if test x"${ac_test_CFLAGS}" = x"set"; then
  CFLAGS=${ac_save_CFLAGS}

else

  if   test x"${FC_ACT}" = x"ifort"; then
    CFLAGS=${I_BASE}" -O3 "${I_ARCH}" -Wall"
    
  elif test x"${FC_ACT}" = x"gfortran"; then
    CFLAGS=${G_BASE}" -O3 "${G_ARCH}" -Wall"
    
  elif test x"${FC_ACT}" = x"mpifrtpx"; then
    CFLAGS=${F_BASE}" -Kfast,ocl,simd,swp "${F_ARCH}
    
  elif test x"${FC_ACT}" = x"mpixlf90_r"; then
    CFLAGS="-qsmp=omp -qthreaded -O3 -qessl -qstrict -qarch=qp -qtune=qp -qsimd=auto -qassert=contig -qdirective"
  
  elif test x"${FC_ACT}" = x"sxmpif90"; then
    CFLAGS="-Popenmp -Cvsafe"
  
  fi

fi


# CFLAGS0 (special flags)

if   test x"${FC_ACT}" = x"ifort"; then
  CFLAGS0=${I_BASE}" -O0 -Wall"
  
elif test x"${FC_ACT}" = x"gfortran"; then
  CFLAGS0=${G_BASE}" -O0 -Wall"
  
elif test x"${FC_ACT}" = x"mpifrtpx"; then
  CFLAGS0=${F_BASE}
  
elif test x"${FC_ACT}" = x"mpixlf90_r"; then
  CFLAGS0="-qsmp=omp -qthreaded -O0 -qessl -qstrict"

elif test x"${FC_ACT}" = x"sxmpif90"; then
  CFLAGS0="-Popenmp"

fi

AC_SUBST(CFLAGS0)


# LDFLAGS
if test -z "${LDFLAGS}"; then
  if   test x"${FC_ACT}" = x"ifort"; then
    LDFLAGS=${I_BASE}" -W1"
  
  elif test x"${FC_ACT}" = x"gfortran"; then
    LDFLAGS=${G_BASE}
  
  elif test x"${FC_ACT}" = x"mpifrtpx"; then
    LDFLAGS=${F_BASE}

  elif test x"${FC_ACT}" = x"mpixlf90_r"; then
    LDFLAGS="-qsmp=omp -qthreaded"
  
  elif test x"${FC_ACT}" = x"sxmpif90"; then
    LDFLAGS="-Popenmp -lmpi"
  
  fi
fi


# SHARED_LDFLAGS
if   test x"${FC_ACT}" = x"ifort"; then
  if test x"${OS_MACOS}" = x"yes"; then
    SHARED_LDFLAGS=" -shared"
  else
    SHARED_LDFLAGS=" -shared -rdynamic -Wl,-soname,libEigenExa.$SUFFIX_DYLIB -W1"
  fi
  
elif test x"${FC_ACT}" = x"gfortran"; then
  SHARED_LDFLAGS=" -shared -rdynamic -Wl,-soname,libEigenExa.$SUFFIX_DYLIB"
  
elif test x"${FC_ACT}" = x"mpifrtpx"; then
  SHARED_LDFLAGS=" -shared -rdynamic -Wl,-soname,libEigenExa.$SUFFIX_DYLIB"

elif test x"${FC_ACT}" = x"mpixlf90_r"; then
  SHARED_LDFLAGS=" -qmkshrobj "
  
elif test x"${FC_ACT}" = x"sxmpif90"; then
  SHARED_LDFLAGS=""
  
fi

AC_SUBST(SHARED_LDFLAGS)


# Test flag

if   test x"${FC_ACT}" = x"ifort"; then
  TEST_FLAG=""
  
elif test x"${FC_ACT}" = x"gfortran"; then
  TEST_FLAG=""
  
elif test x"${FC_ACT}" = x"mpifrtpx"; then
  TEST_FLAG=""
  
elif test x"${FC_ACT}" = x"mpixlf90_r"; then
  TEST_FLAG=""

elif test x"${FC_ACT}" = x"sxmpif90"; then
  TEST_FLAG="-f2003"

fi

AC_SUBST(TEST_FLAG)


# Macro "-D" prefix

if   test x"${FC_ACT}" = x"ifort"; then
  MACRO_D_PREFIX="-D"

elif test x"${FC_ACT}" = x"gfortran"; then
  MACRO_D_PREFIX="-D"

elif test x"${FC_ACT}" = x"mpifrtpx"; then
  MACRO_D_PREFIX="-D"

elif test x"${FC_ACT}" = x"mpixlf90_r"; then
  MACRO_D_PREFIX="-WF,-D"

elif test x"${FC_ACT}" = x"sxmpif90"; then
  MACRO_D_PREFIX="-D"

fi

AC_SUBST(MACRO_D_PREFIX)

# Macro Fortran C-pre-processor

if   test x"${FC_ACT}" = x"ifort"; then
  MACRO_CPP="-cpp"
  MACRO_FPP="-fpp"

elif test x"${FC_ACT}" = x"gfortran"; then
  MACRO_CPP="-cpp"
  MACRO_FPP="-cpp"

elif test x"${FC_ACT}" = x"mpifrtpx"; then
  MACRO_CPP=""
  MACRO_FPP="-Cpp"

elif test x"${FC_ACT}" = x"mpixlf90_r"; then
  MACRO_CPP=""
  MACRO_FPP=""

elif test x"${FC_ACT}" = x"sxmpif90"; then
  MACRO_CPP=""
  MACRO_FPP=""

fi

AC_SUBST(MACRO_FPP)


# Macro "-fPIC"

if   test x"${FC_ACT}" = x"ifort"; then
  if test x"${OS_MACOS}" = x"yes"; then
    MACRO_PIC="-fminshared"
  else
    MACRO_PIC="-fPIC"
  fi

elif test x"${FC_ACT}" = x"gfortran"; then
  MACRO_PIC="-fPIC"

elif test x"${FC_ACT}" = x"mpifrtpx"; then
  MACRO_PIC="-KPIC"

elif test x"${FC_ACT}" = x"mpixlf90_r"; then
  MACRO_PIC=""

elif test x"${FC_ACT}" = x"sxmpif90"; then
  MACRO_PIC=""

fi

AC_SUBST(MACRO_PIC)


# Macro "-module "

if   test x"${FC_ACT}" = x"ifort"; then
  MACRO_MODULE_PATH="-module "

elif test x"${FC_ACT}" = x"gfortran"; then
  MACRO_MODULE_PATH="-J"

elif test x"${FC_ACT}" = x"mpifrtpx"; then
  MACRO_MODULE_PATH="-M"

elif test x"${FC_ACT}" = x"mpixlf90_r"; then
  MACRO_MODULE_PATH=""

elif test x"${FC_ACT}" = x"sxmpif90"; then
  MACRO_MODULE_PATH=""

fi

AC_SUBST(MACRO_MODULE_PATH)


## setting for sxmpif90 (NEC SX-ACE)

if test x"${FC_ACT}" = x"sxmpif90"; then
  RANLIB=touch
  AR=sxar
else
  RANLIB=ranlib
  AR=ar
fi

AC_SUBST(AR)

AM_CONDITIONAL(SKIP_SHARED, test x"${FC_ACT}" = x"sxmpif90")
AM_CONDITIONAL(SKIP_SHARED, test ${USE_SHARED_LIB} -eq 0)
AM_CONDITIONAL(ORENAME_SHARED, test x"${FC_ACT}" = x"mpifrtpx")
AM_CONDITIONAL(SKIP_DEBUG, test x"${DEBUG}" = x0)
AM_CONDITIONAL(OS_MACOS, test x"${OS_MACOS}" = x"yes")

## eigen_s version
#

AC_ARG_ENABLE([eigen_s],,
    [
      if test x"${enableval}" = x"yes"; then
	    EIGEN_S_ARG=0
      elif test x"${enableval}" = x"0" -o \
                x"${enableval}" = x"1"; then
        EIGEN_S_ARG=${enableval}
      else
        AC_MSG_ERROR([Not supported Eigen_S version ${enableval}.])
      fi
    ],
	[
      EIGEN_S_ARG=0
    ])
AC_SUBST(EIGEN_S_ARG)

##

AC_CONFIG_FILES([Makefile
                src/Makefile
                C/Makefile
		benchmark_h/Makefile
		benchmark/Makefile])

AC_OUTPUT
