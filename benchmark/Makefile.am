#--------*---------*---------*---------*---------*---------*---------*-*
#
#  File     : Makefile.am
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

AUTOMAKE_OPTIONS = foreign nostdinc

if SKIP_SHARED
noinst_PROGRAMS = eigenexa_benchmark eigenexa_gev_benchmark
else
noinst_PROGRAMS = eigenexa_benchmark eigenexa_gev_benchmark \
		eigenexa_benchmark_shlib eigenexa_gev_benchmark_shlib
endif

eigenexa_benchmark_SOURCES = \
	mat_set.f \
	ev_test.f \
	w_test.f \
	main2.f

eigenexa_gev_benchmark_SOURCES = \
	mat_set.f \
	KMATH_EIGEN_GEV_check.f \
	KMATH_EIGEN_GEV_main.f

eigenexa_benchmark_EXTRAS = \
	eigen_libs_mod.mod \
	eigen_libs0_mod.mod \
	fs_libs_mod.mod \
	eigen_devel_mod.mod \
	eigen_blacs_mod.mod

if !SKIP_SHARED
eigenexa_benchmark_shlib_SOURCES = $(eigenexa_benchmark_SOURCES)
eigenexa_gev_benchmark_shlib_SOURCES = $(eigenexa_gev_benchmark_SOURCES)
endif
eigenexa_benchmark_OBJS = $(eigenexa_benchmark_SOURCES:.f=.o)
eigenexa_gev_benchmark_OBJS = $(eigenexa_gev_benchmark_SOURCES:.f=.o)

DEFS = 

if !SKIP_SHARED
eigenexa_benchmark_shlib_DEPENDENCIES = libEigenExa.so \
	 $(am__DEPENDENCIES_1)
eigenexa_gev_benchmark_shlib_DEPENDENCIES = libEigenExa.so \
	 $(am__DEPENDENCIES_1)
endif

eigenexa_benchmark_LDFLAGS = -Wl,-rpath,. -L./
eigenexa_gev_benchmark_LDFLAGS = -Wl,-rpath,. -L./
if !SKIP_SHARED
eigenexa_benchmark_shlib_LDFLAGS = $(eigenexa_benchmark_LDFLAGS)
eigenexa_gev_benchmark_shlib_LDFLAGS = $(eigenexa_benchmark_LDFLAGS)
endif

eigenexa_benchmark_LDADD = $(top_builddir)/src/libEigenExa.a $(OPT_LD_LAPACK)
eigenexa_gev_benchmark_LDADD = $(top_builddir)/src/libEigenExa.a $(OPT_LD_LAPACK)
if !SKIP_SHARED
eigenexa_benchmark_shlib_LDADD = -lEigenExa $(OPT_LD_LAPACK)
eigenexa_gev_benchmark_shlib_LDADD = -lEigenExa $(OPT_LD_LAPACK)
endif

#AM_FFLAGS = -g $(MACRO_FPP) \
#		$(MACRO_D_PREFIX)TIMER_PRINT=$(OPT_TIMER_PRINT) \
#		$(TEST_FLAG) -I$(top_srcdir)/src \
#		$(MACRO_MODULE_PATH)$(top_srcdir)/src
AM_FFLAGS = -g $(MACRO_FPP) \
		$(MACRO_D_PREFIX)TIMER_PRINT=$(OPT_TIMER_PRINT) \
		$(TEST_FLAG) -I./ \
		$(MACRO_MODULE_PATH)./

CLEANFILES = *.mod *.o libEigenExa.a libEigenExa.so

main2.o: main2.f $(eigenexa_benchmark_EXTRAS)
	$(F77) $(AM_FFLAGS) $(FFLAGS0) -c -o $@ $<

KMATH_EIGEN_GEV_main.o : KMATH_EIGEN_GEV_main.f
	$(F77) $(AM_FFLAGS) $(FFLAGS0) -c -o $@ $<

$(eigenexa_benchmark_OBJS) : \
	libEigenExa.a $(eigenexa_benchmark_EXTRAS)

$(eigenexa_gev_benchmark_OBJS) : \
	libEigenExa.a $(eigenexa_benchmark_EXTRAS)

if !SKIP_SHARED
$(eigenexa_benchmark_shlib_OBJS) : \
	libEigenExa.so $(eigenexa_benchmark_EXTRAS)

$(eigenexa_gev_benchmark_shlib_OBJS) : \
	libEigenExa.so $(eigenexa_benchmark_EXTRAS)
endif

$(eigenexa_benchmark_EXTRAS):
	ln -s $(top_srcdir)/src/$@ .

libEigenExa.a: $(top_builddir)/src/libEigenExa.a
	ln -s $(top_srcdir)/src/$@ .
if !SKIP_SHARED
libEigenExa.so: $(top_builddir)/src/libEigenExa.so
	ln -s $(top_srcdir)/src/$@ .
endif

