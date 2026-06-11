#--------*---------*---------*---------*---------*---------*---------*-*
#
#  File     : Makefile.am
#
# Copyright(C) 2012-2021 RIKEN.
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
noinst_PROGRAMS = eigenexah_benchmark eigenexa_benchmark
else
noinst_PROGRAMS = eigenexah_benchmark \
		eigenexah_benchmark_shlib \
		eigenexa_benchmark \
		eigenexa_benchmark_shlib
endif

eigenexah_benchmark_SOURCES = \
	mat_set_h.f \
	bench_eigen_h.f
eigenexa_benchmark_SOURCES = \
	mat_set.f \
	bench_eigen.f


eigenexah_benchmark_EXTRAS = \
	eigen_libs_mod.mod \
	eigen_libs0_mod.mod \
	fs_libs_mod.mod \
	eigen_devel_mod.mod \
	eigen_blacs_mod.mod
eigenexa_benchmark_EXTRAS = $(eigenexah_benchmark_EXTRAS)


if !SKIP_SHARED
eigenexah_benchmark_shlib_SOURCES = $(eigenexah_benchmark_SOURCES)
eigenexa_benchmark_shlib_SOURCES = $(eigenexa_benchmark_SOURCES)
endif
eigenexah_benchmark_OBJS = $(eigenexah_benchmark_SOURCES:.f=.o)
eigenexa_benchmark_OBJS = $(eigenexa_benchmark_SOURCES:.f=.o)

DEFS = 

if !SKIP_SHARED
eigenexah_benchmark_shlib_DEPENDENCIES = libEigenExa.so \
	 $(am__DEPENDENCIES_1)
eigenexa_benchmark_shlib_DEPENDENCIES = libEigenExa.so \
	 $(am__DEPENDENCIES_1)
endif

eigenexah_benchmark_LDFLAGS = -Wl,-rpath,. -L./
eigenexa_benchmark_LDFLAGS = -Wl,-rpath,. -L./
if !SKIP_SHARED
eigenexah_benchmark_shlib_LDFLAGS = $(eigenexah_benchmark_LDFLAGS)
eigenexa_benchmark_shlib_LDFLAGS = $(eigenexa_benchmark_LDFLAGS)
endif

eigenexah_benchmark_LDADD = $(top_builddir)/src/libEigenExa.a $(OPT_LD_LAPACK)
eigenexa_benchmark_LDADD = $(top_builddir)/src/libEigenExa.a $(OPT_LD_LAPACK)
if !SKIP_SHARED
eigenexah_benchmark_shlib_LDADD = -lEigenExa $(OPT_LD_LAPACK)
eigenexa_benchmark_shlib_LDADD = -lEigenExa $(OPT_LD_LAPACK)
endif

AM_FFLAGS = -g $(MACRO_FPP) \
		$(MACRO_D_PREFIX)TIMER_PRINT=$(OPT_TIMER_PRINT) \
		$(TEST_FLAG) -I./ \
		$(MACRO_MODULE_PATH)./

CLEANFILES = *.mod *.o libEigenExa.a libEigenExa.so

bench_eigen_h.o: bench_eigen_h.f $(eigenexah_benchmark_EXTRAS)
	$(F77) $(AM_FFLAGS) $(FFLAGS0) -c -o $@ $<
bench_eigen.o: bench_eigen.f $(eigenexa_benchmark_EXTRAS)
	$(F77) $(AM_FFLAGS) $(FFLAGS0) -c -o $@ $<

$(eigenexah_benchmark_OBJS) : \
	libEigenExa.a $(eigenexah_benchmark_EXTRAS)
$(eigenexa_benchmark_OBJS) : \
	libEigenExa.a $(eigenexa_benchmark_EXTRAS)

if !SKIP_SHARED
$(eigenexah_benchmark_shlib_OBJS) : \
	libEigenExa.so $(eigenexah_benchmark_EXTRAS)
$(eigenexa_benchmark_shlib_OBJS) : \
	libEigenExa.so $(eigenexa_benchmark_EXTRAS)
endif

$(eigenexah_benchmark_EXTRAS): 
	ln -s $(top_srcdir)/src/$@ .

libEigenExa.a: $(top_builddir)/src/libEigenExa.a
	ln -s $(top_srcdir)/src/$@ .
if !SKIP_SHARED
libEigenExa.so: $(top_builddir)/src/libEigenExa.so
	ln -s $(top_srcdir)/src/$@ .
endif

