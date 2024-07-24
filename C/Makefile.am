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

noinst_PROGRAMS = c_test

c_test_CSOURCES = \
	EigenExa.c
c_test_FSOURCES = \
	eigen_exa_interfaces.F90

c_test_OBJS = $(c_test_FSOURCES:.F90=.o) \
		$(c_test_CSOURCES:.c=.o)

DEFS = 

c_test_LDFLAGS = -Wl,-rpath,$(PWD)/../src -L$(PWD)/../src

c_test_LDADD = $(c_test_OBJS) \
		-lEigenExa \
		$(OPT_LD_LAPACK)

AM_FFLAGS = -g $(MACRO_FPP) \
		$(TEST_FLAG) -I./ \
		$(MACRO_MODULE_PATH)./ \
		$(TEST_FLAG) -I../src \
		$(MACRO_MODULE_PATH)../src

CLEANFILES = *.o

c_test.o: c_test.c
	$(CC) -c -o $@ $<
EigenExa.o: EigenExa.c EigenExa.h EigenExa.fh
	$(CC) -c -o $@ $< 
eigen_exa_interfaces.o: eigen_exa_interfaces.F90
	$(F77) -c -o $@ $< $(AM_FFLAGS)

