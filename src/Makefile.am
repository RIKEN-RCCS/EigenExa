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

instdir = ${prefix}/lib

inst_LIBRARIES = libEigenExa.a
if !SKIP_SHARED
inst_PROGRAMS  = $(inst_LIBRARIES:.a=.so)
endif

libEigenExa_a_SOURCES = \
	MPI_Allreduce_group.F90 \
	CSTAB.F \
	eigen_house.F \
	eigen_house_h.F \
	eigen_dc.F \
	eigen_trbak.F \
	eigen_hrbak.F \
	eigen_blacs.F \
	eigen_devel.F \
	comm.F \
	comm_h.F \
	eigen_libs0.F \
	eigen_libs.F \
	FS_const.F90 \
	FS_libs.F90 \
	FS_prof.F90 \
	dlacpy.F \
	lapack_eigen.F \
	my_pdlaedz.F \
	my_pdlaed3.F \
	my_pdlaed2.F \
	my_pdlaed1.F \
	my_pdlaed0.F \
	my_pdlasrt.F \
	my_pdsxedc.F \
	mx_pdlaedz.F \
	mx_pdlaed3.F \
	mx_pdlaed2.F \
	mx_pdlaed1.F \
	mx_pdlaed0.F \
	mx_pdstedc.F \
	dc_redist1.F \
	dc_redist1_r.F \
	dc_redist2.F \
	dc_redist2_r.F \
	dc2.F \
	dcx.F \
	bisect.F \
	bisect2.F \
	eigen_t1.F \
	eigen_trd_t2.F \
	eigen_trd_t4.F \
	eigen_trd_t5.F \
	eigen_trd_t5x.F \
	eigen_trd_t6_3.F \
	eigen_trd_t7.F \
	eigen_trd_t8.F \
	eigen_trd.F \
	eigen_prd_t2.F \
	eigen_prd_t4x.F \
	eigen_prd_t5.F \
	eigen_prd_t6_3.F \
	eigen_prd_t7.F \
	eigen_prd_t8.F \
	eigen_prd.F \
	trbakwy4_body.F \
	trbakwy4.F \
	eigen_scaling.F \
	eigen_identity.F \
	eigen_sx.F \
	eigen_s.F  \
	FS_dividing.F90 \
	FS_EDC.F90 \
	FS_PDLAED0.F90 \
	FS_PDLAED1.F90 \
	FS_PDLAED2.F90 \
	FS_PDLAED3.F90 \
	FS_PDLAEDZ.F90 \
	FS_MERGE_D.F90 \
	FS_REDUCE_ZD.F90 \
	FS_PDLASRT.F90 \
	FS2eigen_PDLASRT.F90 \
	dc2_FS.F \
	eigen_FS.F \
	FS_memory_internal.F \
	eigen_hrd_t1.F \
	eigen_hrd_t2.F \
	eigen_hrd_t4.F \
	eigen_hrd_t5.F \
	eigen_hrd_t5x.F \
	eigen_hrd_t6_3.F \
	eigen_hrd_t7.F \
	eigen_hrd_t8.F \
	eigen_hrd.F \
	eigen_h.F  \
	hrbakwy4_body.F \
	hrbakwy4.F \
	eigen_scaling_h.F \
	KMATH_EIGEN_GEV.F

if !SKIP_SHARED
libEigenExa_a_OBJECTS_shared = $(libEigenExa_a_OBJECTS:.o=.lo)
$(libEigenExa_a_OBJECTS_shared): $(inst_LIBRARIES)
endif

if !SKIP_SHARED
libEigenExa_so_SOURCES =
libEigenExa.so$(EXEEXT): $(libEigenExa_a_OBJECTS_shared)
endif
if !SKIP_SHARED
	$(F77) $(SHARED_LDFLAGS) -o $@ $(libEigenExa_a_OBJECTS_shared)
	chmod -x $@
endif

CLEANFILES = *.mod *.lo

DEFS = 

EXT_FFLAGS = $(MACRO_D_PREFIX)AT_BCAST_OVERLAP=1 $(MACRO_D_PREFIX)AT_REDUCE_OVERLAP=0 $(MACRO_D_PREFIX)OVERLAP_DECISION_TYPE=3

AM_FFLAGS = $(MACRO_D_PREFIX)TIMER_PRINT=$(OPT_TIMER_PRINT) $(MACRO_MODULE_PATH)./ $(MACRO_D_PREFIX)BOOST_BY_CACHE_ALIGNMENT=$(OPT_CACHE_ALIGN)

AM_FCFLAGS = $(MACRO_D_PREFIX)TIMER_PRINT=$(OPT_TIMER_PRINT) $(MACRO_MODULE_PATH)./ $(MACRO_D_PREFIX)BOOST_BY_CACHE_ALIGNMENT=$(OPT_CACHE_ALIGN)

FS_FFLAGS = $(MACRO_D_PREFIX)CODE_AKASHI=${EIGEN_S_ARG}

#include_HEADERS = eigen_libs_mod.mod \
#	eigen_libs0_mod.mod \
#	eigen_blacs_mod.mod \
#	comm_mod.mod \
#	dc2_fs_mod.mod \
#	fs2eigen.mod \
#	fs_const_mod.mod \
#	fs_dividing_mod.mod \
#	fs_edc_mod.mod \
#	fs_libs_mod.mod \
#	fs_prof_mod.mod
#
include_HEADERS = \
	eigen_devel_mod.mod \
	eigen_libs_mod.mod \
	eigen_libs0_mod.mod \
	fs_libs_mod.mod \
	eigen_blacs_mod.mod \
	comm_mod.mod


.F.o:
	$(AM_V_PPF77)$(PPF77COMPILE) -c -o $@ $<
.F90.o:
	$(AM_V_PPF77)$(PPF77COMPILE) -c -o $@ $<


if !SKIP_SHARED
.F.lo:
if !ORENAME_SHARED
	$(AM_V_PPF77)$(PPF77COMPILE) $(MACRO_PIC) -c -o $@ $<
else
	$(AM_V_PPF77)$(PPF77COMPILE) $(MACRO_PIC) -c -o tmp_$(@:.lo=.o) $<
	mv tmp_$(@:.lo=.o) $@
endif
endif

if !SKIP_SHARED
.F90.lo:
if !ORENAME_SHARED
	$(AM_V_PPF77)$(PPF77COMPILE) $(MACRO_PIC) -c -o $@ $<
else
	$(AM_V_PPF77)$(PPF77COMPILE) $(MACRO_PIC) -c -o tmp_$(@:.lo=.o) $<
	mv tmp_$(@:.lo=.o) $@
endif
endif



trbakwy4.o : trbakwy4.F
	$(AM_V_PPFC)$(PPFCCOMPILE) $(EXT_FFLAGS) -c -o $@ $<
trbakwy4_body.o : trbakwy4_body.F
	$(AM_V_PPFC)$(PPFCCOMPILE) $(EXT_FFLAGS) -c -o $@ $<
FS_libs.o : FS_libs.F90
	$(F77) $(AM_FFLAGS) $(FFLAGS0) -c -o $@ $<
eigen_libs.o : eigen_libs.F
	$(F77) $(AM_FFLAGS) $(FFLAGS0) -c -o $@ $<
eigen_devel.o : eigen_devel.F
	$(F77) $(AM_FFLAGS) $(FFLAGS0) -c -o $@ $<
comm.o : comm.F
	$(F77) $(AM_FFLAGS) $(FFLAGS0) -c -o $@ $<
comm_h.o : comm_h.F
	$(F77) $(AM_FFLAGS) $(FFLAGS0) -c -o $@ $<
FS_dividing.o : FS_dividing.F90
	$(F77) $(AM_FFLAGS) $(FFLAGS0) -c -o $@ $<
FS_const.o : FS_const.F90
	$(F77) $(AM_FFLAGS) $(FFLAGS0) -c -o $@ $<
FS_prof.o : FS_prof.F90
	$(F77) $(AM_FFLAGS) $(FFLAGS0) -c -o $@ $<
miscC.o : miscC.c
	$(CC) $(AM_CFLAGS) $(CFLAGS0) -c -o $@ $<


if !SKIP_SHARED
trbakwy4.lo : trbakwy4.F
	$(AM_V_PPFC)$(PPFCCOMPILE) $(EXT_FFLAGS) $(MACRO_PIC) -c -o tmp_$(@:.lo=.o) $<
	mv tmp_$(@:.lo=.o) $@
trbakwy4_body.lo : trbakwy4_body.F
	$(AM_V_PPFC)$(PPFCCOMPILE) $(EXT_FFLAGS) $(MACRO_PIC) -c -o tmp_$(@:.lo=.o) $<
	mv tmp_$(@:.lo=.o) $@
FS_libs.lo : FS_libs.F90
	$(F77) $(AM_FFLAGS) $(FFLAGS0) $(MACRO_PIC) -c -o tmp_$(@:.lo=.o) $<
	mv tmp_$(@:.lo=.o) $@
eigen_libs.lo : eigen_libs.F
	$(F77) $(AM_FFLAGS) $(FFLAGS0) $(MACRO_PIC) -c -o tmp_$(@:.lo=.o) $<
	mv tmp_$(@:.lo=.o) $@
eigen_devel.lo : eigen_devel.F
	$(F77) $(AM_FFLAGS) $(FFLAGS0) $(MACRO_PIC) -c -o tmp_$(@:.lo=.o) $<
	mv tmp_$(@:.lo=.o) $@
comm.lo : comm.F
	$(F77) $(AM_FFLAGS) $(FFLAGS0) $(MACRO_PIC) -c -o tmp_$(@:.lo=.o) $<
	mv tmp_$(@:.lo=.o) $@
comm_h.lo : comm_h.F
	$(F77) $(AM_FFLAGS) $(FFLAGS0) $(MACRO_PIC) -c -o tmp_$(@:.lo=.o) $<
	mv tmp_$(@:.lo=.o) $@
FS_dividing.lo : FS_dividing.F90
	$(F77) $(AM_FFLAGS) $(FFLAGS0) $(MACRO_PIC) -c -o tmp_$(@:.lo=.o) $<
	mv tmp_$(@:.lo=.o) $@
FS_const.lo : FS_const.F90
	$(F77) $(AM_FFLAGS) $(FFLAGS0) $(MACRO_PIC) -c -o tmp_$(@:.lo=.o) $<
	mv tmp_$(@:.lo=.o) $@
FS_prof.lo : FS_prof.F90
	$(F77) $(AM_FFLAGS) $(FFLAGS0) $(MACRO_PIC) -c -o tmp_$(@:.lo=.o) $<
	mv tmp_$(@:.lo=.o) $@
miscC.lo : miscC.c
	$(CC) $(AM_CFLAGS) $(CFLAGS0) -fPIC -c -o $@ $<
endif




if !SKIP_SHARED
$(include_HEADERS) : $(inst_LIBRARIES)
else
$(include_HEADERS) : $(inst_LIBRARIES) $(inst_PROGRAMS)
endif

eigen_libs.o	: eigen_devel.o eigen_libs0.o FS_libs.o
FS_libs.o	: FS_const.o eigen_libs0.o
eigen_libs0.o	: comm.o comm_h.o eigen_blacs.o eigen_devel.o CSTAB.o
eigen_blacs.o	: eigen_devel.o
KMATH_EIGEN_GEV.o	: eigen_libs.o
comm.o          : eigen_devel.o
comm_h.o        : comm.o eigen_devel.o
FS_prof.o       : FS_libs.o
FS_dividing.o   : FS_prof.o

libEigenExa_core_objs = \
	eigen_devel.o \
	eigen_libs.o \
	eigen_libs0.o \
	FS_libs.o \
	eigen_blacs.o \
	comm.o \
	comm_h.o \
	FS_prof.o \
	FS_dividing.o \
	FS_const.o \
	MPI_Allreduce_group.o \
	KMATH_EIGEN_GEV.o \
	CSTAB.o
libEigenExa_non_core_objs = \
	$(foreach obj, $(libEigenExa_a_OBJECTS), \
		$(word $(words $(libEigenExa_core_objs) ), \
			$(foreach core, $(libEigenExa_core_objs), \
				$(filter-out $(core), $(obj) ) \
			) \
		) \
	)
$(libEigenExa_non_core_objs) : $(libEigenExa_core_objs)


if !SKIP_SHARED
eigen_libs.lo	: eigen_devel.lo eigen_libs0.lo FS_libs.lo
FS_libs.lo	: FS_const.lo eigen_libs0.lo
eigen_libs0.lo	: comm.lo comm_h.lo eigen_blacs.lo eigen_devel.lo CSTAB.lo
eigen_blacs.lo	: eigen_devel.lo
KMATH_EIGEN_GEV.lo	: eigen_libs.lo
comm.lo		: eigen_devel.lo
comm_h.lo	: comm.lo eigen_devel.lo
FS_prof.lo	: FS_libs.lo
FS_dividing.lo	: FS_prof.lo
$(libEigenExa_non_core_objs:.o=.lo) : $(libEigenExa_core_objs:.o=.lo)
endif



eigen_sx.o eigen_s.o eigen_FS.o eigen_h.o: \
	eigen_scaling.o eigen_scaling_h.o \
	eigen_trd.o eigen_prd.o eigen_hrd.o \
	bisect.o \
	dc2.o dcx.o dc2_FS.o \
	trbakwy4.o hrbakwy4.o 

TRD_OBJS = \
	eigen_trd_t2.o \
	eigen_trd_t4.o \
	eigen_trd_t5.o eigen_trd_t5x.o \
	eigen_trd_t6_3.o \
	eigen_trd_t7.o \
	eigen_trd_t8.o
PRD_OBJS = \
	eigen_prd_t2.o \
	eigen_prd_t4x.o \
	eigen_prd_t5.o \
	eigen_prd_t6_3.o \
	eigen_prd_t7.o \
	eigen_prd_t8.o
HRD_OBJS = \
	eigen_hrd_t2.o \
	eigen_hrd_t4.o \
	eigen_hrd_t5.o eigen_hrd_t5x.o \
	eigen_hrd_t6_3.o \
	eigen_hrd_t7.o \
	eigen_hrd_t8.o

eigen_trd.o : $(TRD_OBJS) eigen_t1.o
eigen_prd.o : $(PRD_OBJS) eigen_t1.o
eigen_hrd.o : $(HRD_OBJS) eigen_hrd_t1.o
$(TRD_OBJS) $(PRD_OBJS) eigen_t1.o : eigen_house.o
$(HRD_OBJS) eigen_hrd_t1.o : eigen_house_h.o

MY_DC_OBJS = \
	my_pdlaed0.o \
	my_pdlaed1.o \
	my_pdlaed2.o \
	my_pdlaed3.o \
	my_pdlaedz.o \
	my_pdsxedc.o \
	my_pdlasrt.o
MX_DC_OBJS = \
	mx_pdlaed0.o \
	mx_pdlaed2.o \
	mx_pdlaedz.o \
	mx_pdlaed1.o \
	mx_pdlaed3.o \
	mx_pdstedc.o
FS_DC_OBJS = \
	FS_PDLAED0.o \
	FS_PDLAED1.o \
	FS_PDLAED2.o \
	FS_PDLAED3.o \
	FS_PDLASRT.o \
	FS_MERGE_D.o \
	FS_REDUCE_ZD.o \
	FS2eigen_PDLASRT.o \
	FS_EDC.o
DC_COMMON_OBJS = \
	eigen_dc.o \
	lapack_eigen.o \
	dlacpy.o

dcx.o: $(MY_DC_OBJS) $(DC_COMMON_OBJS)
dc2.o: $(MX_DC_OBJS) $(DC_COMMON_OBJS)
dc2_FS.o: $(FS_DC_OBJS)
$(MY_DC_OBJS) $(MX_DC_OBJS) $(FS_DC_OBJS): $(DC_COMMON_OBJS)

my_pdsxedc.o: my_pdlasrt.o my_pdlaed0.o
my_pdlaed0.o: my_pdlaed1.o my_pdlaedz.o
my_pdlaed1.o: my_pdlaed2.o my_pdlaed3.o
mx_pdstedc.o: my_pdlasrt.o mx_pdlaed0.o 
mx_pdlaed0.o: mx_pdlaed1.o
mx_pdlaed1.o: mx_pdlaed2.o mx_pdlaed3.o mx_pdlaedz.o

trbakwy4.o: trbakwy4_body.o eigen_trbak.o
trbakwy4_body.o: eigen_trbak.o
hrbakwy4.o: hrbakwy4_body.o eigen_hrbak.o
hrbakwy4_body.o: eigen_hrbak.o

if !SKIP_SHARED
eigen_sx.lo eigen_s.lo eigen_FS.lo eigen_h.lo: \
	eigen_scaling.lo eigen_scaling_h.lo \
	eigen_trd.lo eigen_prd.lo eigen_hrd.lo \
	bisect.lo \
	dc2.lo dcx.lo dc2_FS.lo \
	trbakwy4.lo hrbakwy4.lo 
eigen_trd.lo: $(TRD_OBJS:.o=.lo) eigen_t1.lo
eigen_prd.lo: $(PRD_OBJS:.o=.lo) eigen_t1.lo
eigen_hrd.lo: $(HRD_OBJS:.o=.lo)
$(TRD_OBJS:.o=.lo) $(PRD_OBJS:.o=.lo) eigen_t1.lo: eigen_house.lo
$(HRD_OBJS:.o=.lo) eigen_hrd_t1.lo: eigen_house_h.lo
dcx.lo: $(MY_DC_OBJS:.o=.lo) $(DC_COMMON_OBJS:.o=.lo)
dc2.lo: $(MX_DC_OBJS:.o=.lo) $(DC_COMMON_OBJS:.o=.lo)
dc2_FS.lo: $(FS_DC_OBJS:.o=.lo)
$(MY_DC_OBJS:.o=.lo) $(MX_DC_OBJS:.o=.lo) $(FS_DC_OBJS:.o=.lo): $(DC_COMMON_OBJS:.o=.lo)
my_pdsxedc.lo: my_pdlasrt.lo my_pdlaed0.lo
my_pdlaed0.lo: my_pdlaed1.lo my_pdlaedz.lo
my_pdlaed1.lo: my_pdlaed2.lo my_pdlaed3.lo
mx_pdstedc.lo: my_pdlasrt.lo mx_pdlaed0.lo 
mx_pdlaed0.lo: mx_pdlaed1.lo
mx_pdlaed1.lo: mx_pdlaed2.lo mx_pdlaed3.lo mx_pdlaedz.lo
trbakwy4.lo: trbakwy4_body.lo eigen_trbak.lo
trbakwy4_body.lo: eigen_trbak.lo
hrbakwy4.lo: hrbakwy4_body.lo eigen_hrbak.lo
hrbakwy4_body.lo: eigen_hrbak.lo
endif

