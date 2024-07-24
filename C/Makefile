MPIFC = mpiifort -nofor-main
MPICC = mpiicx -Wall
#MPIFC = mpif90
#MPICC = mpigcc

a.out: test.o EigenExa.o
	$(MPIFC)  test.o EigenExa.o \
	-Wl,-rpath,$(PWD)/../src -L$(PWD)/../src -lEigenExa \
	-L/mnt/nasforlpnct/opt/intel/oneapi/mkl/2023.2.0/lib/intel64 \
	-lmkl_scalapack_lp64 \
	-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
	-lmkl_blacs_intelmpi_lp64 \
	-liomp5 -lpthread -lm -ldl

test.o: test.c
	$(MPICC) -c test.c

EigenExa.o: EigenExa.c EigenExa.h
	$(MPICC) -c EigenExa.c -Wno-implicit-function-declaration


clean:
	\rm a.out *.o
