#!/bin/bash

cat << EOS > IN-Random
!    N   nvec bx  by m t s e
  4096   4096 48 128 1 2 0 1
  4096   4096 48 128 1 2 1 1
  4096   4096 48 128 1 2 0 1
  4096   4096 48 128 1 2 1 1
  8192   8192 48 128 1 2 0 1
  8192   8192 48 128 1 2 1 1
  8192   8192 48 128 1 2 0 1
  8192   8192 48 128 1 2 1 1
 10000  10000 48 128 1 2 0 1
 10000  10000 48 128 1 2 1 1
 10000  10000 48 128 1 2 0 1
 10000  10000 48 128 1 2 1 1
 16384  16384 48 128 1 2 0 1
 16384  16384 48 128 1 2 1 1
 16384  16384 48 128 1 2 0 1
 16384  16384 48 128 1 2 1 1
 20000  20000 48 128 1 2 0 1
 20000  20000 48 128 1 2 1 1
 20000  20000 48 128 1 2 0 1
 20000  20000 48 128 1 2 1 1
 32768  32768 48 128 1 2 0 1
 32768  32768 48 128 1 2 1 1
 32768  32768 48 128 1 2 0 1
 32768  32768 48 128 1 2 1 1
 40000  40000 48 128 1 2 0 1
 40000  40000 48 128 1 2 1 1
 40000  40000 48 128 1 2 0 1
 40000  40000 48 128 1 2 1 1
 65536  65536 48 128 1 2 0 1
 65536  65536 48 128 1 2 1 1
 65536  65536 48 128 1 2 0 1
 65536  65536 48 128 1 2 1 1
131072 131072 48 128 1 2 0 1
131072 131072 48 128 1 2 1 1
131072 131072 48 128 1 2 0 1
131072 131072 48 128 1 2 1 1
-1 0 48 128 1 2 1 1
EOS
cat << EOS > IN-Frank
!    N   nvec bx  by m t s e
  4096   4096 48 128 1 0 0 1
  4096   4096 48 128 1 0 1 1
  4096   4096 48 128 1 0 0 1
  4096   4096 48 128 1 0 1 1
  8192   8192 48 128 1 0 0 1
  8192   8192 48 128 1 0 1 1
  8192   8192 48 128 1 0 0 1
  8192   8192 48 128 1 0 1 1
 10000  10000 48 128 1 0 0 1
 10000  10000 48 128 1 0 1 1
 10000  10000 48 128 1 0 0 1
 10000  10000 48 128 1 0 1 1
 16384  16384 48 128 1 0 0 1
 16384  16384 48 128 1 0 1 1
 16384  16384 48 128 1 0 0 1
 16384  16384 48 128 1 0 1 1
 20000  20000 48 128 1 0 0 1
 20000  20000 48 128 1 0 1 1
 20000  20000 48 128 1 0 0 1
 20000  20000 48 128 1 0 1 1
 32768  32768 48 128 1 0 0 1
 32768  32768 48 128 1 0 1 1
 32768  32768 48 128 1 0 0 1
 32768  32768 48 128 1 0 1 1
 40000  40000 48 128 1 0 0 1
 40000  40000 48 128 1 0 1 1
 40000  40000 48 128 1 0 0 1
 40000  40000 48 128 1 0 1 1
 65536  65536 48 128 1 0 0 1
 65536  65536 48 128 1 0 1 1
 65536  65536 48 128 1 0 0 1
 65536  65536 48 128 1 0 1 1
131072 131072 48 128 1 0 0 1
131072 131072 48 128 1 0 1 1
131072 131072 48 128 1 0 0 1
131072 131072 48 128 1 0 1 1
-1 0 48 128 1 2 1 1
EOS

for d0 in 1ppn 2ppn 4ppn; do
	if [ ! -d $d0 ]; then
		mkdir $d0
	fi
	cd $d0

	if [ -f eigenexah_benchmark ]; then
		\rm eigenexah_benchmark
	fi
	if [ -f eigenexa_benchmark ]; then
		\rm eigenexa_benchmark
	fi
	if [ -f IN-Random ]; then
		\rm IN-Random
	fi
	if [ -f IN-Frank ]; then
		\rm IN-Frank
	fi

	#ln -s ../eigenexah_benchmark .
	ln -s ../eigenexa_benchmark .
	ln -s ../IN-Random .
	ln -s ../IN-Frank .

	for d1 in 00001 00002 00004 00008 00016 00032 00064 00128 00256 00512 01024 02048 04096 08192; do
		if [ ! -d $d1 ]; then
			mkdir $d1
		fi
		cd $d1

		if [ -f eigenexah_benchmark ]; then
			\rm eigenexah_benchmark
		fi
		if [ -f eigenexa_benchmark ]; then
			\rm eigenexa_benchmark
		fi
		if [ -f IN-Random ]; then
			\rm IN-Random
		fi
		if [ -f IN-Frank ]; then
			\rm IN-Frank
		fi

		#ln -s ../eigenexah_benchmark .
		ln -s ../eigenexa_benchmark .
		ln -s ../IN-Random .
		ln -s ../IN-Frank .
		
		if [ $d1 -le 384 ]; then
			RSC=small
			ELP='24:00:00'
		else
			RSC=large
			ELP='24:00:00'
		fi

		if [ $d0 = '1ppn' ]; then
			PPN=1
			THR=48
		fi
		if [ $d0 = '2ppn' ]; then
			PPN=2
			THR=24
		fi
		if [ $d0 = '4ppn' ]; then
			PPN=4
			THR=12
		fi

		NOD=`echo $d1 | awk '{c=int(sqrt($0)); if(c*c!=$0){c=int(sqrt($0/2))}print $0/c"x"c}'	`
		PE=`expr $d1 "*" $PPN`

		FILE=run-$d1-$d0-large.sh

		cat << EOS > $FILE
#!/bin/sh
#PJM -L "node=$NOD"
#PJM -L "rscunit=rscunit_ft01"
#PJM -L "rscgrp=$RSC"
#PJM -L "elapse=$ELP"
#PJM -L "freq=2200,eco_state=2"
#PJM --mpi "max-proc-per-node=$PPN"
#PJM -s


export PLE_MPI_STD_EMPTYFILE=off
export OMP_NUM_THREADS=$THR

export NOD=$NOD
export PE=$PE
export PPN=$PPN

export OMPI_MCA_plm_ple_memory_allocation_policy=interleave_local
export OMPI_MCA_plm_ple_numanode_assign_policy=simplex

printenv

date
/usr/bin/llio_transfer ./IN-Frank
/usr/bin/llio_transfer ./IN-Random
/usr/bin/llio_transfer ./eigenexa_benchmark
date
mpiexec \
        -stdout-proc ./output.%j/%/1000r/stdout \
        -stderr-proc ./output.%j/%/1000r/stderr \
        ./eigenexa_benchmark -f IN-Frank
date
mpiexec \
        -stdout-proc ./output.%j/%/1000r/stdout \
        -stderr-proc ./output.%j/%/1000r/stderr \
        ./eigenexa_benchmark -f IN-Random
date
/usr/bin/llio_transfer --purge ./eigenexa_benchmark
date

EOS
		cd ..
	done
	cd ..
done

