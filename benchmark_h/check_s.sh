#!/bin/bash
#PJM -L "node=10"               # 4ノード
#PJM -L "rscgrp=small"         # リソースグループの指定
#PJM -L "elapse=01:00:00"      # ジョブの経過時間制限値
#PJM -g hp230279           # 課題のグループ指定
#PJM -x PJM_LLIO_GFSCACHE=/vol0004 # ジョブで使用するデータ領域のvolume
#PJM --mpi "max-proc-per-node=4" # 1ノードあたりに生成するMPIプロセス数の上限値


INTEL_MPI="yes"
INTEL_MPI1=`ldd eigenexah_benchmark | awk '/libmpifort/{ print $1}'`
INTEL_MPI2=`which mpiexec.hydra | grep " no "`
if [ x$INTEL_MPI1 = x ]; then
	INTEL_MPI="no"
fi
if [ x$INTEL_MPI2 != x ]; then
	INTEL_MPI="no"
fi

if [ $INTEL_MPI = "yes" ]; then
	export MKL_DYNAMIC=FALSE
	export OMP_DYNAMIC=FALSE
	export I_MPI_CBWR=2
	export FI_SOCKETS_IFACE=enp4s0f0
	export FI_PROVIDER=sockets
	export I_MPI_FABRICS=shm
fi

\rm LOG-*-*

for P in 10 20 30 40 50 60 70 80; do
for T in 12; do
echo "P="$P" T="$T
export OMP_NUM_THREADS=$T
if [ $INTEL_MPI = "yes" ]; then
mpiexec.hydra \
        -np $P -genv OMP_NUM_THREADS $T \
        ./eigenexa_benchmark < /dev/null |& tee LOG-$P-$T
else

mpirun \
	-np $P -x OMP_NUM_THREADS=$T \
	./eigenexa_benchmark < /dev/null |& tee LOG-$P-$T
#mpirun \
#	-np $P -x OMP_NUM_THREADS=$T \
#	~/default/EigenExa-2.12/benchmark_h/eigenexa_benchmark < /dev/null |& tee LOG-$P-$T

fi
done
done

