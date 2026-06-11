#!/bin/bash
#PJM -L "node=400"               # 4ノード
#PJM -L "rscgrp=large"         # リソースグループの指定
#PJM -L "elapse=02:00:00"      # ジョブの経過時間制限値
#PJM -g hp230279           # 課題のグループ指定
#PJM -x PJM_LLIO_GFSCACHE=/vol0004 # ジョブで使用するデータ領域のvolume
#PJM --mpi "max-proc-per-node=4" # 1ノードあたりに生成するMPIプロセス数の上限値


if [ -f IN-check ]; then
  \rm IN-check
fi

#echo "!   N  nvec bx  by m t s e" > IN-check
#awk 'BEGIN{ for(N=3;N<=256;N++){ \
#	print N" "N" 48 128 1 0 0 1"; \
#	print N" "N" 48 128 1 0 1 1"; \
#	print N" "N" 48 128 1 2 0 1"; \
#	print N" "N" 48 128 1 2 1 1"; \
#   } exit;}END{}' >> IN-check
#60000 46967
for N in 60029; do
  echo $N" "$N" 48 128 1 0 1 1" >> IN-check #Frank 64
  echo $N" "$N" 48 128 1 2 1 1" >> IN-check #Randam 64
  echo $N" "$N" 48 128 1 2 2 1" >> IN-check
  echo $N" "$N" 48 128 1 0 2 1" >> IN-check
done
echo "-1 0 0 0 0 0 0 0" >> IN-check

INTEL_MPI="yes"
INTEL_MPI1=`ldd eigenexa_benchmark | awk '/libmpifort/{ print $1}'`
INTEL_MPI2=`which mpiexec.hydra | grep " no "`
if [ x$INTEL_MPI1 = x ]; then
	INTEL_MPI="no"
fi
if [ x$INTEL_MPI2 != x ]; then
	INTEL_MPI="no"
fi

if [ $INTEL_MPI = "yes" ]; then
	export I_MPI_CBWR=2
	export FI_SOCKETS_IFACE=enp4s0f0
	export FI_PROVIDER=sockets
	export I_MPI_FABRICS=shm
fi

\rm LOG-$P-$T
  
for P in 1600; do
for T in 12; do
  export OMP_NUM_THREADS=$T
  if [ $INTEL_MPI = "yes" ]; then
# intel MPI
    mpiexec.hydra \
	-np $P -genv OMP_NUM_THREADS $T \
	./eigenexa_benchmark -f IN-check < /dev/null |& tee LOG-$P-$T
  else
# openMPI
    mpirun \
	-np $P -x OMP_NUM_THREADS=$T \
	./eigenexa_benchmark -f IN-check < /dev/null |& tee LOG-$P-$T

  fi
done
done

killall -9 mpiexec.hydra

