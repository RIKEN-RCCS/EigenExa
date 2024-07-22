#!/bin/bash

if [ -f IN-check ]; then
  \rm IN-check
fi

echo "!   N  nvec bx  by m t s e" > IN-check
awk 'BEGIN{ for(N=3;N<=256;N++){ \
	print N" "N" 48 128 1 0 0 1"; \
	print N" "N" 48 128 1 0 1 1"; \
	print N" "N" 48 128 1 2 0 1"; \
	print N" "N" 48 128 1 2 1 1"; \
    } exit;}END{}' >> IN-check
for N in 511 512 513 1023 1024 1025; do
  echo $N" "$N" 48 128 1 0 0 1" >> IN-check
  echo $N" "$N" 48 128 1 0 1 1" >> IN-check
  echo $N" "$N" 48 128 1 2 0 1" >> IN-check
  echo $N" "$N" 48 128 1 2 1 1" >> IN-check
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
  
for P in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24; do
for T in 1 2 3 4; do
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

