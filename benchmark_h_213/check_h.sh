#!/bin/sh


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

for P in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24; do
for T in 1 2 3 4; do
echo "P="$P" T="$T
export OMP_NUM_THREADS=$T
if [ $INTEL_MPI = "yes" ]; then
mpiexec.hydra \
        -np $P -genv OMP_NUM_THREADS $T \
        ./eigenexah_benchmark < /dev/null |& tee LOG-$P-$T
else
mpirun \
	-np $P -x OMP_NUM_THREADS=$T \
	./eigenexah_benchmark < /dev/null |& tee LOG-$P-$T
fi
done
done

