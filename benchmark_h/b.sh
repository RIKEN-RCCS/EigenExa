#!/bin/sh
#------ pjsub option --------#
#PJM -L rscgrp=n22247a
#PJM -L node=1
#PJM -L elapse=1:00:00
#PJM -g n22247
#
#PJM -j
#------- Program execution -------#

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
\rm -rf log_json
mkdir log_json

for P in 20; do
for T in 4; do
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
fi
done
done

