#!/bin/sh
#PJM -L "rscgrp=ea"
#PJM -L "node=4"
#PJM --mpi "proc=160"
#PJM -L "elapse=100"
#PJM -j
#PJM -g a30102

\rm LOG-*-*
\rm -rf log_json
mkdir log_json
P=${PJM_MPI_PROC}
T=1
echo "P="$P" T="$T
export OMP_NUM_THREADS=$T
mpiexec.hydra \
        -np $P \
        ./eigenexa_benchmark < /dev/null |& tee LOG-$P-$T