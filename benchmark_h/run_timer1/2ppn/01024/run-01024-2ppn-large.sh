#!/bin/sh
#PJM -L "node=32x32"
#PJM -L "rscunit=rscunit_ft01"
#PJM -L "rscgrp=large"
#PJM -L "elapse=24:00:00"
#PJM -L "freq=2200,eco_state=2"
#PJM --mpi "max-proc-per-node=2"
#PJM -s


export PLE_MPI_STD_EMPTYFILE=off
export OMP_NUM_THREADS=24

export NOD=32x32
export PE=2048
export PPN=2

export OMPI_MCA_plm_ple_memory_allocation_policy=interleave_local
export OMPI_MCA_plm_ple_numanode_assign_policy=simplex

printenv

date
/usr/bin/llio_transfer ./IN-Frank
/usr/bin/llio_transfer ./IN-Random
/usr/bin/llio_transfer ./eigenexa_benchmark
date
mpiexec         -stdout-proc ./output.%j/%/1000r/stdout         -stderr-proc ./output.%j/%/1000r/stderr         ./eigenexa_benchmark -f IN-Frank
date
mpiexec         -stdout-proc ./output.%j/%/1000r/stdout         -stderr-proc ./output.%j/%/1000r/stderr         ./eigenexa_benchmark -f IN-Random
date
/usr/bin/llio_transfer --purge ./eigenexa_benchmark
date

