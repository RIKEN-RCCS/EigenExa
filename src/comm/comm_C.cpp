
#include "comm.hpp"
#include <mpi.h>

namespace comm {

  extern "C" void comm_barrier(MPI_Fint* fcomm) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    barrier(comm);
  }

  extern "C" void comm_bcast_dbl(double* buf, int* n, int* root, int* col_id, MPI_Fint* fcomm) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    bcast<double>(buf, static_cast<int64_t>(*n), *root, static_cast<int64_t>(*col_id), comm);
  }

  extern "C" void comm_bcastw_dbl(double* buf, int* n, int* root, int* lda, int* lpx,
                     double* work, int* col_id, MPI_Fint* fcomm) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    bcastw<double>(buf, static_cast<int64_t>(*n), static_cast<int64_t>(*root), static_cast<int64_t>(*lda), static_cast<int64_t>(*lpx), work, static_cast<int64_t>(*col_id), comm);
  }
  extern "C" void comm_reduce_dbl(double* buf, double* work, int* n, int* col_id, MPI_Fint* fcomm) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    reduce<double>(buf, work, static_cast<int64_t>(*n), static_cast<int64_t>(*col_id), comm);
  }

  extern "C" void comm_allgather_dbl(double* sendbuf, double* recvbuf, int* n, int* col_id, MPI_Fint* fcomm) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    allgather<double>(sendbuf, recvbuf, static_cast<int64_t>(*n), static_cast<int64_t>(*col_id), comm);
  }

}
