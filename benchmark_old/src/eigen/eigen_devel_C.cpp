#include "eigen_devel.hpp"
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace eigen_devel {
  extern "C" void set_TRD_parameters_C(int inod, int nnod, MPI_Fint fcomm) {
    TRD_COMM_WORLD = MPI_Comm_f2c(fcomm);
    TRD_inod = static_cast<int64_t>(inod);
    TRD_nnod = static_cast<int64_t>(nnod);
  }

#ifdef _OPENMP
  extern "C" void set_MPI_THREAD_MODE_C(int fcomm) {
//    MPI_THREAD_MODE = MPI_Comm_f2c(fcomm);
    MPI_THREAD_MODE = fcomm;
  }
#endif

  extern "C" void set_nod_info_C(
             int x_inod_, int x_nnod_, MPI_Fint xfcomm,
             int y_inod_, int y_nnod_, MPI_Fint yfcomm,
             int z_inod_, int z_nnod_, MPI_Fint zfcomm,
             int w_inod_, int w_nnod_, MPI_Fint wfcomm) {
    x_inod = static_cast<int64_t>(x_inod_);
    y_inod = static_cast<int64_t>(y_inod_);
    z_inod = static_cast<int64_t>(z_inod_);
    w_inod = static_cast<int64_t>(w_inod_);
    x_nnod = static_cast<int64_t>(x_nnod_);
    y_nnod = static_cast<int64_t>(y_nnod_);
    z_nnod = static_cast<int64_t>(z_nnod_);
    w_nnod = static_cast<int64_t>(w_nnod_);
    x_COMM_WORLD = MPI_Comm_f2c(xfcomm);
    y_COMM_WORLD = MPI_Comm_f2c(yfcomm);
    z_COMM_WORLD = MPI_Comm_f2c(zfcomm);
    w_COMM_WORLD = MPI_Comm_f2c(wfcomm);
  }

  extern "C" void set_repro_reduce_C(bool repro) {
    repro_reduce = repro;
  }

}
