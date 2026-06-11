#pragma once
#ifndef EIGEN_LIBS_HPP
#define EIGEN_LIBS_HPP
#include <mpi.h>
namespace eigen_libs0 {
extern "C" {

namespace eigen_comm_int {
void eigen_init0(int comm, char order);
void eigen_get_comm(int &eigen_comm, int &eigen_x_comm, int &eigen_y_comm);
}  // namespace eigen_comm_int
void eigen_get_procs(int &procs, int &x_procs, int &y_procs);
void eigen_get_id(int &id, int &x_id, int &y_id);
void eigen_free0();
void eigen_get_grid_major(char &major);
}

inline void eigen_get_comm(MPI_Comm &eigen_comm, MPI_Comm &eigen_x_comm,
                           MPI_Comm &eigen_y_comm) {
  int comm, x_comm, y_comm;
  eigen_comm_int::eigen_get_comm(comm, x_comm, y_comm);
  eigen_comm = MPI_Comm_f2c(comm);
  eigen_x_comm = MPI_Comm_f2c(x_comm);
  eigen_y_comm = MPI_Comm_f2c(y_comm);
}
inline void eigen_init0(MPI_Comm comm, char order){
  eigen_comm_int::eigen_init0(MPI_Comm_c2f(comm), order);
}

inline int eigen_translate_g2l(int ictr, int nnod) { return ictr / nnod; }
inline int eigen_owner_node(int ictr, int nnod) { return ictr % nnod; }

}  // namespace eigen_libs0
#endif
