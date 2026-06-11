#pragma once
#ifndef EIGEN_BLACS_HPP
#define EIGEN_BLACS_HPP
#include <mpi.h>
#include <vector>
#include <string>
#include <stdexcept>

namespace eigen_blacs {

extern "C" {
  void blacs_get_c(int context, int what, int* ictxt);
  void blacs_gridmap_c(int ictxt, const int* usermap, int ldup, int nprow, int npcol);
  void blacs_gridexit_c(int ictxt);
}

inline int BLACS_ICONTXT_FOR_EIGENEXA = -1;

inline void eigen_blacs_init(MPI_Comm TRD_COMM_WORLD,
                             int x_nnod, int y_nnod,
                             char GRID_major) {
  int ierr = 0;

  // BLACS context を取得
  int zero = 0;
  blacs_get_c(0, 0, &BLACS_ICONTXT_FOR_EIGENEXA);

  // グリッドマップ用配列
  std::vector<int> kk0(x_nnod), kk1(x_nnod);
  std::vector<int> tmpgrid(x_nnod * y_nnod);

  MPI_Group group0, group1;
  MPI_Comm_group(MPI_COMM_WORLD, &group0);
  MPI_Comm_group(TRD_COMM_WORLD, &group1);

  if (GRID_major == 'R') {
    for (int j = 0; j < y_nnod; ++j) {
      for (int i = 0; i < x_nnod; ++i) {
        kk1[i] = j + i * y_nnod;
      }
      int k = x_nnod;
      MPI_Group_translate_ranks(group1, k, kk1.data(), group0, kk0.data());
      for (int i = 0; i < x_nnod; ++i) {
        tmpgrid[i + j * x_nnod] = kk0[i];
      }
    }
  } else {
    for (int j = 0; j < y_nnod; ++j) {
      for (int i = 0; i < x_nnod; ++i) {
        kk1[i] = i + j * x_nnod;
      }
      int k = x_nnod;
      MPI_Group_translate_ranks(group1, k, kk1.data(), group0, kk0.data());
      for (int i = 0; i < x_nnod; ++i) {
        tmpgrid[i + j * x_nnod] = kk0[i];
      }
    }
  }

  MPI_Barrier(TRD_COMM_WORLD);

  blacs_gridmap_c(BLACS_ICONTXT_FOR_EIGENEXA,
                 tmpgrid.data(), x_nnod, x_nnod, y_nnod);

  MPI_Barrier(TRD_COMM_WORLD);

  MPI_Group_free(&group0);
  MPI_Group_free(&group1);
}

inline void eigen_blacs_exit() {
  blacs_gridexit_c(BLACS_ICONTXT_FOR_EIGENEXA);
  // BLACS_EXIT(1); // optional, often omitted
}

inline int eigen_get_blacs_context() {
  return BLACS_ICONTXT_FOR_EIGENEXA;
}

}
#endif
