#pragma once
#ifndef DC2_FS_HPP
#define DC2_FS_HPP

#include <mpi.h>

#include <cstdio>
#include <memory>

#include "../eigen/eigen_dc.hpp"
#include "../eigen/eigen_devel.hpp"
#include "../eigen/eigen_libs0.hpp"
#include "FS_EDC.hpp"
#include "FS_libs.hpp"

namespace eigen_FS {
using eigen_devel::eigen_abort;
using eigen_devel::eigen_timer_print;
using eigen_devel::eigen_timer_reset;
using eigen_libs0::eigen_get_comm;
using eigen_libs0::eigen_get_id;
using eigen_libs0::eigen_get_procs;
using FS_libs::FS_WorkSize;
using std::unique_ptr;

template <class Float>
void dc2_FS(int n, int nvec, Float d[], Float e[], Float z[], int ldz,
            int *info, Float *ret) {
  eigen_dc::flops = 0;
  eigen_dc::dgemm_time = 0;
  eigen_dc::p_time0 = 0;
  eigen_dc::p_timer = 0;
  eigen_dc::p_time2 = 0;
  eigen_dc::p_time3 = 0;
  eigen_dc::p_times = 0;
  eigen_dc::p_timez = 0;
  eigen_timer_reset(1, 0, 0, 0);

  int nprocs, nprow, npcol;
  eigen_get_procs(nprocs, nprow, npcol);
  int iam, myrow, mycol;
  eigen_get_id(iam, myrow, mycol);
  iam -= 1;
  myrow -= 1;
  mycol -= 1;

  MPI_Comm eigen_comm, eigen_x_comm, eigen_y_comm;
  eigen_get_comm(eigen_comm, eigen_x_comm, eigen_y_comm);
  int lwork_, liwork_;
  FS_WorkSize(n, lwork_, liwork_);

  const int FS_COMM_MEMBER = FS_libs::is_comm_member();
  if (!FS_COMM_MEMBER) {
    lwork_ = 0;
    liwork_ = 0;
  }
  int lwork;
  int liwork;
  MPI_Allreduce(&lwork_, &lwork, 1, MPI_INT, MPI_MAX, eigen_comm);
  MPI_Allreduce(&liwork_, &liwork, 1, MPI_INT, MPI_MAX, eigen_comm);

  try {
    unique_ptr<Float[]> work(new Float[lwork]);
    unique_ptr<int[]> iwork(new int[liwork]);

#if defined(__INTEL_COMPILER) && USE_MKL
    const auto mkl_mode = mkl_get_Dynamic();
    MKL_Set_Dynamic(0);
#endif

    FS_prof prof = {};
#if TIMER_PRINT

    prof.init();
#endif
    *info = FS_EDC::FS_EDC<int, Float>(n, d, e, z, ldz, work.get(), lwork,
                                       iwork.get(), liwork, &prof);
#if TIMER_PRINT
    prof.finalize();
#endif
#if TIMER_PRINT > 1
    eigen_dc::p_time0 = prof.region_time[21];
    eigen_dc::p_timer = prof.region_time[70];
    eigen_dc::p_time2 = prof.region_time[50];
    eigen_dc::p_time3 = prof.region_time[60];
    eigen_dc::dgemm_time = prof.region_time[67];
    eigen_dc::p_timez = prof.region_time[40];
#endif
#if defined(__INTEL_COMPILER) && USE_MKL
    MKL_Set_Dynamic(mkl_mode);
#endif

#if TIMER_PRINT > 1
    if (iam == 0) {
      printf("FS_EDC     %f\n", prof.region_time[10]);
    }
#endif

  } catch (std::bad_alloc) {
    eigen_abort();
  }

#if TIMER_PRINT > 1
  if (iam == 0) {
    printf("FS_dividing %f\n", eigen_dc::p_time0);
    printf("FS_pdlasrt  %f\n", eigen_dc::p_timer);
    printf("FS_pdlaed2  %f\n", eigen_dc::p_time2);
    printf("FS_pdlaed3  %f\n", eigen_dc::p_time3);
    printf("FS_pdlaedz  %f\n", eigen_dc::p_timez);
    printf("DGEMM       %f\n", eigen_dc::dgemm_time);
  }
#endif

  MPI_Allreduce(&eigen_dc::flops, ret, 1, FS_const::MPI_TYPE<double>, MPI_SUM,
                eigen_comm);

  return;
}
}  // namespace eigen_FS

#endif
