#pragma once
#ifndef FS_PROF_HPP
#define FS_PROF_HPP
#include <mpi.h>

#include <algorithm>
#include <cstring>

#include "FS_libs.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace FS_prof {
constexpr int FS_max_region = 70 + 1;

class FS_prof {
 public:
  char region_name[FS_max_region][64];
  double region_time[FS_max_region];
  double region_start[FS_max_region];
  int region_ecount[FS_max_region];
#ifdef COUNT_CHECK
  int region_scount[FS_max_region];
#endif
 public:
  void init();
  void start(int id);
  void end(int id);
  void add(const FS_prof &prof_add);
  void finalize();
};

}  // namespace FS_prof

namespace FS_prof {
inline void FS_prof::init() {
  std::strncpy(region_name[1], "total", 32);
  std::strncpy(region_name[10], "FS_EDC", 32);
  std::strncpy(region_name[11], "  DSTEDC(NPROCS=1)", 32);
  std::strncpy(region_name[20], "  FS_PDLAED0", 32);
  std::strncpy(region_name[21], "    FS_dividing", 32);
  std::strncpy(region_name[22], "      FS_create_merge_comm", 32);
  std::strncpy(region_name[23], "      barrier", 32);
  std::strncpy(region_name[24], "      FS_create_mergeXY_group", 32);
  std::strncpy(region_name[25], "      FS_create_merge_comm_recursive", 32);
  std::strncpy(region_name[28], "    DSTEDC", 32);
  std::strncpy(region_name[30], "    FS_PDLAED1", 32);
  std::strncpy(region_name[31], "      FS_MERGE_D", 32);
  std::strncpy(region_name[40], "      FS_PDLAEDZ", 32);
  std::strncpy(region_name[45], "      FS_REDUCE_ZD", 32);
  std::strncpy(region_name[46], "        barrier", 32);
  std::strncpy(region_name[47], "        allreduce", 32);
  std::strncpy(region_name[50], "      FS_PDLAED2", 32);
  std::strncpy(region_name[60], "      FS_PDLAED3", 32);
  std::strncpy(region_name[61], "        DLAED4(Z,D)", 32);
  std::strncpy(region_name[62], "        allreduce x3", 32);
  std::strncpy(region_name[63], "        COPY Q2", 32);
  std::strncpy(region_name[64], "        WAIT/UNPACK", 32);
  std::strncpy(region_name[65], "        PACK/ISEND/IRECV", 32);
  // std::strncpy(region_name[66], "        DLAED4(DELTA)+DGEMM", 32);
  std::strncpy(region_name[66], "        RECALC DELTA+DGEMM", 32);
  std::strncpy(region_name[67], "          DGEMM", 32);
  std::strncpy(region_name[70], "    FS_PDLASRT", 32);
  std::fill_n(this->region_time, FS_max_region, 0);
  std::fill_n(this->region_start, FS_max_region, 0);
  std::fill_n(this->region_ecount, FS_max_region, 0);
#ifdef COUNT_CHECK
  std::fill_n(this->region_scount, FS_max_region, 0);
#endif
}
inline void FS_prof::start(int id) {
#ifdef _OPENMP
  const auto tt = omp_get_wtime();
#else
  const auto tt = MPI_Wtime();
#endif
  this->region_start[id] = tt;
#ifdef COUNT_CHECK
  this->region_scount[id] += 1;
#endif
  return;
}
inline void FS_prof::end(int id) {
#ifdef _OPENMP
  const auto tt = omp_get_wtime();
#else
  const auto tt = MPI_Wtime();
#endif
  this->region_time[id] += (tt - this->region_start[id]);
  this->region_ecount[id] += 1;
  return;
}
inline void FS_prof::add(const FS_prof &prof_add) {
  for (auto i = 0; i < FS_max_region; i++) {
    this->region_time[i] += prof_add.region_time[i];
    this->region_ecount[i] += prof_add.region_ecount[i];
#ifdef COUNT_CHECK
    this->region_scount[i] += prof_add.region_scount[i];
#endif
  }
}
inline void FS_prof::finalize() {
  int nprocs_out;
  MPI_Comm_size(FS_libs::FS_COMM_WORLD, &nprocs_out);
#ifndef TIMER_ALLPROCS
  nprocs_out = 1;
#endif

  // 集計と出力
  if (FS_libs::FS_get_myrank() == 0) {
    printf(
        " ================================================================\n");
#ifdef _BLOCKING_DGEMM
    printf("  TIMING INFO (DGEMM BLOCKING)\n");
#else
    printf("  TIMING INFO (DGEMM NON-BLOCKING)\n");
#endif

    // 各ランクから情報を受け取りながら出力
    double tmp_region_time[FS_max_region];
#ifdef COUNT_CHECK
    int tmp_region_scount[FS_max_region];
#endif
    int tmp_region_ecount[FS_max_region];
    MPI_Status stat;
    for (auto n = 0; n < nprocs_out; n++) {
      if (n == 0) {
        std::copy_n(this->region_time, FS_max_region, tmp_region_time);
#ifdef COUNT_CHECK
        std::copy_n(this->region_scount, FS_max_region, tmp_region_scount);
#endif
        std::copy_n(this->region_ecount, FS_max_region, tmp_region_ecount);
      } else {
        MPI_Recv(tmp_region_time, FS_max_region, MPI_DOUBLE, n, 0,
                 FS_libs::FS_COMM_WORLD, &stat);
#ifdef COUNT_CHECK
        MPI_Recv(tmp_region_scount, FS_max_region, MPI_INT, n, 0,
                 FS_libs::FS_COMM_WORLD, &stat);
#endif
        MPI_Recv(tmp_region_ecount, FS_max_region, MPI_INT, n, 0,
                 FS_libs::FS_COMM_WORLD, &stat);
      }

      printf("  RANK = %d\n", n);
      printf(
          " -ID-+----region name ----------------+----time [s] "
          "---------+-count-\n");
      for (auto i = 0; i < FS_max_region; i++) {
        if (tmp_region_ecount[i] > 0) {
          printf("%d %-30s %f %d\n", i, region_name[i], tmp_region_time[i],
                 tmp_region_ecount[i]);
#ifdef COUNT_CHECK
          if (tmp_region_scount[i] != tmp_region_ecount[i]) {
            printf("  Warning : start/end count are different in [%s]\n",
                   region_name[i]);
          }
#endif
        }
      }
      printf(
          " ==================================================================="
          "=\n");
    }
  } else {
    // FS_get_myrank() != 0
    // ランク0に情報を送信
#ifdef TIMER_ALLPROCS
    if (nprocs_out > 1) {
      MPI_Send(this->region_time, FS_max_region, MPI_DOUBLE, 0, 0,
               FS_libs::FS_COMM_WORLD);
#ifdef COUNT_CHECK
      MPI_Send(this->region_scount, FS_max_region, MPI_INT, 0, 0,
               FS_libs::FS_COMM_WORLD);
#endif
      MPI_Send(this->region_ecount, FS_max_region, MPI_INT, 0, 0,
               FS_libs::FS_COMM_WORLD);
    }
#endif
  }
}
}  // namespace FS_prof

#endif
