#pragma once
#ifndef EIGEN_LIBS_HPP
#define EIGEN_LIBS_HPP
#include <mpi.h>
#include <string>
#include <cstring>
#include <array>
#include <iostream>
#include <cmath>
#include <string>
#include <optional>
#include <thread>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cstdint>
#include <cctype>
#include "eigen_devel.hpp"
#include "eigen_blacs.hpp"
#include "../comm/comm.hpp"
#include "../CSTAB/CSTAB.hpp"

namespace eigen_libs0 {
using std::int64_t;
using eigen_devel::TRD_inod;
using eigen_devel::TRD_nnod;
using eigen_devel::x_nnod;
using eigen_devel::x_inod;
using eigen_devel::y_nnod;
using eigen_devel::y_inod;
using eigen_devel::z_nnod;
using eigen_devel::z_inod;
using eigen_devel::w_nnod;
using eigen_devel::w_inod;
using eigen_devel::eigen_get_initialized;
using eigen_devel::eigen_set_initialized;
using eigen_devel::eigen_unset_initialized;
using eigen_devel::eigen_timer_reset;
using eigen_devel::eigen_abort;
using eigen_devel::TRD_COMM_WORLD;
using eigen_devel::x_COMM_WORLD;
using eigen_devel::y_COMM_WORLD;
using eigen_devel::z_COMM_WORLD;
using eigen_devel::w_COMM_WORLD;
using eigen_devel::diag_0;
using eigen_devel::diag_1;
using eigen_devel::n_common;
using eigen_devel::p0_;
using eigen_devel::q0_;
using eigen_devel::repro_reduce;
using eigen_devel::barrier_overhead_x;
using eigen_devel::barrier_overhead_y;
using eigen_devel::reduce_overhead_x;
using eigen_devel::reduce_overhead_y;
using eigen_devel::reduce_cont_overhead_x;
using eigen_devel::reduce_cont_overhead_y;
using eigen_devel::bcast_overhead_x;
using eigen_devel::bcast_overhead_y;
using eigen_devel::bcast_cont_overhead_x;
using eigen_devel::bcast_cont_overhead_y;
using eigen_devel::eigen_get_wtime;
using eigen_devel::eigen_set_grid_major;
using comm::reduce;
using comm::bcast;
using eigen_blacs::eigen_blacs_init;
using CSTAB::get_optdim;

extern "C" {

namespace eigen_comm_int {
void eigen_init0(int comm, char order);
void eigen_get_comm(int &eigen_comm, int &eigen_x_comm, int &eigen_y_comm);
}  // namespace eigen_comm_int
void eigen_get_procs(int &procs, int &x_procs, int &y_procs);
void eigen_get_id(int &id, int &x_id, int &y_id);
//void eigen_free0();
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

inline int64_t eigen_translate_g2l(int ictr, int64_t nnod) { return ictr / nnod; }
inline int64_t eigen_owner_node(int ictr, int64_t nnod) { return ictr % nnod; }

inline constexpr const char* CodeNAME = "EigenExa";

struct version_t {
  int64_t Major_Version;
  int64_t Minor_Version;
  int64_t Patch_Level;
  std::string date;
  std::string vcode;
};

extern version_t Eigen_Version;
//inline version_t Eigen_Version = {
//  2, 12, 0,
//  "October 25, 2022",
//#ifdef CODE_AKASHI
//  "otome/akashi"
//#else
//  "otome"
//#endif
//};

inline constexpr int64_t eigen_NB   = 64;
inline constexpr int64_t eigen_NB_f = 48;
inline constexpr int64_t eigen_NB_b = 128;

inline void eigen_get_version(int& version,
                              std::string* date = nullptr,
                              std::string* vcode = nullptr) {
  version = Eigen_Version.Major_Version * 10000 +
            Eigen_Version.Minor_Version * 100 +
            Eigen_Version.Patch_Level;

  if (date)  *date  = Eigen_Version.date;
  if (vcode) *vcode = Eigen_Version.vcode;
}

inline void eigen_show_version() {
  const char patch_table[] = " abcdefghijklmnopqrstuvwxyz*";
  char patchlevel = patch_table[std::min(static_cast<int64_t>(26), Eigen_Version.Patch_Level) + 1];

  std::string version = std::to_string(Eigen_Version.Major_Version) + "." +
                        std::to_string(Eigen_Version.Minor_Version) +
                        patchlevel;

  if (TRD_inod == 1) {
    std::cout << "## EigenExa version (" << version
              << ") / (" << Eigen_Version.date
              << ") / (" << Eigen_Version.vcode << ")" << std::endl;
  }
}

inline void eigen_initialized(bool& flag) {
  eigen_get_initialized(flag);
}

inline void eigen_free0(std::optional<int> flag = std::nullopt) {
  int64_t ierr = 0;
  bool local_flag = false;
  eigen_get_initialized(local_flag);

  if (!local_flag) return;

  if (TRD_COMM_WORLD != MPI_COMM_NULL) {
    if (n_common > 1 && x_nnod != y_nnod) {
      if (z_COMM_WORLD != MPI_COMM_NULL) MPI_Comm_free(&z_COMM_WORLD);
      if (w_COMM_WORLD != MPI_COMM_NULL) MPI_Comm_free(&w_COMM_WORLD);
    }

    if (x_COMM_WORLD != MPI_COMM_NULL) MPI_Comm_free(&x_COMM_WORLD);
    if (y_COMM_WORLD != MPI_COMM_NULL) MPI_Comm_free(&y_COMM_WORLD);
    if (TRD_COMM_WORLD != MPI_COMM_NULL) MPI_Comm_free(&TRD_COMM_WORLD);

    if (x_nnod != y_nnod) {
      p0_.clear();
      q0_.clear();
    }

#ifdef EIGEN_TIMER_PRINT
    if (flag.has_value() && flag.value() == 1) {
      eigen_timer_print("EigenExa(finalized)");
    }
#endif

    eigen_timer_reset(0, 0, 0, 0);
  }

  TRD_COMM_WORLD = MPI_COMM_WORLD;
  eigen_unset_initialized();
}

inline void eigen_init_comm_setup(std::optional<MPI_Comm> comm_opt) {
  int ierr = 0;
  int flag = 0;
  MPI_Comm comm0 = comm_opt.value_or(MPI_COMM_WORLD);

  ierr = MPI_Comm_test_inter(comm0, &flag);
  if (flag || ierr != MPI_SUCCESS || comm0 == MPI_COMM_NULL) {
    std::cerr << "*************\n"
                 "** CAUTION **\n"
                 "*************\n"
                 "You are going to initialize EigenExa with\n"
                 "an invalid communicator.\n"
                 "EigenExa terminates this run.\n";
    std::cout.flush();
    std::this_thread::sleep_for(std::chrono::seconds(1));
    std::cout.flush();
    MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    return;
  }

  if (comm0 == MPI_COMM_NULL) {
    TRD_COMM_WORLD = MPI_COMM_NULL;
    x_COMM_WORLD = y_COMM_WORLD = z_COMM_WORLD = w_COMM_WORLD = MPI_COMM_NULL;
    TRD_nnod = TRD_inod = x_nnod = x_inod = y_nnod = y_inod = z_nnod = z_inod = w_nnod = w_inod = 0;
  } else {
    MPI_Comm_dup(comm0, &TRD_COMM_WORLD);
#if defined(__INTEL_COMPILER)
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "I_MPI_CBWR", "2");
    MPI_Comm_set_info(TRD_COMM_WORLD, info);
    MPI_Info_free(&info);
#endif
  }
}

inline void eigen_init_omp_setup() {
// #ifdef _OPENMP
//   if (TRD_COMM_WORLD != MPI_COMM_NULL) {
//     int64_t ierr = 0;
//     int64_t local_size = 1;
// 
// #pragma omp parallel
// #pragma omp master
//     {
//       local_size = omp_get_num_threads();
//     }
// 
//     int64_t th0[2] = {local_size, -local_size};
//     int64_t th1[2] = {0, 0};
// 
//     MPI_Allreduce(th0, th1, 2, MPI_INT, MPI_MAX, TRD_COMM_WORLD);
// 
//     int64_t j = th1[0] + th1[1];
//     if (j != 0) {
//       MPI_Barrier(TRD_COMM_WORLD);
//       std::cout.flush();
//       if (TRD_inod == 1) {
//         std::cerr << "*************\n"
//                      "** CAUTION **\n"
//                      "*************\n"
//                      "EigenExa supports only homogeneous thread setting!\n"
//                      "EigenExa terminates this run.\n";
//       }
// 
//       for (int64_t i = 0; i < 2; ++i) {
//         MPI_Barrier(TRD_COMM_WORLD);
//         std::cout.flush();
//       }
// 
//       MPI_Barrier(TRD_COMM_WORLD);
//       std::cout.flush();
//       std::this_thread::sleep_for(std::chrono::seconds(1));
//       std::cout.flush();
//       MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
//     }
//   }
// #endif
}

inline char GRID_major = 'C';

inline void eigen_init_cartesian_check(std::optional<char> order_opt) {
  int64_t ierr = 0;
  if (TRD_COMM_WORLD != MPI_COMM_NULL) {
    // ランクとサイズ取得
    int trd_n=static_cast<int>(TRD_nnod);
    int trd_i=static_cast<int>(TRD_inod);
    ierr = MPI_Comm_size(TRD_COMM_WORLD, &trd_n);
    ierr = MPI_Comm_rank(TRD_COMM_WORLD, &trd_i);
    TRD_nnod = trd_n;
    TRD_inod = trd_i;
    TRD_inod += 1;

    // トポロジー確認
    int topo_type = 0;
    ierr = MPI_Topo_test(TRD_COMM_WORLD, &topo_type);

    int cart_dim = 1;
    if (topo_type == MPI_CART) {
      ierr = MPI_Cartdim_get(TRD_COMM_WORLD, &cart_dim);
    }

    if (cart_dim == 2) {
      // Cartesian トポロジーあり
      int dims[2], coords[2], periods[2];
      ierr = MPI_Cart_get(TRD_COMM_WORLD, cart_dim, dims, periods, coords);
      x_nnod = dims[0];
      y_nnod = dims[1];
      x_inod = coords[0] + 1;
      y_inod = coords[1] + 1;

      GRID_major = 'R';
      if (order_opt && x_inod == 1 && y_inod == 1 &&
          (*order_opt == 'C' || *order_opt == 'c')) {
        std::cerr << "*************\n"
                     "** CAUTION **\n"
                     "*************\n"
                     "The MPI_CART you specified is based on R-major,\n"
                     "but you also specified C-major option.\n"
                     "EigenExa solve this conflict by taking R-major.\n";
      }
    } else {
      // Cartesian トポロジーなし → 自動分割
      x_nnod = static_cast<int>(std::sqrt(static_cast<double>(TRD_nnod)));
      while (x_nnod > 1 && TRD_nnod % x_nnod != 0) --x_nnod;
      y_nnod = TRD_nnod / x_nnod;

      GRID_major = order_opt.value_or('C');
      GRID_major = (GRID_major == 'R' || GRID_major == 'r') ? 'R' : 'C';
      eigen_set_grid_major(GRID_major);

      if (GRID_major == 'R') {
        x_inod = (TRD_inod - 1) / y_nnod + 1;
        y_inod = (TRD_inod - 1) % y_nnod + 1;
      } else {
        x_inod = (TRD_inod - 1) % x_nnod + 1;
        y_inod = (TRD_inod - 1) / x_nnod + 1;
      }
    }

    // x/y方向の通信グループ作成
    MPI_Comm_split(TRD_COMM_WORLD, y_inod, x_inod, &x_COMM_WORLD);
    MPI_Comm_split(TRD_COMM_WORLD, x_inod, y_inod, &y_COMM_WORLD);

    // datacast用の補助配列と通信グループ
    if (x_nnod != y_nnod) {
      int64_t n1 = std::max(x_nnod, y_nnod);
      int64_t n2 = std::min(x_nnod, y_nnod);
      while (n1 != n2) {
        int64_t n3 = n1 - n2;
        n1 = std::max(n2, n3);
        n2 = std::min(n2, n3);
      }
      n_common = n1;

      p0_.assign(std::max(x_nnod, y_nnod), -1);
      q0_.assign(std::max(x_nnod, y_nnod), -1);

      for (int64_t i = 1; i <= x_nnod; ++i) {
        if ((i - 1) % n_common == (y_inod - 1) % n_common) {
          int64_t delta = y_inod - i;
          if (delta >= 0) {
            for (int64_t j = 1; j <= x_nnod; ++j) {
              int64_t k = delta + (j - 1) * y_nnod;
              if (k % x_nnod == 0) {
                p0_[i - 1] = k / x_nnod + 1;
                q0_[i - 1] = j;
                break;
              }
            }
          } else {
            for (int64_t j = 1; j <= y_nnod; ++j) {
              int64_t k = -delta + (j - 1) * x_nnod;
              if (k % y_nnod == 0) {
                q0_[i - 1] = k / y_nnod + 1;
                p0_[i - 1] = j;
                break;
              }
            }
          }
        }
      }

      diag_0 = diag_1 = 0;
      for (int64_t i = 1; i <= y_nnod / n_common; ++i) {
        int64_t j = (i - 1) * y_nnod + y_inod;
        int64_t k = (j - 1) % x_nnod + 1;
        if (k == x_inod) {
          diag_0 = i;
          diag_1 = (j - 1) / x_nnod + 1;
          break;
        }
      }

      w_inod = (x_inod - 1) / n_common + 1;
      z_inod = (x_inod - 1) % n_common + 1;
      w_nnod = x_nnod / n_common;
      z_nnod = n_common;

      MPI_Comm_split(x_COMM_WORLD, w_inod - 1, z_inod - 1, &z_COMM_WORLD);
      MPI_Comm_split(x_COMM_WORLD, z_inod - 1, w_inod - 1, &w_COMM_WORLD);
    } else {
      n_common = x_nnod;
      p0_.clear();
      q0_.clear();
      diag_0 = (y_inod == x_inod) ? 1 : 0;
      diag_1 = diag_0;
      z_COMM_WORLD = MPI_COMM_SELF;
      w_COMM_WORLD = x_COMM_WORLD;
      z_inod = 1; z_nnod = 1;
      w_inod = x_inod; w_nnod = x_nnod;
    }

#if DEBUG
    if (TRD_inod == 1) std::cout << "Cartesian check is OK\n";
#endif

  } else {
    // 無効なコミュニケータの場合
    x_COMM_WORLD = y_COMM_WORLD = z_COMM_WORLD = w_COMM_WORLD = MPI_COMM_NULL;
    x_nnod = x_inod = y_nnod = y_inod = z_nnod = z_inod = w_nnod = w_inod = 0;
    GRID_major = 'C';
    eigen_set_grid_major(GRID_major);
  }
}

inline void eigen_init_blacs_setup() {
#if TIMER_PRINT > 1
  if (TRD_inod <= 1) {
    std::cout << "GRID major " << GRID_major << " is specified." << std::endl;
  }
#endif

  MPI_Comm comm0;

  if (TRD_COMM_WORLD == MPI_COMM_NULL) {
    comm0 = MPI_COMM_SELF;
    x_nnod = 1; x_inod = 1;
    y_nnod = 1; y_inod = 1;
  } else {
    comm0 = TRD_COMM_WORLD;
  }

  eigen_blacs_init(comm0, x_nnod, y_nnod, GRID_major);

  if (TRD_COMM_WORLD == MPI_COMM_NULL) {
    x_nnod = 0; x_inod = 0;
    y_nnod = 0; y_inod = 0;
  }
}

inline void repro_check() {
#ifdef __INTEL_COMPILER
  const char* val;
  int64_t tmp[2], res[2];
  int64_t ierr = 0;

  val = std::getenv("I_MPI_ADJUST_ALLREDUCE");
  if (val && (val[0] == '4' || val[0] == '6')) {
    repro_reduce = false;
  }

  val = std::getenv("I_MPI_ADJUST_REDUCE");
  if (val && (val[0] == '3' || val[0] == '4' || val[0] == '6')) {
    repro_reduce = false;
  }

  if (repro_reduce) {
    tmp[0] = 1;
    tmp[1] = -1;
  } else {
    tmp[0] = 0;
    tmp[1] = 0;
  }

  ierr = MPI_Allreduce(tmp, res, 2, MPI_INT, MPI_MAX, TRD_COMM_WORLD);

  if (res[0] * res[1] == -1) {
    repro_reduce = true;
  } else {
    repro_reduce = false;
  }
#endif
}


inline void eigen_init_collective_comms() {
  int64_t ierr = 0;
  repro_check();

  if (TRD_COMM_WORLD != MPI_COMM_NULL) {
    ierr = MPI_Barrier(TRD_COMM_WORLD);

    std::vector<double> buff1(2048, 0.0), buff2(2048, 0.0);
    double s1 = 0.0, s2 = 0.0, d1 = 0.0, d2 = 0.0;
    double buff[4] = {0.0, 0.0, 0.0, 0.0};

#pragma omp parallel private(d1, d2)
    {
#pragma omp master
      {
        // Barrier overhead X
        s1 = s2 = 0.0;
        for (int64_t i = 1; i <= 10; ++i) {
          ierr = MPI_Barrier(TRD_COMM_WORLD);
          d1 = eigen_get_wtime();
          ierr = MPI_Barrier(x_COMM_WORLD);
          d2 = eigen_get_wtime();
          if (i > 5) s1 += (d2 - d1);
        }

        // Barrier overhead Y
        for (int64_t i = 1; i <= 10; ++i) {
          ierr = MPI_Barrier(TRD_COMM_WORLD);
          d1 = eigen_get_wtime();
          ierr = MPI_Barrier(y_COMM_WORLD);
          d2 = eigen_get_wtime();
          if (i > 5) s2 += (d2 - d1);
        }

        buff[0] = s1 / 5;
        buff[1] = s2 / 5;
        ierr = MPI_Allreduce(buff, buff + 2, 2, MPI_DOUBLE, MPI_SUM, TRD_COMM_WORLD);
        barrier_overhead_x = buff[2] / TRD_nnod;
        barrier_overhead_y = buff[3] / TRD_nnod;
      }

#pragma omp barrier
#pragma omp master
      {
        // Reduce overhead X
        s1 = s2 = 0.0;
        for (int64_t i = 1; i <= 10; ++i) {
          reduce<double>(buff1.data(), buff2.data(), 1024, 0, y_COMM_WORLD);
          ierr = MPI_Barrier(TRD_COMM_WORLD);
          d1 = eigen_get_wtime();
          reduce<double>(buff1.data(), buff2.data(), 2048, 0, x_COMM_WORLD);
          ierr = MPI_Barrier(x_COMM_WORLD);
          d2 = eigen_get_wtime();
          if (i > 5) s1 += (d2 - d1);
        }

        // Reduce overhead Y
        for (int64_t i = 1; i <= 10; ++i) {
          reduce<double>(buff1.data(), buff2.data(), 1024, 0, x_COMM_WORLD);
          ierr = MPI_Barrier(TRD_COMM_WORLD);
          d1 = eigen_get_wtime();
          reduce<double>(buff1.data(), buff2.data(), 2048, 0, y_COMM_WORLD);
          ierr = MPI_Barrier(y_COMM_WORLD);
          d2 = eigen_get_wtime();
          if (i > 5) s2 += (d2 - d1);
        }

        buff[0] = s1 / 5;
        buff[1] = s2 / 5;
        ierr = MPI_Allreduce(buff, buff + 2, 2, MPI_DOUBLE, MPI_SUM, TRD_COMM_WORLD);
        reduce_overhead_x = buff[2] / TRD_nnod - barrier_overhead_x;
        reduce_overhead_y = buff[3] / TRD_nnod - barrier_overhead_y;
      }

#pragma omp barrier
#pragma omp master
      {
        // Reduce continuation overhead
        s1 = s2 = 0.0;
        for (int64_t i = 1; i <= 10; ++i) {
          reduce<double>(buff1.data(), buff2.data(), 1024, 0, y_COMM_WORLD);
          ierr = MPI_Barrier(TRD_COMM_WORLD);
          d1 = eigen_get_wtime();
          reduce<double>(buff1.data(), buff2.data(), 1024, 0, x_COMM_WORLD);
          reduce<double>(buff2.data(), buff1.data(), 1024, 0, x_COMM_WORLD);
          ierr = MPI_Barrier(x_COMM_WORLD);
          d2 = eigen_get_wtime();
          if (i > 5) s1 += (d2 - d1);
        }

        for (int64_t i = 1; i <= 10; ++i) {
          reduce<double>(buff1.data(), buff2.data(), 1024, 0, x_COMM_WORLD);
          ierr = MPI_Barrier(TRD_COMM_WORLD);
          d1 = eigen_get_wtime();
          reduce<double>(buff1.data(), buff2.data(), 1024, 0, y_COMM_WORLD);
          reduce<double>(buff2.data(), buff1.data(), 1024, 0, y_COMM_WORLD);
          ierr = MPI_Barrier(y_COMM_WORLD);
          d2 = eigen_get_wtime();
          if (i > 5) s2 += (d2 - d1);
        }

        buff[0] = s1 / 5;
        buff[1] = s2 / 5;
        ierr = MPI_Allreduce(buff, buff + 2, 2, MPI_DOUBLE, MPI_SUM, TRD_COMM_WORLD);
        reduce_cont_overhead_x = buff[2] / TRD_nnod - barrier_overhead_x - reduce_overhead_x;
        reduce_cont_overhead_y = buff[3] / TRD_nnod - barrier_overhead_y - reduce_overhead_y;
      }

#pragma omp barrier
#pragma omp master
      {
        // Bcast overhead
        s1 = s2 = 0.0;
        for (int64_t i = 1; i <= 10; ++i) {
          bcast<double>(buff1.data(), 1024, 1, 0, y_COMM_WORLD);
          ierr = MPI_Barrier(TRD_COMM_WORLD);
          d1 = eigen_get_wtime();
          bcast<double>(buff1.data(), 2048, 1, 0, x_COMM_WORLD);
          ierr = MPI_Barrier(x_COMM_WORLD);
          d2 = eigen_get_wtime();
          if (i > 5) s1 += (d2 - d1);
        }

        for (int64_t i = 1; i <= 10; ++i) {
          bcast<double>(buff1.data(), 1024, 1, 0, x_COMM_WORLD);
          ierr = MPI_Barrier(TRD_COMM_WORLD);
          d1 = eigen_get_wtime();
          bcast<double>(buff1.data(), 2048, 1, 0, y_COMM_WORLD);
          ierr = MPI_Barrier(y_COMM_WORLD);
          d2 = eigen_get_wtime();
          if (i > 5) s2 += (d2 - d1);
        }

        buff[0] = s1 / 5;
        buff[1] = s2 / 5;
        ierr = MPI_Allreduce(buff, buff + 2, 2, MPI_DOUBLE, MPI_SUM, TRD_COMM_WORLD);
        bcast_overhead_x = buff[2] / TRD_nnod - barrier_overhead_x;
        bcast_overhead_y = buff[3] / TRD_nnod - barrier_overhead_y;
      }

#pragma omp barrier
#pragma omp master
      {
        // Bcast continuation overhead
        s1 = s2 = 0.0;
        for (int64_t i = 1; i <= 10; ++i) {
          bcast<double>(buff1.data(), 1024, 1, 0, y_COMM_WORLD);
          ierr = MPI_Barrier(TRD_COMM_WORLD);
          d1 = eigen_get_wtime();
          bcast<double>(buff1.data(), 1024, 1, 0, x_COMM_WORLD);
          bcast<double>(buff2.data(), 1024, 1, 0, x_COMM_WORLD);
          ierr = MPI_Barrier(x_COMM_WORLD);
          d2 = eigen_get_wtime();
          if (i > 5) s1 += (d2 - d1);
        }

        for (int64_t i = 1; i <= 10; ++i) {
          bcast<double>(buff1.data(), 1024, 1, 0, x_COMM_WORLD);
          ierr = MPI_Barrier(TRD_COMM_WORLD);
          d1 = eigen_get_wtime();
          bcast<double>(buff1.data(), 1024, 1, 0, y_COMM_WORLD);
          bcast<double>(buff2.data(), 1024, 1, 0, y_COMM_WORLD);
          ierr = MPI_Barrier(y_COMM_WORLD);
          d2 = eigen_get_wtime();
          if (i > 5) s2 += (d2 - d1);
        }

        buff[0] = s1 / 5;
        buff[1] = s2 / 5;
        ierr = MPI_Allreduce(buff, buff + 2, 2, MPI_DOUBLE, MPI_SUM, TRD_COMM_WORLD);
        bcast_cont_overhead_x = buff[2] / TRD_nnod - barrier_overhead_x - bcast_overhead_x;
        bcast_cont_overhead_y = buff[3] / TRD_nnod - barrier_overhead_y - bcast_overhead_y;
      } // end master
    } // end parallel
  } // end if TRD_COMM_WORLD
}

inline void eigen_init0(std::optional<MPI_Comm> comm = std::nullopt,
                        std::optional<char> order = std::nullopt,
                        std::optional<int> scalapack_context = std::nullopt,
                        std::optional<std::vector<std::vector<int>>> gridmap = std::nullopt) {
  bool flag = false;

  eigen_timer_reset(0, 0, 0, 0);
  eigen_get_initialized(flag);

  if (flag) {
    if (TRD_inod == 1) {
      std::cout << "*************\n"
                   "** CAUTION **\n"
                   "*************\n"
                   "You are going to initialize EigenExa,\n"
                   "while EigenExa was not freed at last call.\n"
                   "EigenExa restarts again by itself.\n";
    }
    eigen_free0();
  }

  eigen_init_comm_setup(comm);
  eigen_init_omp_setup();
  eigen_init_cartesian_check(order);
  eigen_init_blacs_setup();
#if DEBUG
  if (TRD_inod == 1) std::cout << "BLACS done\n";
#endif
  eigen_init_collective_comms();
#if DEBUG
  if (TRD_inod == 1) std::cout << "Collective sampling done\n";
#endif

  eigen_set_initialized();
#if DEBUG
  if (TRD_inod == 1) std::cout << "initialization done\n";
#endif
}

inline int64_t eigen_loop_start(int64_t istart, int64_t nnod, int64_t inod) {
  int64_t ret = (istart + nnod - 1 - inod) / nnod + 1;
  return ret;
}

inline int64_t eigen_loop_start(int64_t istart, const std::string& pdir, std::optional<int> inod = std::nullopt) {
  char dir = std::tolower(pdir[0]);
  switch (dir) {
    case 'w': case 't':
      return eigen_loop_start(istart, TRD_nnod, inod.value_or(TRD_inod));
    case 'x': case 'r':
      return eigen_loop_start(istart, x_nnod, inod.value_or(x_inod));
    case 'y': case 'c':
      return eigen_loop_start(istart, y_nnod, inod.value_or(y_inod));
    default:
      return 0;
  }
}

inline int64_t eigen_loop_end(int64_t iend, int64_t nnod, int64_t inod) {
  return (iend + nnod - 0 - inod) / nnod + 0;
}

inline int64_t eigen_loop_end(int64_t iend, const std::string& pdir, std::optional<int> inod = std::nullopt) {
  char dir = std::tolower(pdir[0]);

  switch (dir) {
    case 'w': case 't':
      return eigen_loop_end(iend, TRD_nnod, inod.value_or(TRD_inod));
    case 'x': case 'r':
      return eigen_loop_end(iend, x_nnod, inod.value_or(x_inod));
    case 'y': case 'c':
      return eigen_loop_end(iend, y_nnod, inod.value_or(y_inod));
    default:
      return -1;
  }
}

inline int64_t eigen_translate_l2g(int64_t ictr, int64_t nnod, int64_t inod) {
  return (ictr - 1) * nnod + inod;
}

inline int64_t eigen_translate_l2g(int64_t ictr, const std::string& pdir, std::optional<int> inod = std::nullopt) {
  char dir = std::tolower(pdir[0]);

  switch (dir) {
    case 'w': case 't':
      return eigen_translate_l2g(ictr, TRD_nnod, inod.value_or(TRD_inod));
    case 'x': case 'r':
      return eigen_translate_l2g(ictr, x_nnod, inod.value_or(x_inod));
    case 'y': case 'c':
      return eigen_translate_l2g(ictr, y_nnod, inod.value_or(y_inod));
    default:
      return -1;
  }
}

inline int64_t eigen_translate_g2l(int64_t ictr, int64_t nnod, int64_t inod) {
  return  (ictr-1)/nnod+1;
}

inline int64_t eigen_translate_g2l(int64_t ictr, const std::string& pdir, std::optional<int> inod = std::nullopt) {
  char dir = std::tolower(pdir[0]);

  switch (dir) {
    case 'w': case 't':
      return eigen_translate_g2l(ictr, TRD_nnod, inod.value_or(TRD_inod));
    case 'x': case 'r':
      return eigen_translate_g2l(ictr, x_nnod, inod.value_or(x_inod));
    case 'y': case 'c':
      return eigen_translate_g2l(ictr, y_nnod, inod.value_or(y_inod));
    default:
      return -1;
  }
}

// 基本の所有ノード計算
inline int64_t eigen_owner_node(int64_t ictr, int64_t nnod, int64_t inod) {
    return (ictr - 1) % nnod + 1;
}

// pdir に応じて分岐するラッパー
inline int64_t eigen_owner_node(int64_t ictr, const std::string& pdir, std::optional<int> inod = std::nullopt) {
    char dir = std::tolower(pdir[0]);

    switch (dir) {
        case 'w': case 't':
            return eigen_owner_node(ictr, TRD_nnod, inod.value_or(TRD_inod));
        case 'x': case 'r':
            return eigen_owner_node(ictr, x_nnod, inod.value_or(x_inod));
        case 'y': case 'c':
            return eigen_owner_node(ictr, y_nnod, inod.value_or(y_inod));
        default:
            return -1;
    }
}

inline void eigen_get_matdims0(
    int64_t n,
    int& nx,
    int& ny,
    std::optional<int> m_forward = std::nullopt,
    std::optional<int> m_backward = std::nullopt,
    std::optional<char> mode = std::nullopt)
{
  if (n <= 0) {
    nx = -1;
    ny = -1;
    return;
  }

  char mode_ = mode.value_or('O');

  if (mode_ == 'M') {
    nx = (n - 1) / x_nnod + 1;
    ny = (n - 1) / y_nnod + 1;
    return;
  }

  if (mode_ == 'L') {
    nx = (( (n - 1) / x_nnod + 1 ) + 31) / 32 * 32;
    ny = (n - 1) / y_nnod + 1;
    return;
  }

  // mode == 'O'
  int64_t NPROW = x_nnod;
  int64_t NPCOL = y_nnod;

  int64_t n1 = (n - 1) / NPROW + 1;
  int64_t nm;
  get_optdim(n1, 6, 16 * 4, 16 * 4 * 2, nm);

  int64_t m_f = m_forward.value_or(eigen_NB_f);
  int64_t m_b = m_backward.value_or(eigen_NB_b);

#if 0
  int64_t NB = std::max({m_f, m_b, eigen_NB});
#else
  int64_t NB = std::max(m_b, eigen_NB);
#endif

  int64_t nmz = ((n - 1) / NPROW + 1);
  nmz = ((nmz - 1) / NB + 1) * NB + 1;
  int64_t nn = nmz;
  nmz = (n - 1) / NB + 1;
  nmz = ((nmz - 1) / NPROW + 1) * NB;
  nmz = std::max(nn, nmz);

  int64_t nmw = ((n - 1) / NPCOL + 1);
  nmw = ((nmw - 1) / NB + 1) * NB + 1;
  nn = nmw;
  nmw = (n - 1) / NB + 1;
  nmw = ((nmw - 1) / NPCOL + 1) * NB;
  nmw = std::max(nn, nmw);

  int64_t larray = std::max(nmz, nm) * nmw;

  nx = nm;
  ny = (larray - 1) / nm + 1;

#if 1
  NB = std::min(eigen_NB, n);
  int64_t lddz = (n - 1) / std::min(NPCOL, NPROW) + 1;
  lddz = ((lddz - 1) / NB + 1) * NB;
  int64_t nxx = (n - 1) / std::min(NPCOL, NPROW) + 1;

  std::int64_t LX1 = static_cast<std::int64_t>(lddz) * lddz;
  std::int64_t LX2 = static_cast<std::int64_t>(nxx) * nxx;

  if (LX1 >= (1LL << 31) || LX2 >= (1LL << 31)) {
    std::cout << "Warning :: oversized problem !!" << std::endl;
    larray = -1;
    nx = -1;
    ny = -1;
  }
#endif
}

}
#endif
