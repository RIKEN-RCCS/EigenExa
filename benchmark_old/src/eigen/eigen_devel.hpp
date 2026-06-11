#pragma once
#ifndef EIGEN_DEVEL_HPP
#define EIGEN_DEVEL_HPP
#include <complex>
#include <vector>
#include <string>
#include <mpi.h>
#include <iostream>
#include <thread>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <numeric>
#include <iomanip>
#include <cstdint>
#include <array>

namespace eigen_devel {

extern "C" void eigen_abort();

//extern "C" void eigen_timer_reset(int64_t bcast, int64_t reduce, int64_t redist, int64_t gather);
extern "C" void eigen_timer_reset(int bcast, int reduce, int redist, int gather);

extern "C" double eigen_timer_print();

inline constexpr double ZERO  = 0.0;
inline constexpr double HALF  = 0.5;
inline constexpr double ONE   = 1.0;
inline constexpr double TWO   = 2.0;
inline constexpr double THREE = 3.0;
inline constexpr double FOUR  = 4.0;
inline constexpr double FIVE  = 5.0;
inline constexpr double SIX   = 6.0;
inline constexpr double SEVEN = 7.0;
inline constexpr double EIGHT = 8.0;
inline constexpr double NINE  = 9.0;
inline constexpr double TEN   = 10.0;

inline constexpr double MHALF = -0.5;
inline constexpr double MONE  = -1.0;
inline constexpr double MTWO  = -2.0;

inline constexpr std::complex<double> ZEROZ  = {ZERO, ZERO};
inline constexpr std::complex<double> ONEZ   = {ONE,  ZERO};
inline constexpr std::complex<double> MONEZ  = {MONE, ZERO};
inline constexpr std::complex<double> IONEZ  = {ZERO, ONE};
inline constexpr std::complex<double> IMONEZ = {ZERO, MONE};
#ifdef _OPENMP
//inline int64_t MPI_THREAD_MODE = MPI_THREAD_SINGLE;
inline int MPI_THREAD_MODE = MPI_THREAD_SINGLE;
#endif

// === MPI グリッド情報とフラグ ===
inline int64_t TRD_inod       = 0;
inline int64_t TRD_nnod       = 0;
inline MPI_Comm TRD_COMM_WORLD = MPI_COMM_WORLD;
//extern MPI_Comm x_COMM_WORLD, y_COMM_WORLD, z_COMM_WORLD, w_COMM_WORLD;
inline MPI_Comm x_COMM_WORLD = MPI_COMM_NULL;
inline MPI_Comm y_COMM_WORLD = MPI_COMM_NULL;
inline MPI_Comm z_COMM_WORLD = MPI_COMM_NULL;
inline MPI_Comm w_COMM_WORLD = MPI_COMM_NULL;

inline int64_t x_inod = 0, x_nnod = 0;
inline int64_t y_inod = 0, y_nnod = 0;
inline int64_t z_inod = 0, z_nnod = 0;
inline int64_t w_inod = 0, w_nnod = 0;

inline int64_t n_common = 0, diag_0 = 0, diag_1 = 0;
inline int64_t ERROR_INFO = 0;

inline bool Eigen_initialized_flag = false;
inline bool repro_reduce = false;

inline char Process_Grid_Major = 'C'; // 'R' or 'C'

extern std::vector<int> p0_;
extern std::vector<int> q0_;

inline double barrier_overhead_x = 0.0;
inline double barrier_overhead_y = 0.0;
inline double reduce_overhead_x = 0.0;
inline double reduce_overhead_y = 0.0;
inline double bcast_overhead_x = 0.0;
inline double bcast_overhead_y = 0.0;
inline double reduce_cont_overhead_x = 0.0;
inline double reduce_cont_overhead_y = 0.0;
inline double bcast_cont_overhead_x = 0.0;
inline double bcast_cont_overhead_y = 0.0;

// 通信時間の累積（秒）
inline double comm_time_reduction   = 0.0;
inline double comm_time_dc          = 0.0;
inline double comm_time_backtrafo   = 0.0;

// タイマー用の一時変数
inline double timer_t1 = 0.0;
inline double timer_t2 = 0.0;

// キャッシュ・オーバーラップフラグ
inline bool flag_overlap = false;
inline bool flag_oncache = false;

// バックワードハウスホルダー変換用パラメータ
inline constexpr int64_t nsx = 480;
inline constexpr int64_t nsm = 256;
inline constexpr int64_t ns0 = nsm * nsm + 6;
inline constexpr int64_t MBAND = 2;
#ifdef _OPENMP
inline int64_t TRBK_MASK_FULL = 0;
#endif

// チェックポイント数（初期値は0）
inline int64_t items_bcast  = 0;
inline int64_t items_reduce = 0;
inline int64_t items_redist = 0;
inline int64_t items_gather = 0;

inline double time_bcast  = 0.0;
inline double time_reduce = 0.0;
inline double time_redist = 0.0;
inline double time_gather = 0.0;

extern std::vector<double> time_bcast_;
extern std::vector<double> time_reduce_;
extern std::vector<double> time_redist_;
extern std::vector<double> time_gather_;

extern std::vector<int64_t> counter_bcast_;
extern std::vector<int64_t> counter_reduce_;
extern std::vector<int64_t> counter_redist_;
extern std::vector<int64_t> counter_gather_;

extern std::vector<int64_t> messages_bcast_;
extern std::vector<int64_t> messages_reduce_;
extern std::vector<int64_t> messages_redist_;
extern std::vector<int64_t> messages_gather_;

// --- エラー終了 ---
//inline void eigen_abort(const std::string& message, int64_t code) {
//  std::cerr << message << std::endl;
//  std::cerr.flush();
//  std::this_thread::sleep_for(std::chrono::seconds(1));
//  std::cerr.flush();
//  int64_t ierr;
//  MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
//}

// --- 経過時間取得 ---
inline double eigen_get_wtime() {
#ifdef _OPENMP
  return omp_get_wtime();
#else
  return MPI_Wtime();
#endif
}

// --- 初期化フラグの設定 ---
inline void eigen_set_initialized() {
  Eigen_initialized_flag = true;
}

inline void eigen_unset_initialized() {
  Eigen_initialized_flag = false;
}

inline void eigen_get_initialized(bool& flag) {
  flag = Eigen_initialized_flag;
}

// --- グリッドメジャーの設定と取得 ---
inline void eigen_set_grid_major(char major) {
  Process_Grid_Major = major;
}

inline char eigen_get_grid_major() {
  return Process_Grid_Major;
}


inline void eigen_timer_reset(int64_t bcast = 0, int64_t reduce = 0, int64_t redist = 0, int64_t gather = 0) {
  items_bcast  = bcast;
  items_reduce = reduce;
  items_redist = redist;
  items_gather = gather;

#ifdef EIGEN_TIMER_PRINT
  int64_t total_bcast  = items_bcast + items_redist;
  int64_t total_gather = items_gather + items_redist;

  time_bcast  = ZERO;
  time_reduce = ZERO;
  time_redist = ZERO;
  time_gather = ZERO;

  time_bcast_.assign(total_bcast, ZERO);
  time_reduce_.assign(items_reduce, ZERO);
  time_redist_.assign(items_redist, ZERO);
  time_gather_.assign(total_gather, ZERO);

  counter_bcast_.assign(total_bcast, 0);
  counter_reduce_.assign(items_reduce, 0);
  counter_redist_.assign(items_redist, 0);
  counter_gather_.assign(total_gather, 0);

  messages_bcast_.assign(total_bcast, 0);
  messages_reduce_.assign(items_reduce, 0);
  messages_redist_.assign(items_redist, 0);
  messages_gather_.assign(total_gather, 0);
#endif
}

inline double eigen_timer_print(const std::string& message) {
  double total_time = 0.0;

#ifdef EIGEN_TIMER_PRINT
  if (TRD_inod == 1) {
    std::cout << "COMM_STAT / [ " << message << " ]\n";

    auto print_section = [](const std::string& label, double time,
                            const std::vector<double>& times,
                            const std::vector<int64_t>& counters,
                            const std::vector<int64_t>& messages,
                            int64_t start, int64_t end, const std::string& prefix) {
      int64_t total_bytes = std::accumulate(counters.begin() + start, counters.begin() + end, 0LL) * 8;
      double throughput = (total_bytes > 0 && time > 0.0) ? total_bytes / time * 1e-9 : 0.0;
      std::cout << std::setw(10) << label << " :: "
                << std::scientific << std::setprecision(16) << time
                << "  " << throughput << " [GB/s]\n";
      for (int64_t i = start; i < end; ++i) {
        std::cout << "          " << prefix << " "
                  << std::setw(25) << times[i]
                  << std::setw(14) << counters[i]
                  << std::setw(10) << messages[i] << "\n";
      }
    };

    print_section("BCAST", time_bcast, time_bcast_, counter_bcast_, messages_bcast_,
                  0, items_bcast, "//");
    print_section("BCAST", time_bcast, time_bcast_, counter_bcast_, messages_bcast_,
                  items_bcast, items_bcast + items_redist, ";;");

    print_section("REDUCE", time_reduce, time_reduce_, counter_reduce_, messages_reduce_,
                  0, items_reduce, "//");

    print_section("GATHER", time_gather, time_gather_, counter_gather_, messages_gather_,
                  0, items_gather, "//");
    print_section("GATHER", time_gather, time_gather_, counter_gather_, messages_gather_,
                  items_gather, items_gather + items_redist, ";;");

    std::cout << "   REDIST :: " << std::scientific << std::setprecision(16) << time_redist << "\n";
    for (int64_t i = 0; i < items_redist; ++i) {
      std::cout << "          // "
                << std::setw(25) << time_redist_[i]
                << std::setw(14) << counter_redist_[i]
                << std::setw(10) << messages_redist_[i] << "\n";
    }

    total_time = time_bcast + time_reduce + time_redist + time_gather;
    std::cout << "   Total  :: " << total_time << "\n";
  }

  items_bcast = items_reduce = items_redist = items_gather = 0;
#endif

  return total_time;
}

#ifdef _OPENMP
// 初期化
inline void sync_other_than_master_init(omp_lock_t& TRBK_lock, int64_t TRBK_mask[2]) {
  int64_t local_size = omp_get_num_threads();
  TRBK_MASK_FULL = 0;

  if (local_size > 64) {
    TRBK_MASK_FULL = local_size - 1;
  } else {
    for (int64_t i = 1; i < local_size; ++i) {
      TRBK_MASK_FULL |= (1LL << i);
    }
  }

  omp_init_lock(&TRBK_lock);

  omp_set_lock(&TRBK_lock);
  TRBK_mask[0] = 0;
  TRBK_mask[1] = 0;
  omp_unset_lock(&TRBK_lock);
}

// 終了処理
inline void sync_other_than_master_finalize(omp_lock_t& TRBK_lock) {
  omp_destroy_lock(&TRBK_lock);
}

// 同期本体
inline void sync_other_than_master(omp_lock_t& TRBK_lock, int64_t TRBK_mask[2]) {
  int64_t local_size = omp_get_num_threads();
  int64_t local_rank = omp_get_thread_num();

  if (local_size == 1 || local_rank == 0) return;

  int64_t T;

  // 参加登録
  omp_set_lock(&TRBK_lock);
  if (local_size > 64) {
    T = TRBK_mask[0] + 1;
  } else {
    T = TRBK_mask[0] | (1LL << local_rank);
  }
  if (T == TRBK_MASK_FULL) TRBK_mask[1] = T;
  TRBK_mask[0] = T;
  omp_unset_lock(&TRBK_lock);

  // 全スレッドの参加を待つ
  while (true) {
    omp_set_lock(&TRBK_lock);
    T = TRBK_mask[0];
    omp_unset_lock(&TRBK_lock);
    if (T == TRBK_MASK_FULL) break;
    std::this_thread::yield();
  }

  // 離脱処理
  omp_set_lock(&TRBK_lock);
  if (local_size > 64) {
    T = TRBK_mask[1] - 1;
  } else {
    T = TRBK_mask[1] & ~(1LL << local_rank);
  }
  if (T == 0) TRBK_mask[0] = 0;
  TRBK_mask[1] = T;
  omp_unset_lock(&TRBK_lock);

  // 全スレッドの離脱を待つ
  while (true) {
    omp_set_lock(&TRBK_lock);
    T = TRBK_mask[1];
    omp_unset_lock(&TRBK_lock);
    if (T == 0) break;
    std::this_thread::yield();
  }
}
#endif

}  // namespace eigen_devel

#endif
