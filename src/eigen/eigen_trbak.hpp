#pragma once
#ifndef EIGEN_TRBAK_HPP
#define EIGEN_TRBAK_HPP
#ifdef _OPENMP
#include <omp.h>
#endif
#include <memory>
namespace eigen_trbak {
inline double trbk_time_bcast = 0.0;
inline double trbk_time_reduc = 0.0;
inline double trbk_time_fr = 0.0;
inline double trbk_time_trbk1 = 0.0;
inline double trbk_time_trbk1_ = 0.0;
inline double trbk_time_trbk1x = 0.0;
inline double trbk_time_trbk1x_ = 0.0;
inline double trbk_time_trbk2 = 0.0;
inline double trbk_time_reduc_overhead_x = 0.0;
#ifdef _OPENMP
inline omp_lock_t trbk_lock;
inline int64_t trbk_mask[2];
#endif

inline int64_t do_overlap_bcast_level = 0;
inline int64_t trbk_time_counter = 0;
inline int64_t trbk_time_interval = 0;
inline int64_t trbk_time_next = 0;
inline int64_t trbk_switched = 0;

//inline double *trbk_buf = nullptr;
std::unique_ptr<double[]> trbk_buf = nullptr;

}
#endif
