#pragma once
#ifndef TRBAKWY4_HPP
#define TRBAKWY4_HPP
#include <mpi.h>
#include <algorithm>
#include <omp.h>
#include <iostream>
#include <exception>
#include <cstdint>
#include "../eigen/eigen_trbak.hpp"
#include "../eigen/eigen_devel.hpp"
#include "../eigen/eigen_libs0.hpp"
#include "../comm/comm.hpp"
#include "../blas.hpp"

namespace trbakwy4 {
using std::int64_t;
using eigen_trbak::trbk_time_bcast;
using eigen_trbak::trbk_time_reduc;
using eigen_trbak::trbk_time_fr;
using eigen_trbak::trbk_time_trbk1;
using eigen_trbak::trbk_time_trbk1_;
using eigen_trbak::trbk_time_trbk1x;
using eigen_trbak::trbk_time_trbk1x_;
using eigen_trbak::trbk_time_trbk2;
using eigen_trbak::trbk_time_reduc_overhead_x;
using eigen_trbak::trbk_switched;
using eigen_trbak::do_overlap_bcast_level;
using eigen_trbak::trbk_time_counter;
using eigen_trbak::trbk_time_interval;
using eigen_trbak::trbk_time_next;
using eigen_trbak::trbk_buf;
using eigen_trbak::trbk_lock;
using eigen_trbak::trbk_mask;
using eigen_devel::TRD_COMM_WORLD;
using eigen_devel::x_COMM_WORLD;
using eigen_devel::y_COMM_WORLD;
using eigen_devel::TRD_nnod;
using eigen_devel::TRD_inod;
using eigen_devel::nsm;
using eigen_devel::ns0;
using eigen_devel::x_nnod;
using eigen_devel::x_inod;
using eigen_devel::y_nnod;
using eigen_devel::y_inod;
using eigen_devel::ONE;
using eigen_devel::ZERO;
using eigen_devel::eigen_get_wtime;
using eigen_devel::sync_other_than_master_init;
using eigen_devel::sync_other_than_master;
using eigen_devel::sync_other_than_master_finalize;
using eigen_devel::reduce_cont_overhead_x;
using eigen_devel::comm_time_backtrafo;
using eigen_devel::eigen_timer_reset;
using eigen_devel::eigen_timer_print;
using eigen_libs0::eigen_loop_start;
using eigen_libs0::eigen_loop_end;
using eigen_libs0::eigen_owner_node;
using eigen_libs0::eigen_translate_g2l;
using eigen_libs0::eigen_get_id;
using comm::reduce;
using comm::bcast;
using comm::bcastw;
using comm::allgather;
using comm::barrier;
using lapacke::gemm;
using lapacke::trsm;
using lapacke::dot;
using lapacke::gemv;
using lapacke::scal;
using lapacke::axpy;
using lapacke::ger;
using CSTAB::get_optdim;
using CSTAB::adjust_base;
using CSTAB::round_offset;
using CSTAB::n_columns;
using CSTAB::L1_WINDOW;
using CSTAB::L1_LSIZE;
using CSTAB::L2_LSIZE;

#ifndef USE_BCASTW
#define USE_BCASTW 1
#endif

void trbk_decide_overlap_level(int64_t i)
{
    int64_t local_size = omp_get_num_threads();
    int64_t local_rank = omp_get_thread_num();

    if (local_size > 1) {
        trbk_time_counter += 1;

        if (trbk_time_counter >= trbk_time_next) {
            double f0 = 1.0*(local_size - 1) / local_size;
            double f1 = 1.0 / f0;

            std::array<double, 3> bcast_time{};

            if (do_overlap_bcast_level == 2) {
                bcast_time[0] = trbk_time_bcast + trbk_time_trbk2 * f0 + trbk_time_trbk1 * f0;
                bcast_time[1] = std::max(trbk_time_bcast, trbk_time_trbk2) + trbk_time_trbk1 * f0;
                bcast_time[2] = std::max(trbk_time_bcast, trbk_time_trbk2 + trbk_time_trbk1);
            } else if (do_overlap_bcast_level == 1) {
                bcast_time[0] = trbk_time_bcast + trbk_time_trbk2 * f0 + trbk_time_trbk1;
                bcast_time[1] = std::max(trbk_time_bcast, trbk_time_trbk2) + trbk_time_trbk1;
                bcast_time[2] = std::max(trbk_time_bcast, trbk_time_trbk2 + trbk_time_trbk1 * f1);
            } else {
                bcast_time[0] = trbk_time_bcast + trbk_time_trbk2 + trbk_time_trbk1;
                bcast_time[1] = std::max(trbk_time_bcast, trbk_time_trbk2 * f1) + trbk_time_trbk1;
                bcast_time[2] = std::max(trbk_time_bcast, trbk_time_trbk2 * f1 + trbk_time_trbk1 * f1);
            }

            int64_t ll0 = do_overlap_bcast_level;
            int64_t ll1 = 0;
            for (int64_t ll2 = 1; ll2 <= 2; ++ll2) {
                if (bcast_time[ll1] >= bcast_time[ll2]) {
                    ll1 = ll2;
                }
            }

            std::array<int64_t, 6> ll{};
            int64_t ierr = 0;

            if (y_nnod >= 1024) {
                ll[ll1] = 1;
                MPI_Allreduce(&ll[0], &ll[3], 3, MPI_INT, MPI_SUM, y_COMM_WORLD);
            } else {
                int64_t packed = ll[0] + (ll[1] + ll[2] * 1024) * 1024;
                MPI_Allreduce(&packed, &ll[3], 1, MPI_INT, MPI_SUM, y_COMM_WORLD);
                int64_t ll1_ = ll[3];
                ll[3] = ll1_ % 1024;
                ll1_ /= 1024;
                ll[4] = ll1_ % 1024;
                ll[5] = ll1_ / 1024;
            }

            ll1 = 0;
#if OVERLAP_DECISION_TYPE == 3
            for (int64_t ll2 = 1; ll2 <= 2; ++ll2) {
                if (ll[3 + ll2] > ll[3 + ll1]) {
                    ll1 = ll2;
                }
            }
#elif OVERLAP_DECISION_TYPE == 2
            for (int64_t ll2 = 0; ll2 <= 2; ++ll2) {
                if (ll[3 + ll2] > 0) {
                    ll1 = ll2;
                }
            }
#elif OVERLAP_DECISION_TYPE == 1
            for (int64_t ll2 = 2; ll2 >= 0; --ll2) {
                if (ll[3 + ll2] > 0) {
                    ll1 = ll2;
                }
            }
#endif

#if _DEBUG_
            if (trd_inod == 1) {
                std::cout << i << " Overlap decision " << ll0 << " -> " << ll1
                          << "  BCAST_TIME = [" << bcast_time[0] << ", "
                          << bcast_time[1] << ", " << bcast_time[2] << "]\n";
            }
#endif

            if (ll1 == 0) {
                trbk_switched += 1;
            }

            if (ll0 != ll1) {
                trbk_time_interval = (trbk_switched >= 8) ? 32 : 1;
            } else {
                trbk_time_interval *= 2;
            }

            trbk_time_next += trbk_time_interval;
            do_overlap_bcast_level = ll1;
        }
    } else {
        do_overlap_bcast_level = 0;
    }
}

template <typename T>
void eigen_trbakwy_block_body1(
    const T* z, int64_t nmz,
    const T* v, int64_t nm, int64_t m,
    T* ss, T* sm,
    int64_t i_2, int64_t i_3, int64_t j_2, int64_t j_3)
{
    constexpr int64_t blas_chunk1 = 16;
    constexpr int64_t blas_chunk2 = 16;
    constexpr int64_t blas_chunk3 = 64;

    int64_t local_size = 1;
    int64_t local_rank = 0;
#ifdef _OPENMP
    local_size = omp_get_num_threads();
    local_rank = omp_get_thread_num();
#endif

    int64_t ll_size = local_size;
    int64_t ll_rank = local_rank;

#ifdef _OPENMP
    if ((do_overlap_bcast_level <= 1 && i_2 == 1) ||
        local_size == 1 || local_rank >= 1) {
        if ((do_overlap_bcast_level == 2 && local_size > 1) ||
            (i_2 > 1 && local_size > 1)) {
            ll_size = local_size - 1;
            ll_rank = local_rank - 1;
        } else {
            ll_size = local_size;
            ll_rank = local_rank;
	}
#else
            ll_size = local_size;
            ll_rank = local_rank;
#endif

        if (i_2 == 1) {
            int64_t j_5 = j_3 - j_2 + 1;
            if (j_5 > 0) {
                int64_t ii_step = 0;
                for (int64_t i_1 = 0; i_1 < m; i_1 += blas_chunk1) {
                    for (int64_t j_1 = 0; j_1 < m; j_1 += blas_chunk2) {
                        int64_t ii_3 = std::min(m, i_1 + blas_chunk1-1);
                        int64_t jj_3 = std::min(m, j_1 + blas_chunk2-1);
                        int64_t blk_size1 = ii_3 - i_1 + 1;
                        int64_t blk_size2 = jj_3 - j_1 + 1;
                        if (blk_size1 > 0 && blk_size2 > 0 ) {
                            if (i_1 + blk_size1 >= j_1 + blk_size2) {
                                int64_t ii_2 = i_1;
                                if (ii_step % ll_size == ll_rank) {
                                    gemm<T>(CblasTrans, CblasNoTrans,
                                        blk_size1, blk_size2, j_5,
                                        static_cast<T>(-1),
                                        &v[(j_2 - 1) + ii_2 * nm], nm,
                                        &v[(j_2 - 1) + j_1 * nm], nm,
                                        static_cast<T>(0),
                                        &sm[ii_2 + j_1 * nsm], nsm);
                                }
                                ++ii_step;
                            }
                        }
                    }
                }
            }
        }

#if 1 || defined(__INTEL_COMPILER)
        int64_t j_5 = j_3 - j_2 + 1;
        int64_t ii_step = 0;

        for (int64_t ii_2 = i_2; ii_2 <= i_3; ii_2 += blas_chunk3) {
            int64_t i_5 = std::min(blas_chunk3, i_3 - ii_2 + 1);
            if (m > 0 && j_5 > 0 && i_5 > 0) {
                if (ii_step % ll_size == ll_rank) {
                    gemm<T>(CblasTrans, CblasNoTrans,
                        m, i_5, j_5,
                        static_cast<T>(1),
                        &v[(j_2 - 1)], nm,
                        &z[(j_2 - 1) + (ii_2 - 1) * nmz], nmz,
                        static_cast<T>(0),
                        &ss[(ii_2 - 1) * m], m);
                }
                ++ii_step;
            }
        }
#else
        int64_t j_5 = j_3 - j_2 + 1;
        int64_t i_5 = (i_3 - i_2) / ll_size + 1;
        if (i_5 % 2 != 0) ++i_5;
        int64_t i_4 = i_5 * ll_rank;
        i_5 = std::min(i_5, i_3 - (i_2 + i_4) + 1);

        if (m > 0 && j_5 > 0 && i_5 > 0) {
            gemm<T>(CblasTrans, CblasNoTrans,
                    m, i_5, j_5,
                    static_cast<T>(1),
                    &v[(j_2 - 1)], nm,
                    &z[(j_2 - 1) + (i_2 + i_4 - 1) * nmz], nmz,
                    static_cast<T>(1),
                    &ss[(i_2 + i_4 - 1) * m], m);
        }
#endif
#ifdef _OPENMP
    }
#endif
}

template <typename T>
void eigen_trbakwy_block_body2(
    T* z, int64_t nmz,
    T* v, int64_t nm, int64_t m,
    T* ss, T* sm,
    int64_t i_2, int64_t i_3, int64_t j_2, int64_t j_3)
{
#if defined(__FUJITSU)
    constexpr int64_t chunk_m = 96;
    constexpr int64_t chunk_n = 40;
#else
    constexpr int64_t chunk_m = 128;
    constexpr int64_t chunk_n = 48;
#endif

    int64_t local_size = 1;
    int64_t local_rank = 0;
#ifdef _OPENMP
    local_size = omp_get_num_threads();
    local_rank = omp_get_thread_num();
#endif

    int64_t j_5 = (j_3 - j_2) / local_size + 1;
    j_5 = ((j_5 - 1) / 16 + 1) * 16;
    int64_t j_4 = j_5 * local_rank;
    j_5 = std::min(j_5, j_3 - (j_2 + j_4) + 1);

    if (j_5 > 0) {
        trsm<T>(CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasNonUnit,
                j_5, m, static_cast<T>(1),
                &sm[0], nsm,
                &v[(j_2 - 1 + j_4)], nm);
    }

#pragma omp barrier

    int64_t ll_size = local_size;
    int64_t ll_rank = local_rank;
#ifdef _OPENMP
    if (do_overlap_bcast_level == 0 || local_size == 1 || local_rank >= 1) {
      if (do_overlap_bcast_level >= 1 && local_size > 1) {
        ll_size = local_size - 1;
        ll_rank = local_rank - 1;
      } else {
        ll_size = local_size;
        ll_rank = local_rank;
      }
#else
        ll_size = local_size;
        ll_rank = local_rank;
#endif
    if (ll_rank >= 0) {
        int64_t ii_step = 0;
        for (int64_t ii_0 = i_2; ii_0 <= i_3; ii_0 += chunk_n) {
            int64_t ii_2 = ii_0;
            int64_t ii_3 = std::min(ii_0 + chunk_n - 1, i_3);

            for (int64_t jj_0 = j_2; jj_0 <= j_3; jj_0 += chunk_m) {
                int64_t jj_2 = jj_0;
                int64_t jj_3 = std::min(jj_0 + chunk_m - 1, j_3);

                int64_t i_blk = ii_3 - ii_2 + 1;
                int64_t j_blk = jj_3 - jj_2 + 1;

                if (m > 0 && j_blk > 0 && i_blk > 0) {
                    if (ii_step % ll_size == ll_rank) {
                        gemm<T>(CblasNoTrans, CblasNoTrans,
                                j_blk, i_blk, m,
                                static_cast<T>(1),
                                &v[(jj_2 - 1)], nm,
                                &ss[(ii_2 - 1) * m], m,
                                static_cast<T>(1),
                                &z[(jj_2 - 1) + (ii_2 - 1) * nmz], nmz);
                    }
                    ++ii_step;
                }
            }
        }
    }
#ifdef _OPENMP
    }
#endif
}

template <typename T>
void eigen_trbakwy_block_body(
    int64_t local_nvec,
    T* z, int64_t nmz,
    T* v, int64_t nm, int64_t m, int64_t i,
    T* ss, T* tt, int64_t nss, int64_t iblk,
    T& dcom, T& dx, T& dy, T& dz)
{
    int64_t i_1, i_2, i_3, i_4;
    int64_t j_2, j_3;
    int64_t i_0, m_0;
    int64_t ii, jj;
    double fr;
    double ds, de;

    int64_t local_rank = 0, local_size = 1;
#ifdef _OPENMP
    local_size = omp_get_num_threads();
    local_rank = omp_get_thread_num();
#endif
    int id, xid, yid;

    i_2 = 1;
    i_3 = local_nvec;

    j_2 = eigen_loop_start(1, "X");
    j_3 = eigen_loop_end(i + m - 1 - iblk, "X");

#ifdef _OPENMP
#if AT_BCAST_OVERLAP
    if (do_overlap_bcast_level == 2) {
        sync_other_than_master(trbk_lock, trbk_mask);
    }
#endif

#if !defined(AT_BCAST_OVERLAP) || (defined(AT_BCAST_OVERLAP) && do_overlap_bcast_level != 2)
#pragma omp barrier
#endif
#endif

    jj = i_3 - i_2 + 1;

#if AT_REDUCE_OVERLAP
    fr = trbk_time_fr;

    if (TRD_nnod > 1 &&
        trbk_time_reduc_overhead_x < trbk_time_trbk1 + trbk_time_trbk1x) {
        ii = std::max(8, static_cast<int>(jj * fr));
        ii = ((ii - 1) / 13 + 1) * 13;
        ii = std::min(jj, ii);
        fr = (trbk_time_trbk1 + trbk_time_trbk1x) /
             (trbk_time_trbk1 + trbk_time_trbk1x + trbk_time_reduc);
    } else {
        fr = ZERO;
        ii = jj;
    }

#if TIMER_PRINT
#pragma omp master
    if (trd_inod == 1) {
        printf("Overlap[on,off] %d %d %f %f %f %f\n",
               ii, (i_3 - i_2 + 1) - ii,
               trbk_time_trbk1 + trbk_time_trbk1x,
               trbk_time_reduc,
               (trbk_time_trbk1 + trbk_time_trbk1x) * (jj - ii) / jj,
               trbk_time_reduc * ii / jj);
    }
#endif

#else
    fr = ZERO;
    ii = i_3 - i_2 + 1;
#endif

    if (local_size == 1 || local_rank == 1) {
        ds = eigen_get_wtime();
    }

    eigen_trbakwy_block_body1<T>(z, nmz, v, nm, m,
                               &ss[ns0], &ss[0],
                               i_2, i_2 + ii - 1, j_2, j_3);

    if (local_size == 1 || local_rank == 1) {
        de = eigen_get_wtime();
        trbk_time_trbk1_ = de - ds;
#if TIMER_PRINT
        dx += (de - ds);
#endif
    }

#pragma omp barrier
#pragma omp master
    {
        ds = eigen_get_wtime();
        ss[ns0 - 6] = trbk_time_reduc;
        ss[ns0 - 5] = trbk_time_bcast;
        ss[ns0 - 4] = trbk_time_trbk1_;
        ss[ns0 - 3] = trbk_time_trbk1x_;
        ss[ns0 - 2] = trbk_time_trbk2;
        ss[ns0 - 1] = fr;

        reduce<T>(&ss[0], &tt[0], ns0 + ii * m, 3, x_COMM_WORLD);

        trbk_time_reduc  = ss[ns0 - 6] / x_nnod;
        trbk_time_bcast  = ss[ns0 - 5] / x_nnod;
        trbk_time_trbk1  = ss[ns0 - 4] / x_nnod;
        trbk_time_trbk1x = ss[ns0 - 3] / x_nnod;
        trbk_time_trbk2  = ss[ns0 - 2] / x_nnod;
        trbk_time_fr     = ss[ns0 - 1] / x_nnod;

        de = eigen_get_wtime();
        trbk_time_reduc = (de-ds);
        dcom += (de - ds);
    }

#if AT_REDUCE_OVERLAP
    if (i_3 >= i_2 + ii) {
#if TIMER_PRINT
        if (local_size == 1 || local_rank == 1) {
            ds = eigen_get_wtime();
        }
#endif
        eigen_trbakwy_block_body1<T>(z, nmz, v, nm, m,
                                   &ss[ns0], &ss[0],
                                   i_2 + ii, i_3, j_2, j_3);
#if TIMER_PRINT
        if (local_size == 1 || local_rank == 1) {
            de = eigen_get_wtime();
            dx += (de - ds);
        }
#endif
    }
#endif

#pragma omp barrier

#if AT_REDUCE_OVERLAP
    if (local_size == 1 || local_rank == 1) {
        if (i_3 >= i_2 + ii) {
            trbk_time_trbk1x_ = de - ds;
        } else {
            trbk_time_trbk1x_ = ZERO;
        }
    }

#pragma omp master
    {
        ds = eigen_get_wtime();
        reduce<T>(&ss[ns0 + ii * m], &tt[ns0 + ii * m],
                     (i_3 - i_2 + 1 - ii) * m, 4, x_COMM_WORLD);
        de = eigen_get_wtime();
        trbk_time_reduc += (de - ds);
        dcom += (de - ds);
    }
#endif

#pragma omp barrier

#if TIMER_PRINT
    ds = eigen_get_wtime();
#endif

#pragma omp for
    for (m_0 = 1; m_0 <= m; ++m_0) {
        i_1 = m_0 + (m_0 - 1) * nsm;
        if (ss[i_1 - 1] == ZERO) {
            ss[i_1 - 1] = 1.0;
        } else {
            ss[i_1 - 1] *= 0.5;
        }
    }

#pragma omp barrier
#pragma omp master
#if AT_BCAST_OVERLAP
#ifdef _OPENMP
    if (trbk_switched < 4) {
        trbk_decide_overlap_level(i);
    }
#endif
#endif
#if TIMER_PRINT
    de = eigen_get_wtime();
    dz += (de - ds);
#endif

    if (local_size == 1 || local_rank == 1) {
        ds = eigen_get_wtime();
    }
#pragma omp barrier

    eigen_trbakwy_block_body2<T>(z, nmz, v, nm, m,
                               &ss[ns0], &ss[0],
                               i_2, i_3, j_2, j_3);

    if (local_size == 1 || local_rank == 1) {
        de = eigen_get_wtime();
        trbk_time_trbk2 = de - ds;
#if TIMER_PRINT
        dy += (de - ds);
#endif
    }

#if AT_BCAST_OVERLAP
#ifdef _OPENMP
    if (trbk_switched < 4 && do_overlap_bcast_level == 0) {
#endif
#pragma omp barrier
    }
#endif
}

template <typename T>
void trbakwy_datacast(
    int64_t iloop_end, int64_t m, int64_t i,
    const T* a, int64_t nma,
    T* v, int64_t nm,
    T* ss,
#if !USE_BCASTW
    T* wk,
#endif
    int64_t iblk)
{
    int64_t j, iy, k0;
    int64_t i_1, j_1, j_4, j_5;
    int64_t jloop_sta, jloop_end;
    int64_t nodes[nsm];

#ifdef _OPENMP
    double ds = omp_get_wtime();
#endif

    std::fill(ss, ss + (iloop_end * m + ns0), static_cast<T>(ZERO));

#if USE_BCASTW

    jloop_end = std::min(eigen_loop_end(i + m - 1 - iblk, "X"), nm);
    if (m % y_nnod == 0) {
        for (j = 0; j < m; ++j) {
            if (y_inod == eigen_owner_node(i + j, "Y")) {
                iy  = eigen_owner_node(i + j, "Y");
                i_1 = eigen_translate_g2l(i + j, "Y", iy);
                k0  = (j / y_nnod) * jloop_end;
                std::copy(&a[(i_1 - 1) * nma], &a[(i_1 - 1) * nma + jloop_end], &v[k0]);
            }
        }

        k0 = (m / y_nnod) * jloop_end;
        allgather<T>(v, static_cast<T*>(trbk_buf.get()), k0, 1, y_COMM_WORLD);

        for (j = 0; j < m; ++j) {
            iy = eigen_owner_node(i + j, "Y");
            k0 = ((iy - 1) * (m / y_nnod) + (j / y_nnod)) * jloop_end;
            std::copy(&trbk_buf[k0], &trbk_buf[k0 + jloop_end], &v[j * nm]);

            jloop_sta = eigen_loop_start(i + j, "X");
            if (jloop_sta <= jloop_end) {
                std::fill(&v[j * nm + jloop_sta - 1], &v[j * nm + jloop_end], static_cast<T>(ZERO));
            }
        }

    } else {
        for (j = 0; j < m; ++j) {
            if (y_inod == eigen_owner_node(i + j, "Y")) {
                iy  = eigen_owner_node(i + j, "Y");
                i_1 = eigen_translate_g2l(i + j, "Y", iy);
                std::copy(&a[(i_1 - 1)* nma], &a[(i_1 - 1) * nma + jloop_end], &v[j * nm]);

                jloop_sta = eigen_loop_start(i + j, "X");
                if (jloop_sta <= jloop_end) {
                    std::fill(&v[j * nm + jloop_sta - 1], &v[j * nm + jloop_end], static_cast<T>(ZERO));
                }
            }
        }

        for (j = 0; j < m; j += y_nnod) {
            iy = eigen_owner_node(i + j, "Y");
            k0 = std::min(m - j, y_nnod);
            bcastw<T>(&v[j * nm], jloop_end, iy, nm, k0,static_cast<T*>(trbk_buf.get()), 2, y_COMM_WORLD);
        }
    }

#else

    if (m > y_nnod && y_nnod > 1) {
        for (j = 0; j < m; ++j) {
            nodes[j] = eigen_owner_node(i + j, "Y");
        }

        for (iy = 1; iy <= y_nnod; ++iy) {
            jloop_sta = eigen_loop_start(1, "X");
            jloop_end = eigen_loop_end(i + m - 1 - iblk, "X");

            k0 = 0;
            for (j = 0; j < m; ++j) {
                if (nodes[j] == iy) {
                    i_1 = eigen_translate_g2l(i + j, "Y", iy);
                    if (y_inod == iy) {
                        for (j_1 = jloop_sta; j_1 <= jloop_end; ++j_1) {
                            wk[k0 + j_1 - jloop_sta] = a[i_1 * nma + j_1 - 1];
                        }
                    }
                    k0 += (jloop_end - jloop_sta + 1);
                }
            }

            bcast<T>(wk, k0, iy, 2, y_COMM_WORLD);

            k0 = 0;
            for (j = 0; j < m; ++j) {
                if (nodes[j] == iy) {
                    for (j_1 = jloop_sta; j_1 <= jloop_end; ++j_1) {
                        v[j * nm + j_1 - 1] = wk[k0 + j_1 - jloop_sta];
                    }
                    k0 += (jloop_end - jloop_sta + 1);
                }
            }

            for (j = 0; j < m; ++j) {
                if (nodes[j] == iy) {
                    j_4 = eigen_loop_start(i + j, "X");
                    j_5 = eigen_loop_end(i + m - 1 - iblk, "X");
                    for (j_1 = j_4; j_1 <= j_5; ++j_1) {
                        v[j * nm + j_1 - 1] = static_cast<T>(ZERO);
                    }
                }
            }
        }

    } else {
        for (j = 0; j < m; ++j) {
            nodes[j] = eigen_owner_node(i + j, "Y");
            if (nodes[j] == y_inod) {
                i_1 = eigen_translate_g2l(i + j, "Y");
                jloop_sta = eigen_loop_start(1, "X");
                jloop_end = eigen_loop_end(i + m - 1 - iblk, "X");
                for (j_1 = jloop_sta; j_1 <= jloop_end; ++j_1) {
                    v[j * nm + j_1 - 1] = a[i_1 * nma + j_1 - 1];
                }

                jloop_sta = eigen_loop_start(i + j, "X");
                jloop_end = eigen_loop_end(i + m - 1 - iblk, "X");
                for (j_1 = jloop_sta; j_1 <= jloop_end; ++j_1) {
                    v[j * nm + j_1 - 1] = static_cast<T>(ZERO);
                }
            }
        }

        jloop_sta = eigen_loop_start(1, "X");
        jloop_end = eigen_loop_end(i + m - 1 - iblk, "X");
        for (j = 0; j < m; ++j) {
            bcast<T>(&v[j * nm + jloop_sta - 1], jloop_end - jloop_sta + 1,
                        nodes[j], 2, y_COMM_WORLD);
        }
    }

#endif

#ifdef _OPENMP
     double de = omp_get_wtime();
     trbk_time_bcast = de - ds;
#endif
}

template <typename T>
void eigen_trbakwy_body(
    int64_t n, int64_t nvec,
    const T* a, int64_t nma,
    T* z, int64_t nmz,
    T* beta,
    T* v1, T* v2, T* v3, int64_t nm,
    int64_t m,
    T* ss1, T* ss2, T* ss3,
    T* tt,
    int64_t iblk,
    int64_t nss
) {

    std::vector<T> wk;
    int64_t nodes[nsm];
    int64_t nx, ierr = 0;
    int64_t lwk;
    int64_t i, L, mode;
    int64_t x_root;
    int64_t iloop_sta, iloop_end;
    int64_t jloop_sta, jloop_end;
    int64_t i_1, i_4, j_1;
    T s0, s1, s2, s3;
    T d0, d1, d2, ds, de, dcom, dx, dy, dz;

#if TIMER_PRINT
    dx = ZERO;
    dy = ZERO;
    dz = ZERO;
#endif

    #pragma omp master
    {
        MPI_Barrier(TRD_COMM_WORLD);
    }
    #pragma omp barrier

    #pragma omp master
    {
        d1 = 0.0;
#if TIMER_PRINT
        d1 = eigen_get_wtime();
#endif
        dcom = ZERO;
    }
    #pragma omp barrier

#if AT_BCAST_OVERLAP
    #pragma omp master
    {
        do_overlap_bcast_level  = 2;
        trbk_time_counter  = 0;
        trbk_time_interval = 1;
        trbk_time_next     = 2;
        trbk_switched      = 0;
        trbk_time_trbk1    = ZERO;
        trbk_time_trbk1_   = ZERO;
        trbk_time_trbk1x   = ZERO;
        trbk_time_trbk1x_  = ZERO;
        trbk_time_trbk2    = ZERO;
        trbk_time_bcast    = ZERO;
        trbk_time_reduc    = ONE;
        trbk_time_fr       = ZERO;

        sync_other_than_master_init(trbk_lock, trbk_mask);
    }
#else
    #pragma omp master
    {
        do_overlap_bcast_level = 0;
    }
#endif
    #pragma omp barrier

    nx = std::min(( (n - (1 + iblk) + 1) % m ) + (1 + iblk) - 1, n);

    #pragma omp master
    {
        lwk = ((m - 1) / y_nnod + 1) * ((n - 1) / x_nnod + 1);
        lwk = std::max<int64_t>(lwk, n);
        wk.resize(lwk);
    }
    #pragma omp barrier

    #pragma omp for
    for (i = 1; i <= nx; i++) {
        int64_t l = i - iblk;

        if (i >= 1 + iblk &&
            ((i - 1) % y_nnod) + 1 == y_inod &&
            ((l - 1) % x_nnod) + 1 == x_inod) {

            i_1 = (i - y_inod) / y_nnod + 1;
            j_1 = (l - x_inod) / x_nnod + 1;

            beta[i - 1] = a[(j_1 - 1) + (i_1 - 1) * nma] * beta[i - 1];

        } else {
            beta[i - 1] = ZERO;
        }
    }

    #pragma omp master
    {
        reduce<T>(beta, wk.data(), nx, 1, TRD_COMM_WORLD);
    }
    #pragma omp barrier

    #pragma omp for
    for (i = 1 + iblk; i <= nx; i++) {
        if (beta[i - 1] == static_cast<T>(ZERO)) {
            s0 = static_cast<T>(ONE);
        } else {
            s0 = static_cast<T>(ONE) / beta[i - 1];
        }
        beta[i - 1] = s0;
    }

    iloop_sta = eigen_loop_start(1, "Y");
    iloop_end = eigen_loop_end(nvec, "Y");

    nx = std::min(((n - (1 + iblk) + 1) % m) + (1 + iblk) - 1, n);

    for (i = 1 + iblk; i <= nx; i++) {

        #pragma omp barrier

        if (beta[i - 1] == ZERO)
            continue;

        jloop_sta = eigen_loop_start(1, "X");
        jloop_end = eigen_loop_end(i - iblk, "X");

        i_4 = ( (iloop_end - iloop_sta + 1) % 4 ) + iloop_sta;

        #pragma omp master
        {
            ds = eigen_get_wtime();
            nodes[0] = eigen_owner_node(i, "Y");
            if (nodes[0] == y_inod) {
                i_1 = eigen_translate_g2l(i, "Y");
                for (j_1 = jloop_sta; j_1 <= jloop_end; j_1++) {
                    v1[(j_1 - 1)] = a[(j_1 - 1) + (i_1 - 1) * nma];
                }
            }

            bcast<T>(
                &v1[(jloop_sta - 1)],
                jloop_end - jloop_sta + 1,
                nodes[0], 1, y_COMM_WORLD
            );

            de = eigen_get_wtime();
            dcom += (de - ds);
        }
        #pragma omp barrier

#ifndef USE_BLAS
#define USE_BLAS 1
#endif

        #pragma omp master
        {
            barrier(y_COMM_WORLD);
            barrier(x_COMM_WORLD);
        }

        if (jloop_end >= jloop_sta) {

            i_4 = ((iloop_end - iloop_sta + 1) % 8) + iloop_sta;

            if (i_4 > iloop_sta) {
                #pragma omp for schedule(static)
                for (i_1 = iloop_sta; i_1 <= i_4 - 1; i_1++) {

#if USE_BLAS
                    ss1[i_1 - 1] =
                        dot<T>(jloop_end - jloop_sta + 1,
                             &z[(jloop_sta - 1) + (i_1 - 1) * nmz], 1,
                             &v1[(jloop_sta - 1)], 1);
#else
                    s0 = ZERO;
                    for (j_1 = jloop_sta; j_1 <= jloop_end; j_1++) {
                        s0 += v1[(j_1 - 1)] *
                              z[(j_1 - 1) + (i_1 - 1) * nmz];
                    }
                    ss1[i_1 - 1] = s0;
#endif
                }
            }

            #pragma omp for schedule(static)
#if USE_BLAS
            for (i_1 = i_4; i_1 <= iloop_end; i_1 += 8) {
                gemv<T>(
                    CblasColMajor,
                    CblasTrans,
                    jloop_end - jloop_sta + 1, 8,
                    static_cast<T>(ONE),
                    &z[(jloop_sta - 1) + (i_1 - 1) * nmz], nmz,
                    &v1[(jloop_sta - 1)], 1,
                    static_cast<T>(ZERO),
                    &ss1[i_1 - 1], 1
                );
            }
#else
            for (i_1 = i_4; i_1 <= iloop_end; i_1 += 4) {

                s0 = ZERO;
                s1 = ZERO;
                s2 = ZERO;
                s3 = ZERO;

                for (j_1 = jloop_sta; j_1 <= jloop_end; j_1++) {
                    T v = v1[(j_1 - 1)];
                    s0 += v * z[(j_1 - 1) + (i_1 + 0 - 1) * nmz];
                    s1 += v * z[(j_1 - 1) + (i_1 + 1 - 1) * nmz];
                    s2 += v * z[(j_1 - 1) + (i_1 + 2 - 1) * nmz];
                    s3 += v * z[(j_1 - 1) + (i_1 + 3 - 1) * nmz];
                }

                ss1[i_1 + 0 - 1] = s0;
                ss1[i_1 + 1 - 1] = s1;
                ss1[i_1 + 2 - 1] = s2;
                ss1[i_1 + 3 - 1] = s3;
            }
#endif

        } else {

            #pragma omp for
            for (i_1 = iloop_sta; i_1 <= iloop_end; i_1++) {
                ss1[i_1 - 1] = static_cast<T>(ZERO);
            }
        }

        #pragma omp barrier

        #pragma omp master
        {
            ds = eigen_get_wtime();

            reduce<T>(
                &ss1[iloop_sta - 1], tt, iloop_end - iloop_sta + 1, 2,
                x_COMM_WORLD
            );

            de = eigen_get_wtime();
            dcom += (de - ds);

#if USE_BLAS
            scal<T>(iloop_end - iloop_sta + 1, beta[i - 1], &ss1[iloop_sta - 1], 1);
#else
            s0 = beta[i - 1];
            for (i_1 = iloop_sta; i_1 <= iloop_end; i_1++) {
                ss1[i_1 - 1] *= s0;
            }
#endif
        }

        #pragma omp barrier
        jloop_sta = eigen_loop_start(1, "X");
        jloop_end = eigen_loop_end(i - iblk, "X");
        i_4 = ((iloop_end - iloop_sta + 1) % 8) + iloop_sta;

        if (jloop_end >= jloop_sta) {
          if (i_4 > iloop_sta) {
             #pragma omp for schedule(static)
             for (i_1 = iloop_sta; i_1 <= i_4 - 1; i_1++) {

#if USE_BLAS
                 axpy<T>(
                   jloop_end - jloop_sta + 1,
                   ss1[i_1 - 1], &v1[(jloop_sta - 1)], 1, &z[(jloop_sta - 1) + (i_1 - 1) * nmz], 1
                 );
#else
                 s0 = ss1[i_1 - 1];
                 for (j_1 = jloop_sta; j_1 <= jloop_end; j_1++) {
                     z[(j_1 - 1) + (i_1 - 1) * nmz] += s0 * v1[(j_1 - 1)];
                 }
#endif
             }
          }
          #pragma omp for schedule(static)
#if USE_BLAS
          for (i_1 = i_4; i_1 <= iloop_end; i_1 += 8) {
            ger<T>(
                jloop_end - jloop_sta + 1, 8,
                ONE, &v1[(jloop_sta - 1)], 1,
                &ss1[i_1 - 1], 1, &z[(jloop_sta - 1) + (i_1 - 1) * nmz], nmz
            );
          }
#else
          for (i_1 = i_4; i_1 <= iloop_end; i_1 += 4) {

              s0 = ss1[i_1 + 0 - 1];
              s1 = ss1[i_1 + 1 - 1];
              s2 = ss1[i_1 + 2 - 1];
              s3 = ss1[i_1 + 3 - 1];

              for (j_1 = jloop_sta; j_1 <= jloop_end; j_1++) {

                  T v = v1[(j_1 - 1)];

                  z[(j_1 - 1) + (i_1 + 0 - 1) * nmz] += s0 * v;
                  z[(j_1 - 1) + (i_1 + 1 - 1) * nmz] += s1 * v;
                  z[(j_1 - 1) + (i_1 + 2 - 1) * nmz] += s2 * v;
                  z[(j_1 - 1) + (i_1 + 3 - 1) * nmz] += s3 * v;
              }
          }
#endif
        }
    }


    #pragma omp barrier

    d2 = eigen_get_wtime();
    d0 = d2 - d1;

    mode = 0;
    i = nx + 1;

    #pragma omp master
    {
        if (nx + 1 <= n) {
            ds = eigen_get_wtime();
            trbakwy_datacast<T>(
                iloop_end, m, i, a, nma, 
                v1, nm, ss1, 
#if !USE_BCASTW
                wk.data(),
#endif
                iblk
            );

            de = eigen_get_wtime();
            dcom += (de - ds);
        }
    }
    #pragma omp barrier

    #pragma omp master
    {
        if (nx + 1 + m <= n) {
            ds = eigen_get_wtime();
            trbakwy_datacast<T>(
                iloop_end, m, i + m, a, nma,
                v2, nm, ss2,
#if !USE_BCASTW
                wk.data(),
#endif
                iblk
            );

            de = eigen_get_wtime();
            dcom += (de - ds);
        }
    }

    for (i = nx + 1; i <= n; i += m) {
        if (mode == 0) {
            eigen_trbakwy_block_body<T>(iloop_end, z, nmz,
                v1, nm, m,
                i, ss1, tt, nss, iblk,
                dcom, dx, dy, dz
            );
        }

        if (mode == 1) {
            eigen_trbakwy_block_body<T>(iloop_end, z, nmz,
                v2, nm, m,
                i, ss2, tt, nss, iblk,
                dcom, dx, dy, dz
            );
        }

        if (mode == 2) {
            eigen_trbakwy_block_body<T>(iloop_end, z, nmz,
                v3, nm, m,
                i, ss3, tt, nss, iblk,
                dcom, dx, dy, dz
            );
        }

        #pragma omp master
        {
            if (i + 2 * m <= n) {
                ds = eigen_get_wtime();
                if (mode == 0) {
                    trbakwy_datacast<T>(
                        iloop_end, m, i + 2 * m, a, nma,
                        v3, nm, ss3,
#if !USE_BCASTW
                        wk.data(),
#endif
                        iblk
                    );
                }

                if (mode == 1) {
                    trbakwy_datacast<T>(
                        iloop_end, m, i + 2 * m, a, nma,
                        v1, nm, ss1,
#if !USE_BCASTW
                        wk.data(),
#endif
                        iblk
                    );
                }

                if (mode == 2) {
                    trbakwy_datacast<T>(
                        iloop_end, m, i + 2 * m, a, nma,
                        v2, nm, ss2,
#if !USE_BCASTW
                        wk.data(),
#endif
                        iblk
                    );
                }

                de = eigen_get_wtime();
                dcom += (de - ds);
            }
        }

        #pragma omp barrier

        mode = (mode + 1) % 3;
    }
    #pragma omp barrier

    if (omp_get_num_threads() > 1) {
        #pragma omp barrier
        if (omp_get_thread_num() == 1) {
            ss1[0] = dx;
            ss1[1] = dy;
        }

        #pragma omp barrier
        if (omp_get_thread_num() == 0) {
            dx = ss1[0];
            dy = ss1[1];
        }
    }

    #pragma omp barrier

#if AT_BCAST_OVERLAP
    #pragma omp master
    {
        sync_other_than_master_finalize(trbk_lock);
    }
#endif

    #pragma omp master
    {
        MPI_Barrier(TRD_COMM_WORLD);
    }
    #pragma omp barrier

#if TIMER_PRINT > 1
    #pragma omp master
    {
        d2 = eigen_get_wtime();
        if (TRD_inod == 1) {
            std::cout << "TRBAK= " << (d2 - d1) << "\n";
            std::cout << "COMM= "  << dcom << "\n";
            std::cout << "   " << (d2 - d1) << " "
                      << (2.0 * nvec * n * n) / (d2 - d1) * 1e-9
                      << " GFLOPS\n";

            if (dx > 0) {
                std::cout << "   " << dx << " "
                          << (1.0 * nvec * n * n) / dx * 1e-9
                          << " GFLOPS\n";
            }

            if (dy > 0) {
                std::cout << "   " << dy << " "
                          << (1.0 * nvec * n * n) / dy * 1e-9
                          << " GFLOPS\n";
            }

            std::cout << "   " << d0 << " " << dz << "\n";
        }
    }
#endif

    return;
}

template <typename T>
void eigen_common_trbakwy(
    int64_t n,
    int64_t nvec,
    const T* a, int64_t nma0,
    T* z, int64_t nmz0,
    T* beta,
    int64_t m0,
    int64_t iblk
) {
    if (nvec == 0) return;

    int64_t nma = nma0;
    int64_t nmz = nmz0;

    int64_t m = std::min(nsm, m0);
    if (m < 1) m = 1;

    eigen_timer_reset(2, 4, 0, 1);

    int64_t na = (n - 1) / y_nnod + 1;
    na = na + ((na - 1) % 2);

    int64_t nm;
    get_optdim(nma, 9, 16 * 4, 16 * 6, nm);

    int64_t len_v  = std::max(nm * m, n);
    int64_t len_ss = na * m + ns0;

#if BOOST_BY_CACHE_ALIGNMENT
    len_v  += n_columns;
    len_ss += n_columns;
#endif

    auto v1  = std::make_unique<T[]>(len_v);
    auto ss1 = std::make_unique<T[]>(len_ss);
    auto v2  = std::make_unique<T[]>(len_v);
    auto ss2 = std::make_unique<T[]>(len_ss);
    auto v3  = std::make_unique<T[]>(len_v);
    auto ss3 = std::make_unique<T[]>(len_ss);
    auto tt  = std::make_unique<T[]>(len_ss);

    std::fill(v1.get(),  v1.get()  + len_v,  T(0));
    std::fill(ss1.get(), ss1.get() + len_ss, T(0));
    std::fill(v2.get(),  v2.get()  + len_v,  T(0));
    std::fill(ss2.get(), ss2.get() + len_ss, T(0));
    std::fill(v3.get(),  v3.get()  + len_v,  T(0));
    std::fill(ss3.get(), ss3.get() + len_ss, T(0));
    std::fill(tt.get(),  tt.get()  + len_ss, T(0));

    int64_t i_v1 = 0, i_s1 = 0;
    int64_t i_v2 = 0, i_s2 = 0;
    int64_t i_v3 = 0, i_s3 = 0;
    int64_t i_t  = 0;

#if BOOST_BY_CACHE_ALIGNMENT
    adjust_base(v1.get(),  z, i_v1);
    adjust_base(ss1.get(), z, i_s1);
    adjust_base(v2.get(),  z, i_v2);
    adjust_base(ss2.get(), z, i_s2);
    adjust_base(v3.get(),  z, i_v3);
    adjust_base(ss3.get(), z, i_s3);
    adjust_base(tt.get(),  z, i_t);

    int64_t kx = (L1_WINDOW / 8)
           + (L1_LSIZE)
           + (L2_LSIZE / 8);

    i_v1 += kx * 5;
    i_s1 += kx * 1;
    i_v2 += kx * 5;
    i_s2 += kx * 1;
    i_v3 += kx * 5;
    i_s3 += kx * 1;
    i_t  += kx * 1;

    round_offset(&i_v1);
    round_offset(&i_s1);
    round_offset(&i_v2);
    round_offset(&i_s2);
    round_offset(&i_v3);
    round_offset(&i_s3);
    round_offset(&i_t);
#endif

    trbk_buf = std::make_unique<double[]>(m * nm);
    trbk_time_reduc_overhead_x = reduce_cont_overhead_x;
    MPI_Barrier(TRD_COMM_WORLD);
#pragma omp parallel
    {
        eigen_trbakwy_body<T>(
            n, nvec,
            a, nma,
            z, nmz,
            beta,
            v1.get() + i_v1, v2.get() + i_v2, v3.get() + i_v3, nm, m,
            ss1.get() + i_s1, ss2.get() + i_s2, ss3.get() + i_s3,
            tt.get()  + i_t, iblk, na
        );
    }

    MPI_Barrier(TRD_COMM_WORLD);
    double comm_time_backtrafo =
        eigen_timer_print("EigenExa (Back-transformation)");
}

}
#endif
