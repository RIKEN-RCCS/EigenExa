/**
 * @file FS_EDC.hpp
 * @brief template int FS_EDC
 */
#pragma once
#ifndef FS_EDC_HPP
#define FS_EDC_HPP

#include <cmath>
#include <cstdio>
#include <type_traits>

#include "../blas.hpp"
#include "../eigen/eigen_libs0.hpp"
#include "FS_const.hpp"
#include "FS_libs.hpp"
#include "FS_pdlaed0.hpp"
#include "FS_prof.hpp"

namespace FS_EDC {
using FS_const::ONE;
using FS_const::ZERO;
using std::printf;
/**
 * @brief @n
 *   Purpose @n
 *   ======= @n
 *   FS_EDC computes all eigenvalues and eigenvectors of a             @n
 *   symmetric tridiagonal matrix in parallel, using the divide and    @n
 *   conquer algorithm.
 *
 *  @return              INFO  @n
 *                       = 0: successful exit   @n
 *                       !=0: error exit
 *
 *
 *
 * @param[in]     n      (global input) int @n
 *                       The order of the tridiagonal matrix T.  N >= 0.
 *
 * @param[in,out] D      (global input/output) float or double array, dimension
 * (N) @n On entry, the diagonal elements of the tridiagonal matrix.  @n On
 * exit, if INFO = 0, the eigenvalues in descending order.
 *
 * @param[in,out] E      (global input/output) float or double array,
 * dimension(N-1) @n On entry, the subdiagonal elements of the tridiagonal
 * matrix. @n On exit, E has been destroyed.
 *
 * @param[in,out] Q      (local output) float or double array, local dimension
 * (LDQ, *)   @n Q contains the orthonormal eigenvectors of the symmetric
 * tridiagonal matrix.   @n On output, Q is distributed across the P processes
 * in non block cyclic format. @n
 *
 * @param[in]     ldq    (local input) int @n
 *                       leading dimension of array Q.
 *
 * @param         work   (local workspace/output) float or double array,
 * dimension (LWORK)
 *
 * @param[in]     lwork  (local input/output) int, the dimension of the array
 * WORK. @n LWORK = 1 + 6*N + 3*NP*NQ + NQ*NQ @n LWORK can be obtained from
 * subroutine FS_WorkSize.
 *
 * @param         iwork  (local workspace/output) int array, dimension (LIWORK)
 *
 * @param[in]     liwork (input) int                   @n
 *                       The dimension of the array IWORK. @n
 *                       LIWORK = 1 + 8*N + 8*NPCOL        @n
 *                       LIWORK can be obtained from subroutine FS_WorkSize.
 *
 * @param[out]    prof   (global output) type(FS_prof) @n
 *                       profiling information of each subroutines.
 *
 * @note This routine is modified from ScaLAPACK PDSTEDC.f
 */

template <class Integer, class Float>
Integer FS_EDC(const Integer n, Float *D, Float *E, Float *Q, const Integer ldq,
               Float *work, Integer lwork, Integer *iwork, const Integer liwork,
               FS_prof::FS_prof *prof) {
  int nnod, x_nnod, y_nnod;
  int inod, x_inod, y_inod;

  FS_prof::FS_prof prof_tmp = {};

  Integer info = 0;
#ifdef _DEBUGLOG
  printf("FS_EDC start\n");
#endif

#if TIMER_PRINT
  if (prof != nullptr) {
    prof_tmp = *prof;
  } else {
    prof_tmp.init();
  }
  prof_tmp.start(10);
#endif

  eigen_libs0::eigen_get_procs(nnod, x_nnod, y_nnod);
  eigen_libs0::eigen_get_id(inod, x_inod, y_inod);

  if (n == 0) {
    // 何もしない
  } else if (n == 1) {
    if (x_inod == 1 && y_inod == 1) {
      Q[0] = ONE<Float>;
    }
  } else if (x_nnod * y_nnod == 1) {
    // If P=NPROW*NPCOL=1, solve the problem with DSTEDC.
#if TIMER_PRINT
    prof_tmp.start(11);
#endif
    info =
        lapacke::stedc<Float>('I', n, D, E, Q, ldq, work, lwork, iwork, liwork);

#if TIMER_PRINT
    prof_tmp.end(11);
#endif
  } else {
    // Scale matrix to allowable range, if necessary.
    auto orgnrm = lapacke::lanst<Integer, Float>('M', n, D, E);
    if (std::isnan(orgnrm)) {
      orgnrm = ZERO<Float>;
    }
    if (orgnrm != ZERO<Float>) {
      info = lapacke::lascl<Float>('G', 0, 0, orgnrm, ONE<Float>, n, 1, D, n);
      if (n - 1 >= 1) {
        info = lapacke::lascl<Float>('G', 0, 0, orgnrm, ONE<Float>, n - 1, 1, E,
                                     n - 1);
      }
    }
    info = eigen_FS::FS_pdlaed0<Integer, Float>(n, D, E, Q, ldq, work, lwork,
                                                iwork, liwork, prof_tmp);
    // Scale back.
    if (info == 0 && orgnrm != ZERO<Float>) {
      info = lapacke::lascl<Float>('G', 0, 0, ONE<Float>, orgnrm, n, 1, D, n);
    }
  }

#ifdef _DEBUGLOG
  printf("FS_EDC end. INFO=%d\n", info);
#endif

#if TIMER_PRINT
  prof_tmp.end(10);
  if (prof != nullptr) {
    *prof = prof_tmp;
  } else {
    prof_tmp.finalize();
  }
#endif
  return info;
}
}  // namespace FS_EDC

#endif
