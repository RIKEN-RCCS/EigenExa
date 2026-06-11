#pragma once
#ifndef FS_PDLAED2_HPP
#define FS_PDLAED2_HPP

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <type_traits>

#include "../blas.hpp"
#include "FS_const.hpp"
#include "FS_dividing.hpp"
#include "FS_libs.hpp"
#include "FS_prof.hpp"
namespace FS_pdlaed2 {
template <class Float>
int get_NPA(int n, int nb, const FS_dividing::bt_node<Float> &subtree,
            int myrow) {
  int npa = 0;
#pragma omp parallel for reduction(+ : npa)
  for (int i = 0; i < n; i += nb) {
    const auto row = subtree.FS_info_G1L('R', i).rocsrc;
    if (row == myrow) {
      for (auto j = 0; j < nb; j++) {
        if (i + j < n) {
          npa += 1;
        }
      }
    }
  }
  return npa;
}

inline void init_ctot(int n, const int coltyp[], const int indcol[], int lctot,
                      int ctot[], int npcol) {
#pragma omp parallel for
  for (auto j = 0; j < 4; j++) {
    std::fill_n(&ctot[j * lctot], npcol, 0);
  }

  for (auto j = 0; j < n; j++) {
    auto ct = coltyp[j];
    auto col = indcol[j];
    ctot[ct * lctot + col] += 1;
  }
}

/**
 * \brief PSM(*) = Position in SubMatrix (of types 0 through 3)
 */
inline void set_psm(int lctot, int psm[], const int ctot[], int npcol) {
#pragma omp parallel for
  for (auto col = 0; col < npcol; col++) {
    psm[0 * lctot + col] = 1;
    psm[1 * lctot + col] = 1 + ctot[0 * lctot + col];
    psm[2 * lctot + col] = psm[1 * lctot + col] + ctot[1 * lctot + col];
    psm[3 * lctot + col] = psm[2 * lctot + col] + ctot[2 * lctot + col];
  }
}

inline void set_ptt(int lctot, int ptt[4], const int ctot[], int npcol) {
  std::fill_n(ptt, 4, 0);
  ptt[0] = 1;
#pragma omp parallel for
  for (auto i = 1; i < 4; i++) {
    int ct = 0;
    for (auto j = 0; j < npcol; j++) {
      ct += ctot[(i - 1) * lctot + j];
    }
    ptt[i] = ct;
  }

  for (auto i = 1; i < 4; i++) {
    ptt[i] += ptt[i - 1];
  }
}

template <class Float>
void set_indxp(int n, int &k2, int nj, int pj, Float d[], Float c, Float s,
               int indxp[]) {
  const auto c_2 = static_cast<Float>(pow(c, 2));
  const auto s_2 = static_cast<Float>(pow(s, 2));
  const auto t = d[pj] * c_2 + d[nj] * s_2;
  d[nj] = d[pj] * s_2 + d[nj] * c_2;
  d[pj] = t;

  int i;
  for (i = 0; k2 + i < n; i++) {
    if (d[pj] < d[indxp[k2 + i]]) {
      indxp[k2 + i - 1] = indxp[k2 + i];
      indxp[k2 + i] = pj;
    } else {
      break;
    }
  }
  k2 -= 1;
  indxp[k2 + i] = pj;
}

template <class Float>
void pdlaed2_comm(int mycol, int ldq, Float q[], int npa,
                  const FS_dividing::bt_node<Float> &subtree, Float qbuf[],
                  int nj, int pj, Float c, Float s, int indcol[]) {
  const auto njj_info = subtree.FS_info_G1L('C', nj);
  const auto &njj = njj_info.l_index;
  const auto &njcol = njj_info.rocsrc;
  const auto pjj_info = subtree.FS_info_G1L('C', pj);
  const auto &pjj = pjj_info.l_index;
  const auto &pjcol = pjj_info.rocsrc;

  if (indcol[pj] == indcol[nj] && mycol == njcol) {
    lapacke::rot<Float>(npa, &q[pjj * ldq + 0], 1, &q[njj * ldq + 0], 1, c, s);
  } else if (mycol == pjcol) {
    MPI_Request req[2];
    MPI_Irecv(qbuf, npa, FS_const::MPI_TYPE<Float>,
              subtree.group_Y_processranklist_[njcol], 1,
              FS_libs::FS_COMM_WORLD, &req[1]);

    MPI_Isend(&q[pjj * ldq + 0], npa, FS_const::MPI_TYPE<Float>,
              subtree.group_Y_processranklist_[njcol], 1,
              FS_libs::FS_COMM_WORLD, &req[0]);
    MPI_Waitall(2, req, MPI_STATUS_IGNORE);
    lapacke::rot<Float>(npa, &q[pjj * ldq + 0], 1, qbuf, 1, c, s);
  } else if (mycol == njcol) {
    MPI_Request req[2];
    MPI_Irecv(qbuf, npa, FS_const::MPI_TYPE<Float>,
              subtree.group_Y_processranklist_[pjcol], 1,
              FS_libs::FS_COMM_WORLD, &req[1]);

    MPI_Isend(&q[njj * ldq + 0], npa, FS_const::MPI_TYPE<Float>,
              subtree.group_Y_processranklist_[pjcol], 1,
              FS_libs::FS_COMM_WORLD, &req[0]);
    MPI_Waitall(2, req, MPI_STATUS_IGNORE);
    lapacke::rot<Float>(npa, qbuf, 1, &q[njj * ldq + 0], 1, c, s);
  }
}

template <class Float>
int FS_pdlaed2(int n, int n1, Float d[], Float q[], int ldq,
               const FS_dividing::bt_node<Float> &subtree, Float &rho,
               Float z[], Float w[], Float dlamda[], int ldq2, Float q2[],
               int indx[], int ctot[], Float qbuf[], int coltyp[], int indcol[],
               int indxc[], int indxp[], int psm[], FS_prof::FS_prof &prof) {
#ifdef _DEBUGLOG
  if (FS_libs::get_myrank() == 0) {
    printf("FS_PDLAED2 start.\n");
  }
#endif
#if TIMER_PRINT
  prof.start(50);
#endif
  // Quick return if possible
  if (n == 0) {
    return 0;
  }

  const auto grid_info = subtree.FS_grid_info();
  const auto npcol = grid_info.npcol;
  const auto myrow = grid_info.myrow;
  const auto mycol = grid_info.mycol;
  
  const auto nb = subtree.FS_get_NB();

  const auto n2 = n - n1;
  std::fill_n(dlamda, n, 0);

  if (rho < FS_const::ZERO<Float>) {
    lapacke::scal<Float>(n2, FS_const::MONE<Float>, &z[n1], 1);
  }
  // Normalize z so that norm(z) = 1.  Since z is the concatenation of
  // two normalized vectors, norm2(z) = sqrt(2).

  const auto T = FS_const::ONE<Float> / std::sqrt(FS_const::TWO<Float>);
  lapacke::scal<Float>(n, T, z, 1);

  rho = std::abs(FS_const::TWO<Float> * rho);

  // Calculate the allowable deflation tolerance
  const auto imax = lapacke::iamax<int, Float>(n, z, 1);
  const auto jmax = lapacke::iamax<int, Float>(n, d, 1);
  const auto abs_zmax = std::abs(z[imax]);
  const auto abs_dmax = std::abs(d[jmax]);
  const auto eps = std::numeric_limits<Float>::epsilon() / 2;
  const auto tol = FS_const::EIGHT<Float> * eps * std::max(abs_dmax, abs_zmax);

  // If the rank-1 modifier is small enough, no more needs to be done
  // except to reorganize Q so that its columns correspond with the
  // elements in D.
  if (rho * abs_zmax <= tol) {
    return 0;
  }

  // If there are multiple eigenvalues then the problem deflates.  Here
  // the number of equal eigenvalues are found.  As each equal
  // eigenvalue is found, an elementary reflector is computed to rotate
  // the corresponding eigensubspace so that the corresponding
  // components of Z are zero in this new basis.
  std::iota(indx, indx + n, 0);
  std::sort(indx, indx + n, [&d](int i1, int i2) { return d[i1] < d[i2]; });
#pragma omp parallel
  {
#pragma omp for
    for (auto i = 0; i < n1; i++) {
      coltyp[i] = 0;
    }
#pragma omp for
    for (auto i = n1; i < n; i++) {
      coltyp[i] = 2;
    }
#pragma omp for
    for (auto i = 0; i < n; i += nb) {
      const auto col = subtree.FS_info_G1L('C', i).rocsrc;
      for (auto j = 0; j < nb; j++) {
        if (i + j < n) {
          indcol[i + j] = col;
        }
      }
    }
  }

  const auto npa = get_NPA(n, nb, subtree, myrow);
  int k = 0;
  int k2 = n;
  int pj = 0;
  bool set_pj = true;
  for (int j = 0; j < n; j++) {
    const auto nj = indx[j];
    if (rho * std::abs(z[nj]) <= tol) {
      // Deflate due to small z component.
      k2 -= 1;
      coltyp[nj] = 3;
      indxp[k2] = nj;
    } else if (set_pj) {
      set_pj = false;
      pj = nj;
    } else {
      // Check if eigenvalues are close enough to allow deflation.
      const auto s_pj = z[pj];
      const auto c_nj = z[nj];
      // Find sqrt(a**2+b**2) without overflow or
      // destructive underflow.
      const auto tau = lapacke::lapy2<Float>(c_nj, s_pj);

      const auto t = d[nj] - d[pj];
      const auto c = c_nj / tau;
      const auto s = -s_pj / tau;
      if (std::abs(t * c * s) <= tol) {
        // Deflation is possible.
        z[nj] = tau;
        z[pj] = FS_const::ZERO<Float>;
        if (coltyp[nj] != coltyp[pj]) {
          coltyp[nj] = 1;
        }
        coltyp[pj] = 3;
#pragma omp parallel
        {
#pragma omp master
          {
            pdlaed2_comm(mycol, ldq, q, npa, subtree, qbuf, nj, pj, c, s,
                         indcol);
          }
#pragma omp single nowait
          { set_indxp(n, k2, nj, pj, d, c, s, indxp); }
        }
        pj = nj;
      } else {
        dlamda[k] = d[pj];
        w[k] = z[pj];
        indxp[k] = pj;
        k += 1;
        pj = nj;
      }
    }
  }

  // Record the last eigenvalue.
  dlamda[k] = d[pj];
  w[k] = z[pj];
  indxp[k] = pj;
  k += 1;
  // Count up the total number of the various types of columns, then
  // form a permutation which positions the four column types into
  // four uniform groups (although one or more of these groups may be
  // empty).

  int ptt[4];
  const auto lctot = subtree.y_nnod_;
  init_ctot(n, coltyp, indcol, lctot, ctot, npcol);
  // PSM[*] =  Position in SubMatrix (of types 0 through 3)
  set_psm(lctot, psm, ctot, npcol);
  set_ptt(lctot, ptt, ctot, npcol);
  // Fill out the INDXC array so that the permutation which it induces
  // will place all type-1 columns first, all type-2 columns next,
  // then all type-3's, and finally all type-4's.
  for (auto j = 0; j < n; j++) {
    auto js = indxp[j];
    auto col = indcol[js];
    auto ct = coltyp[js];
    auto i = subtree.FS_index_L2G('C', psm[col + ct * lctot] - 1, col);
    indx[j] = i;
    indxc[ptt[ct] - 1] = i + 1;
    psm[col + ct * lctot] += 1;
    ptt[ct] += 1;
  }
#pragma omp parallel for
  for (auto j = 0; j < n; j++) {
    const auto js = indxp[j];
    const auto col = indcol[js];
    if (col == mycol) {
      const auto jjs = subtree.FS_index_G2L('C', js);
      const auto i = indx[j];
      const auto jjq2 = subtree.FS_index_G2L('C', i);
      std::copy_n(&q[jjs * ldq], npa, &q2[jjq2 * ldq2]);
    }
  }

  std::copy_n(d, n, z);
#pragma omp parallel for
  for (auto j = k; j < n; j++) {
    const auto js = indxp[j];
    const auto i = indx[j];
    d[i] = z[js];
  }

#if TIMER_PRINT
  prof.end(50);
#endif

#ifdef _DEBUGLOG
  if (FS_libs::get_myrank() == 0) {
    printf("FS_PDLAED2 end.\n");
  }
#endif

  return k;
}
}  // namespace FS_pdlaed2
#endif
