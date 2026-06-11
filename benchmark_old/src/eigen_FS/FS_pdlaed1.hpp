#pragma once
#ifndef FS_PDLAED1_HPP
#define FS_PDLAED1_HPP

#include <mpi.h>
#include <cstdio>

#include <algorithm>

#include "FS_const.hpp"
#include "FS_dividing.hpp"
#include "FS_libs.hpp"
#include "FS_merge_d.hpp"
#include "FS_pdlaed2.hpp"
#include "FS_pdlaed3.hpp"
#include "FS_pdlaedz.hpp"
#include "FS_prof.hpp"
#include "FS_reduce_zd.hpp"
namespace FS_pdlaed1 {
using eigen_FS::FS_reduce_zd;
using std::printf;
template<class Float>
int FS_pdlaed1(int n, int n1, Float d[], Float q[], int ldq,
               const FS_dividing::bt_node<Float> &subtree, Float rho, Float work[],
               int iwork[], FS_prof::FS_prof &prof) {
#ifdef _DEBUGLOG
  if (FS_libs::get_myrank() == 0) {
    printf("FS_pdlaed1 start.\n");
  }
#endif
#if TIMER_PRINT
  prof.start(30);
#endif

  // Quick return if possible
  if (n == 0) {
    return 0;
  }
  // The following values are  integer pointers which indicate
  // the portion of the workspace used by a particular array
  // in FS_PDLAED2 and FS_PDLAED3.
  const auto grid_info = subtree.FS_grid_info();

  const auto nblk = subtree.FS_get_NBLK();
  const auto nb = subtree.FS_get_NB();
  const auto np = (nblk / grid_info.nprow) * nb;
  const auto nq = (nblk / grid_info.npcol) * nb;
  const auto ldq2 = std::max(np, 1);
  const auto ldu = nq;
  const auto lsendq2 = ldq2 * nq;
  //
  const int iz = 0;
  const auto z = &work[iz];
  const auto idlamda = iz + n;
  const auto dlamda = &work[idlamda];
  const auto iw = idlamda + n;
  const auto w = &work[iw];
  const auto ipq2 = iw + n;
  const auto q2 = &work[ipq2];
  const auto ipu = ipq2 + ldq2 * nq;
  const auto u = &work[ipu];
  const auto ibuf = ipu + ldu * nq;
  const auto buf = &work[ibuf];
  const auto isendq2 = ibuf + 4 * n;
  const auto sendq2 = &work[isendq2];
  const auto irecvq2 = isendq2 + lsendq2;
  const auto recvq2 = &work[irecvq2];
  //
  const int ictot = 0;
  const auto ipsm = ictot + grid_info.npcol * 4;
  const auto indx = ipsm + grid_info.npcol * 4;
  const auto indxc = indx + n;
  const auto indxp = indxc + n;
  const auto indcol = indxp + n;
  const auto coltyp = indcol + n;
  const auto indrow = coltyp + n;
  const auto indxr = indrow + n;
  const auto indxcb = indxr + n;
  //
  // for FS_merge_d, FS_pdlaedz, FS_reduce_zd
  const auto izwork = iz + n * 2;
  const auto idwork = izwork + n;
  //
  // merge d
  //
#if TIMER_PRINT
  prof.start(31);
#endif
  eigen_FS::FS_merge_d<Float>(n, d, subtree, &work[idwork]);
#if TIMER_PRINT
  prof.end(31);
#endif
  //
  // Form the z-vector which consists of the last row of Q_1 and the
  // first row of Q_2.
  //
  FS_pdlaedz::FS_pdlaedz<Float>(n, n1, q, ldq, subtree, &work[izwork], prof);
  //
  // MPI_allreduce d and z
  //
  FS_reduce_zd<Float>(n, subtree, &work[izwork], &work[iz], d, prof);
  //
  // Deflate eigenvalues.
  //
  const auto k = FS_pdlaed2::FS_pdlaed2(
      n, n1, d, q, ldq, subtree, rho, z, w, dlamda, ldq2, q2,
      &iwork[indx], &iwork[ictot], buf, &iwork[coltyp], &iwork[indcol],
      &iwork[indxc], &iwork[indxp], &iwork[ipsm], prof);


  //
  // Solve Secular Equation.
  //
  int lctot = subtree.y_nnod_;
  int info = 0;
  if (k != 0) {
    info = FS_pdlaed3::FS_pdlaed3_C(k, n, n1, d, rho, dlamda, w, ldq, q, subtree, ldq2, q2,
                           ldu, u, &iwork[indx], lctot, &iwork[ictot],
                           sendq2, recvq2, z, buf, &iwork[indrow],
                           &iwork[indcol], &iwork[indxc], &iwork[indxr],
                           &iwork[indxcb], prof);
  }
#if TIMER_PRINT
  prof.end(30);
#endif
#ifdef _DEBUGLOG
  if (FS_libs::get_myrank() == 0) {
    printf("FS_pdlaed1 end. INFO=%d\n", info);
  }
#endif
  return info;
}
}  // namespace FS_pdlaed1

#endif
