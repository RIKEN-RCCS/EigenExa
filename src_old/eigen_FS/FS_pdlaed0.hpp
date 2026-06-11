#pragma once
#ifndef FS_PDLAED0_HPP
#define FS_PDLAED0_HPP

#include <mpi.h>

#include <cstdio>

#include "../blas.hpp"
#include "../eigen/eigen_libs0.hpp"
#include "FS2eigen_pdlasrt.hpp"
#include "FS_const.hpp"
#include "FS_dividing.hpp"
#include "FS_libs.hpp"
#include "FS_pdlaed1.hpp"
#include "FS_pdlasrt.hpp"
#include "FS_prof.hpp"

namespace eigen_FS {

template <typename Integer, typename Float>
Integer FS_pdlaed0(Integer n, Float d[], Float e[], Float q[], Integer ldq,
                   Float work[], Integer lwork, Integer iwork[], Integer liwork,
                   FS_prof &prof) {
#ifdef _DEBUGLOG
  if (FS_libs::get_myrank() == 0) {
    printf("FS_PDLAED0 start.\n");
  }
#endif
  FS_dividing::bt_node<Float> root_node = {};
  const auto info = [&]() mutable -> int {
    auto nnod = (FS_libs::is_comm_member()) ? FS_libs::get_procs()
                                            : FS_libs::Nod{1, 1, 1};

    if (FS_libs::is_comm_member()) {
#if TIMER_PRINT
      prof.start(20);
#endif
      std::unique_ptr<bool[]> hint(new bool[nnod.nod]);
      FS_dividing::FS_create_hint(hint.get());
      const auto info_dividing =
          root_node.FS_dividing(n, d, e, std::move(hint), prof);
      if (info_dividing != 0) {
        return info_dividing;
      }
      const auto NBLK = root_node.FS_get_NBLK();
      const auto NB = root_node.FS_get_NB();
      const Integer NP = (NBLK / root_node.x_nnod_) * NB;
      const Integer NQ = (NBLK / root_node.y_nnod_) * NB;

#pragma omp parallel for collapse(2)
      for (Integer j = 0; j < NQ; j++) {
        for (Integer i = 0; i < NP; i++) {
          q[i + j * ldq] = FS_const::ZERO<Float>;
        }
      }

      const auto getleaf = root_node.FS_dividing_getleaf(0);
      const auto leafnode = getleaf.first;
      if (getleaf.second != 0) {
        return getleaf.second;
      }

#if TIMER_PRINT
      prof.start(28);
#endif

      const auto id = leafnode->nstart_;
      const auto mat_size = leafnode->FS_get_N_active();
      if (mat_size > 0) {
        const auto pq = root_node.FS_info_G2L(id, id);
        int info = 0;
        if (id < n - 1) {
          info = lapacke::stedc<Float>('I', mat_size, &d[id], &e[id],
                                       &q[pq.row + pq.col * ldq], ldq, work,
                                       lwork, iwork, liwork);
        } else {
          q[pq.row + pq.col * ldq] = 1.0;
          info = 0;
        }

        if (info != 0) {
          return info;
        }
      }

#if TIMER_PRINT
      prof.end(28);
#endif
      for (auto node = leafnode; node->parent_node_ != nullptr; node = node->parent_node_) {
        const auto parent_node = node->parent_node_;
        const auto id = parent_node->nstart_;
        const auto n0 = parent_node->FS_get_N_active();
        const auto n1 = parent_node->sub_bt_node_[0].FS_get_N_active();

        if (n0 == n1) {
          continue;
        }

        const auto rho = e[id + n1 - 1];

        const auto q_top = parent_node->FS_get_QTOP();

#ifdef _DEBUGLOG
        const auto nb = parent_node->FS_get_NB();
        if (FS_libs::get_myrank() == 0) {
          printf("+---------------------\n");
          printf("FS_PDLAED0 merge loop\n");
          printf(" layer = %d\n", parent_node->layer_);
          printf(" N     = %d\n", n0);
          printf(" NB    = %d\n", nb);
        }
#endif

        FS_prof prof_layer = {};
#if TIMER_PRINT
        prof_layer.init();
#endif
        int info_pdlaed1 = FS_pdlaed1::FS_pdlaed1<Float>(
            n0, n1, &d[id], &q[q_top.row + q_top.col * ldq], ldq, *parent_node,
            rho, work, iwork, prof_layer);
#if TIMER_PRINT > 2
        prof_layer.finalize();
#endif
#if TIMER_PRINT
        prof.add(prof_layer);
#endif
        if (info_pdlaed1 != 0) {
          return info_pdlaed1;
        }
        continue;
      }


#ifdef _DEBUGLOG
      if (FS_libs::get_myrank() == 0) {
        printf("+---------------------\n");
      }
#endif
    }

#ifdef _DEBUGLOG
    if (FS_libs::get_myrank() == 0) {
      printf("START Bcast\n");
    }
#endif

    MPI_Comm eigen_comm, x_eigen_comm, y_eigen_comm;
    eigen_libs0::eigen_get_comm(eigen_comm, x_eigen_comm, y_eigen_comm);

    // TODO MPI_INTがtemplateで書き換えるべき部分
    MPI_Bcast(&nnod.x, 1, MPI_INT, 0, eigen_comm);
    MPI_Bcast(&nnod.y, 1, MPI_INT, 0, eigen_comm);

    const auto NBLK = root_node.FS_get_NBLK();
    const auto NB = root_node.FS_get_NB();
    const Integer NP = (NBLK / nnod.x) * NB;
    const Integer NQ = (NBLK / nnod.y) * NB;

    const auto ldq2 = NP;
    const Integer ipq2 = 0;

    const auto i_send_q = ipq2 + ldq2 * NQ;
    const auto i_recv_q = i_send_q + ldq2 * NQ;
    const auto i_buffer = i_recv_q + ldq2 * NQ;

    const Integer index_row = 0;
    const auto index_col = index_row + n;
    const auto index = index_col + n;
    const auto index_recv = index + n;

    int eigen_np, x_eigen_np, y_eigen_np;
    eigen_libs0::eigen_get_procs(eigen_np, x_eigen_np, y_eigen_np);

    if (nnod.nod == eigen_np) {
      FS_pdlasrt::FS_pdlasrt(n, d, q, ldq, root_node, &work[ipq2], ldq2,
                             &work[i_send_q], &work[i_recv_q], &work[i_buffer],
                             &iwork[index_row], &iwork[index_col],
                             &iwork[index], &iwork[index_recv], prof);
    } else {
      eigen_FS::FS2eigen::FS2eigen_pdlasrt(
          n, d, ldq, q, root_node, (int *)work, work,
          (eigen_FS::FS2eigen::GpositionValue<Float>
               *)&work[std::max(0, (NP * NQ / 2))],
          iwork, prof);
    }
    return 0;
  }();

  if (FS_libs::is_comm_member()) {
    root_node.FS_dividing_free();
  }

#if TIMER_PRINT
  prof.end(20);
#endif
#ifdef _DEBUGLOG
  if (FS_libs::get_myrank() == 0) {
    printf("FS_PDLAED0 end. info = %d\n", info);
  }
#endif
  return info;
}

}  // namespace eigen_FS

#endif
