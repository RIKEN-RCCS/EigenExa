#pragma once
#ifndef FS_PDLASRT_HPP
#define FS_PDLASRT_HPP

#include <mpi.h>

#include <algorithm>
#include <cstdio>
#include <numeric>

#include "../blas.hpp"
#include "FS_dividing.hpp"
#include "FS_libs.hpp"
#include "FS_prof.hpp"
namespace FS_pdlasrt {
using std::printf;
template <class Float>
void FS_pdlasrt(int n, Float d[], Float q[], int ldq,
                const FS_dividing::bt_node<Float> &subtree, Float q2[],
                int ldq2, Float sendq[], Float recvq[], Float buf[],
                int indrow[], int indcol[], int indx[], int indrcv[],
                FS_prof::FS_prof &prof) {
#ifdef _DEBUGLOG
  if (FS_libs::get_myrank() == 0) {
    printf("FS_pdlasrt start.");
  }
#endif
#if TIMER_PRINT
  prof.start(70);
#endif
  //

  const auto grid_info = subtree.FS_grid_info();
  const auto nprow = grid_info.nprow;
  const auto npcol = grid_info.npcol;
  const auto myrow = grid_info.myrow;
  const auto mycol = grid_info.mycol;
  const auto nblk = subtree.FS_get_NBLK();
  const auto nb = subtree.FS_get_NB();
  const auto np = (nblk / nprow) * nb;
  const auto nq = (nblk / npcol) * nb;
//
#pragma omp parallel for schedule(static, 1)
  for (auto i = 0; i < n; i += nb) {
    const auto row = subtree.FS_info_G1L('R', i).rocsrc;
    const auto col = subtree.FS_info_G1L('C', i).rocsrc;
    for (auto j = 0; j < nb; j++) {
      if (i + j < n) {
        indrow[i + j] = row;
        indcol[i + j] = col;
      }
    }
  }
  //
  // Sort the eigenvalues in D
  //
  std::iota(indx, indx + n, 0);
  std::sort(indx, indx + n, [&d](int i1, int i2) { return d[i1] < d[i2]; });

#pragma omp parallel
  {
#pragma omp for
    for (auto i = 0; i < n; i++) {
      buf[i] = d[indx[i]];
    }
#pragma omp for
    for (auto i = 0; i < n; i++) {
      d[i] = buf[i];
    }
  }

  // 列方向の入れ替え
  // 固有値昇順に合わせたソートと非ブロックサイクリック化
  for (auto pj = 0; pj < npcol; pj++) {
    // 入れ替え対象のプロセス列
    // デッドロックしないように逆順に回す
    const auto pjcol = (npcol - 1 - mycol - pj + npcol) % npcol;

    //
    int nsend = 0;
    int nrecv = 0;
    for (int j = 0; j < n; j++) {
      // 元のグローバルインデクス
      const auto gi = indx[j];

      // 元のグローバルインデクスを保持するプロセス列
      const auto col = indcol[gi];

      // 入れ替え先インデクスJを保持するプロセス列
      const auto jcol = j % npcol;
      // PJCOLに送信する列
      if (col == mycol && jcol == pjcol) {
        // ローカルインデクス
        const auto jl = subtree.FS_index_G2L('C', gi);

        // 送信バッファに格納
        std::copy_n(&q[jl * ldq], np, &sendq[nsend * np]);
        nsend += 1;
      }

      // PJCOLから受信する列
      if (col == pjcol && jcol == mycol) {
        // 格納先のローカルインデクス
        indrcv[nrecv] = j / npcol;

        // 受信数
        nrecv += 1;
      }
    }

    // irecv
    MPI_Request req;
    if (nrecv > 0) {
      MPI_Irecv(recvq, np * nrecv, FS_const::MPI_TYPE<Float>,
                subtree.group_Y_processranklist_[pjcol], 1,
                FS_libs::FS_COMM_WORLD, &req);
    }

    // send
    if (nsend > 0) {
      MPI_Send(sendq, np * nsend, FS_const::MPI_TYPE<Float>,
               subtree.group_Y_processranklist_[pjcol], 1,
               FS_libs::FS_COMM_WORLD);
    }

    // waitと展開
    MPI_Status stat;
    if (nrecv > 0) {
      MPI_Wait(&req, &stat);
#pragma omp parallel for
      for (auto j = 0; j < nrecv; j++) {
        const auto jl = indrcv[j];
        std::copy_n(&recvq[j * np], np, &q2[jl * ldq2]);
      }
    }
  }
  //======================================================================
  // 列方向の入れ替え
  // 非ブロックサイクリック化
  for (auto pj = 0; pj < nprow; pj++) {
    // 入れ替え対象のプロセス行
    // デッドロックしないように逆順に回す
    const auto pjrow = (nprow - 1 - myrow - pj + nprow) % nprow;

    //
    int nsend = 0;
    int nrecv = 0;
    for (int i = 0; i < n; i++) {
      // 元のグローバルインデクスを保持するプロセス行
      const auto row = indrow[i];
      // 入れ替え先インデクスIを保持するプロセス行
      const auto irow = i % nprow;

      // PJROWに送信する行
      if (row == myrow && irow == pjrow) {
        // ローカルインデクス
        const auto il = subtree.FS_index_G2L('R', i);

        // 送信バッファに格納
        lapacke::copy<Float>(nq, &q2[il], ldq2, &sendq[nsend * nq], 1);
        nsend += 1;
      }

      // PJROWから受信する例
      if (row == pjrow && irow == myrow) {
        // 格納先のローカルインデクス
        indrcv[nrecv] = i / nprow;

        // 受信数
        nrecv += 1;
      }
    }
    MPI_Request req;
    // irecv
    if (nrecv > 0) {
      MPI_Irecv(recvq, nrecv * nq, FS_const::MPI_TYPE<Float>,
                subtree.group_X_processranklist_[pjrow], 1,
                FS_libs::FS_COMM_WORLD, &req);
    }

    // send
    if (nsend > 0) {
      MPI_Send(sendq, nsend * nq, FS_const::MPI_TYPE<Float>,
               subtree.group_X_processranklist_[pjrow], 1,
               FS_libs::FS_COMM_WORLD);
    }

    // waitと展開
    if (nrecv > 0) {
      MPI_Status stat;
      MPI_Wait(&req, &stat);

#pragma omp parallel for
      for (int i = 0; i < nrecv; i++) {
        const auto il = indrcv[i];
        lapacke::copy<Float>(nq, &recvq[i * nq], 1, &q[il], ldq);
      }
    }
  }

#if TIMER_PRINT
  prof.end(70);
#endif

#ifdef _DEBUGLOG
  if (FS_libs::get_myrank() == 0) {
    printf("FS_pdlasrt end.");
  }
#endif
}
}  // namespace FS_pdlasrt

#endif
