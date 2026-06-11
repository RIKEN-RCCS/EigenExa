#pragma once
#ifndef FS_PDLAED3_HPP
#define FS_PDLAED3_HPP

#include <mpi.h>
#include <omp.h>

#include <cmath>
#include <memory>


#include "../blas.hpp"
#include "../eigen/eigen_dc.hpp"
#include "FS_const.hpp"
#include "FS_prof.hpp"
#include "MPI_Allreduce_group.hpp"

namespace FS_pdlaed3 {
using eigen_FS::MPI_Group_Allreduce;
using FS_libs::FS_COMM_WORLD;
using std::abs;
using std::max;
using std::min;
using std::sqrt;

int get_pdc(int lctot, const int ctot[], int npcol, int mycol) {
  int pdc = 0;
  for (auto col = 0; col != mycol; col = (col + 1) % npcol) {
    pdc +=
        ctot[0 * lctot + col] + ctot[1 * lctot + col] + ctot[2 * lctot + col];
  }
  return pdc;
}

int get_pdr(int pdc, int klr, int mykl, int nprow, int myrow) {
  auto pdr = pdc;
  auto kl = klr + (mykl % nprow);
  for (auto row = 0; row != myrow; row = (row + 1) % nprow) {
    pdr += kl;
    kl = klr;
  }
  return pdr;
}

/**
 * \brief pjcolのrowインデックスリストを作成
 */
int get_klr(int k, const int indx[], const int indcol[], int pjcol, int indxr[]) {
  int klr = 0;
  for (auto i = 0; i < k; i++) {
    const auto gi = indx[i];
    const auto row = indcol[gi];
    if (row == pjcol) {
      indxr[klr] = i;
      klr += 1;
    }
  }
  return klr;
}

/**
 * \brief 自身のCOLインデクスリストを作成
 */
void set_indxc(int k, const int indx[], const int indcol[], int mycol, int indxc[]) {
  int klc = 0;
  for (auto i = 0; i < k; i++) {
    const auto gi = indx[i];
    const auto col = indcol[gi];
    if (col == mycol) {
      indxc[klc] = i;
      klc += 1;
    }
  }
}

class ComputeArea {
 public:
  int np1;  // 行列の上側
  int np2;  // 行列の下側
};
ComputeArea get_np12(int n, int n1, int np, int myrow, const int indrow[] ) {
  int minrow = n - 1;
  int maxrow = 0;
  int npa = 0;
  int np1 = 0;
  int np2 = 0;
#pragma omp parallel for reduction(min : minrow) reduction(max : maxrow) \
    reduction(+ : npa)
  for (auto i = 0; i < n; i++) {
    if (indrow[i] == myrow) {
      minrow = min(minrow, i);
      maxrow = max(maxrow, i);
      npa += 1;
    }
  }
  if (minrow >= n1) {
    // 上側を持たない
    np1 = 0;
    np2 = npa;
  } else if (maxrow < n1) {
    // 下側を持たない
    np1 = npa;
    np2 = 0;
  } else {
    // 両側を持つ
    np1 = np / 2;
    np2 = npa - np / 2;
  }
  return ComputeArea{np1, np2};
}

template <class Float>
int FS_pdlaed3_C(int k, int n, int n1, Float d[], Float &rho, Float dlamda[],
                 Float w[], int ldq, Float q[],
                 const FS_dividing::bt_node<Float> &subtree, int ldq2,
                 Float q2[], int ldu, Float u[], int indx[], int lctot,
                 const int ctot[], Float q2buf1[], Float q2buf2[], Float z[],
                 Float buf[], int indrow[], int indcol[], int indxc[],
                 int indxr[], int indxcb[], FS_prof::FS_prof &prof) {
#ifdef _DEBUGLOG
  if (FS_libs::FS_get_myrank() == 0) {
    printf("FS_pdlaed3 start\n");
  }
#endif

#if TIMER_PRINT
  prof.start(60);
#endif
  if (k == 0) {
    return 0;
  }
  
  const auto grid_info = subtree.FS_grid_info();
  const auto nprow = grid_info.nprow;
  const auto npcol = grid_info.npcol;
  const auto myrow = grid_info.myrow;
  const auto mycol = grid_info.mycol;
  const auto nblk = subtree.FS_get_NBLK();
  const auto nb = subtree.FS_get_NB();
  const auto np = (nblk / nprow) * nb;
  const auto nq = (nblk / npcol) * nb;

#pragma omp parallel for schedule(static, 1)
  for (auto i = 0; i < n; i += nb) {
    auto row = subtree.FS_info_G1L('R', i).rocsrc;
    auto col = subtree.FS_info_G1L('C', i).rocsrc;
    for (auto j = 0; j < nb; j++) {
      if (i + j < n) {
        indrow[i + j] = row;
        indcol[i + j] = col;
      }
    }
  }

  const auto mykl = ctot[0 * lctot + mycol] + ctot[1 * lctot + mycol] +
                    ctot[2 * lctot + mycol];
  const auto klr = mykl / nprow;
  const auto myklr = (myrow == 0) ? klr + (mykl % nprow) : klr;

  const auto pdc = get_pdc(lctot, ctot, npcol, mycol);
  const auto pdr = get_pdr(pdc, klr, mykl, nprow, myrow);

#pragma omp parallel for
  for (auto i = 0; i < k; i++) {
    dlamda[i] = lapacke::lamc3<Float>(dlamda[i], dlamda[i]) - dlamda[i];
  }

#pragma omp parallel for
  for (auto i = 0; i < 4 * k; i++) {
    buf[i] = FS_const::ZERO<Float>;
  }

#if TIMER_PRINT > 1
  prof.start(61);
#endif
  int sinfo = 0;
  int info = 0;
  if (myklr > 0) {
#pragma omp parallel reduction(max : sinfo)
    {
      std::unique_ptr<Float[]> sz(new Float[k]);
      std::fill_n(sz.get(), k, FS_const::ONE<Float>);
      std::unique_ptr<Float[]> sbuf(new Float[k]);
#pragma omp for schedule(static, 1)
      for (auto i = 0; i < myklr; i++) {
        const int kk = pdr + i;
        const auto buf_index = (pdr - pdc + i);
        Float aux;
        int iinfo =
            lapacke::laed4<Float>(k, kk + 1, dlamda, w, sbuf.get(), rho, aux);
        if (k == 1 || k == 2) {
          buf[buf_index] = FS_const::ZERO<Float>;
          buf[mykl + buf_index] = aux;
        } else {
          auto sbufD_min = sbuf[kk];
          auto sbufB_min = dlamda[kk];
          if ((kk + 1) < k) {
            if (abs(sbuf[kk + 1]) < abs(sbufD_min)) {
              sbufD_min = sbuf[kk + 1];
              sbufB_min = dlamda[kk + 1];
            }
          }
          
          buf[buf_index] = sbufD_min;
          buf[mykl + buf_index] = sbufB_min;
        }
        if (iinfo != 0) {
          sinfo = kk;
          printf("error");
        }
        // ..Compute part of z
#pragma loop nofp_relaxed nofp_contract noeval
        for (auto j = 0; j < k; j++) {
          auto temp = dlamda[j] - dlamda[kk];
          if (j == kk) {
            temp = FS_const::ONE<Float>;
          } else {
            temp = temp;
          }
          sbuf[j] = sbuf[j] / temp;
        }
        for (auto j = 0; j < k; j++) {
          sz[j] *= sbuf[j];
        }
      }
#pragma omp master
      {
        std::copy_n(sz.get(), k, z);
        // count up the flops on the Loewner law's update
        eigen_dc::flops += Float(myklr) * Float(k * 3);
      }
#pragma omp barrier
      for (auto i = 1; i < omp_get_num_threads(); i++) {
        if (omp_get_thread_num() == i) {
          for (auto i = 0; i < k; i++) {
            z[i] *= sz[i];
          }
        }
#pragma omp barrier
      }
    }
    info = sinfo;
  } else {
#pragma omp parallel for
    for (auto i = 0; i < k; i++) {
      z[i] = FS_const::ONE<Float>;
    }
  }
#if TIMER_PRINT > 1
  prof.end(61);
#endif

#if TIMER_PRINT > 1
  prof.start(62);
#endif

#pragma omp parallel for
  for (auto i = 0; i < k; i++) {
    buf[2 * k + i] = z[i];
  }

  MPI_Group_Allreduce<Float>(&buf[2 * k], z, k, FS_const::MPI_TYPE<Float>,
                             MPI_PROD, FS_COMM_WORLD, subtree.MERGE_GROUP_);

#pragma omp parallel for
  for (auto i = 0; i < k; i++) {
    z[i] = abs(sqrt(-z[i])) * ((w[i] >= 0) ? 1 : -1);
  }

  if (mykl > 0) {
    MPI_Group_Allreduce<Float>(buf, &buf[2 * k], 2 * mykl,
                               FS_const::MPI_TYPE<Float>, MPI_SUM,
                               FS_COMM_WORLD, subtree.MERGE_GROUP_X_);
  }

#pragma omp parallel
  {
#pragma omp for
    for (auto i = 0; i < 2 * mykl; i++) {
      buf[i] = FS_const::ZERO<Float>;
    }
#pragma omp for
    for (auto i = 0; i < mykl; i++) {
      buf[pdc + i] = buf[2 * k + i];
      buf[k + pdc + i] = buf[2 * k + mykl + i];
    }
  }

  MPI_Group_Allreduce<Float>(buf, &buf[2 * k], 2 * k, FS_const::MPI_TYPE<Float>,
                             MPI_SUM, FS_COMM_WORLD, subtree.MERGE_GROUP_Y_);

  // Copy of D at the good place

  Float *sbufd = &buf[2 * k];
  Float *sbufb = &buf[3 * k];
  if (k == 1 || k == 2) {
    for (auto i = 0; i < k; i++) {
      const auto gi = indx[i];
      d[gi] = sbufb[i];
    }
  } else {
#pragma omp parallel for
    for (auto i = 0; i < k; i++) {
      const auto gi = indx[i];
      d[gi] = sbufb[i] - sbufd[i];
    }
  }

#if TIMER_PRINT > 1
  prof.end(62);
#endif

  set_indxc(k, indx, indcol, mycol, indxc);

#ifdef _BLOCKING_DGEMM
  // ブロッキングのために列方向昇順の逆引きリストを作成
#pragma omp parallel for
  for (auto j = 0; j < mykl; j++) {
    const auto kk = indxc[j];
    const auto ju = indx[kk];
    const auto jju = subtree.FS_index_G2L('C', ju);
    indxcb[jju] = j;
  }
#endif

  const auto compute_area = get_np12(n, n1, np, myrow, indrow);
  const auto np1 = compute_area.np1;
  const auto np2 = compute_area.np2;

  // Qの初期化
  // デフレーションされた列はQ2からコピー
#if TIMER_PRINT > 1
  prof.start(63);
#endif
#pragma omp parallel
#pragma omp for collapse(2)
  for (auto j = 0; j < mykl; j++) {
    for (auto i = 0; i < np; i++) {
      q[j * ldq + i] = FS_const::ZERO<Float>;
    }
  }
#pragma omp for collapse(2)
  for (auto j = mykl; j < nq; j++) {
    for (auto i = 0; i < np; i++) {
      q[j * ldq + i] = q2[j * ldq2 + i];
    }
  }
#if TIMER_PRINT > 1
  prof.end(63);
#endif

  // Compute eigenvectors of the modified rank-1 modification.

  MPI_Status stat;
  MPI_Request req[2];
  int nrecv = 0;
  int nsend = 0;
  for (auto pj = 0; pj < npcol; pj++) {
    Float *sendq2 = nullptr;
    Float *recvq2 = nullptr;

    // 送受信バッファのポインタ
    if (pj % 2 == 0) {
      sendq2 = q2buf1;
      recvq2 = q2buf2;
    } else {
      sendq2 = q2buf2;
      recvq2 = q2buf1;
    }

    // 演算対象のプロセス列
    const auto pjcol = (mycol + npcol + pj) % npcol;

    // 処理する列数
    const auto mykl = ctot[0 * lctot + mycol] + ctot[1 * lctot + mycol] +
                      ctot[2 * lctot + mycol];
    const auto n12 = ctot[0 * lctot + pjcol] + ctot[1 * lctot + pjcol];
    const auto n23 = ctot[1 * lctot + pjcol] + ctot[2 * lctot + pjcol];

    // pjcolのrowインデックスリストを作成
    const auto klr = get_klr(k, indx, indcol, pjcol, indxr);

    // 前ループの送受信のwaitと展開
    // wait -> recvbuf -> q2
    if (pj != 0 && npcol > 1) {
#if TIMER_PRINT > 1
      prof.start(64);
#endif
      // wait for irecv
      if (nrecv > 0) {
        MPI_Wait(&req[1], &stat);
      }

      // copy recvbuf -> q2
      // 上側
      if (np1 > 0) {
#pragma omp parallel for
        for (auto jq2 = 0; jq2 < n12; jq2++) {
          const auto js = jq2 * np1;
          std::copy_n(&recvq2[js], np1, &q2[jq2 * ldq2 + 0]);
        }
      }
      // 下側
      if (np2 > 0) {
#pragma omp parallel for
        for (auto j = 0; j < n23; j++) {
          const auto jq2 = j + ctot[0 * lctot + pjcol];
          const auto js = n12 * np1 + j * np2;
          std::copy_n(&recvq2[js], np2, &q2[jq2 * ldq2 + np1]);
        }
      }

      // wait for isend
      if (nsend > 0) {
        MPI_Wait(&req[0], &stat);
      }
#if TIMER_PRINT > 1
      prof.end(64);
#endif
    }

    // 送受信バッファのポインタの入れ替え
    if (pj % 2 == 0) {
      sendq2 = q2buf2;
      recvq2 = q2buf1;
    } else {
      sendq2 = q2buf1;
      recvq2 = q2buf2;
    }

    // 次ループの送受信
    // q2 -> sendbuf -> isend
    if (pj != npcol - 1 && npcol > 1) {
#if TIMER_PRINT > 1
      prof.start(65);
#endif
      // copy Q2 -> sendbuf
      // recvbufとsendbufを切り替えながら処理するためsendbufへの格納は初回のみ
      if (pj == 0) {
        // 上側
        if (np1 > 0) {
#pragma omp parallel for
          for (auto jq2 = 0; jq2 < n12; jq2++) {
            const auto js = jq2 * np1;
            std::copy_n(&q2[jq2 * ldq2 + 0], np1, &sendq2[js]);
          }
        }
        // 下側
        if (np2 > 0) {
#pragma omp parallel for
          for (auto j = 0; j < n23; j++) {
            const auto jq2 = j + ctot[0 * lctot + pjcol];
            const auto js = n12 * np1 + j * np2;
            std::copy_n(&q2[jq2 * ldq2 + np1], np2, &sendq2[js]);
          }
        }
      }
      // isend
      nsend = np1 * n12 + np2 * n23;
      if (nsend > 0) {
        const auto dstcol = (mycol + npcol - 1) % npcol;  // 左側に送信
        const auto dest = subtree.group_Y_processranklist_[dstcol];
        MPI_Isend(sendq2, nsend, FS_const::MPI_TYPE<Float>, dest, 1,
                  FS_libs::FS_COMM_WORLD, &req[0]);
      }

      // irecv
      const auto pjcoln = (mycol + npcol + pj + 1) % npcol;
      nrecv = np1 * (ctot[0 * lctot + pjcoln] + ctot[1 * lctot + pjcoln]) +
              np2 * (ctot[1 * lctot + pjcoln] + ctot[2 * lctot + pjcoln]);
      if (nrecv > 0) {
        const auto srccol = (mycol + npcol + 1) % npcol;  // 右側から受信
        const auto source = subtree.group_Y_processranklist_[srccol];
        MPI_Irecv(recvq2, nrecv, FS_const::MPI_TYPE<Float>, source, 1,
                  FS_libs::FS_COMM_WORLD, &req[1]);
      }
#ifdef _MPITEST
      bool flag;
      if (nsend > 0) {
        MPI_Test(&req[0], &flag, &stat);
      }
      if (nrecv > 0) {
        MPI_Test(&req[1], &flag, &stat);
      }
#endif

#if TIMER_PRINT > 1
      prof.end(65);
#endif
    }

    if (mykl != 0) {
#if TIMER_PRINT > 1
      prof.start(66);
#endif
#pragma omp parallel
      {
        std::unique_ptr<Float[]> sbuf(new Float[k]);
#pragma omp for
        for (auto j = 0; j < mykl; j++) {
          const auto kk = indxc[j];
          const auto ju = indx[kk];
          const auto jju = subtree.FS_index_G2L('C', ju);
          Float aux;
          if (k == 1 || k == 2) {
            lapacke::laed4<Float>(k, kk + 1, dlamda, w, sbuf.get(), rho, aux);
          } else {
            for (auto i = 0; i < k; i++) {
              sbuf[i] = dlamda[i] - sbufb[kk];
              sbuf[i] += sbufd[kk];
            }
          }

          if (k == 1 || k == 2) {
            for (auto i = 0; i < klr; i++) {
              const auto kk = indxr[i];
              const auto iu = indx[kk];
              const auto iiu = subtree.FS_index_G2L('C', iu);
              u[jju * ldu + iiu] = sbuf[kk];
            }
            continue;
          }
          for (auto i = 0; i < k; i++) {
            sbuf[i] = z[i] / sbuf[i];
          }
          const auto temp = lapacke::nrm2<Float>(k, sbuf.get(), 1);

          for (auto i = 0; i < klr; i++) {
            const auto kk = indxr[i];
            const auto iu = indx[kk];
            const auto iiu = subtree.FS_index_G2L('C', iu);
            u[jju * ldu + iiu] = sbuf[kk] / temp;
          }
        }
      }

      // Compute the updated eigenvectors.

      // 上側のdgemm
      if (np1 != 0 && n12 != 0) {
#if TIMER_PRINT > 1
        prof.start(67);
#endif
        lapacke::gemm<Float>(CblasNoTrans, CblasNoTrans, np1, mykl, n12,
                             FS_const::ONE<Float>, q2, ldq2, u, ldu,
                             FS_const::ONE<Float>, q, ldq);
        eigen_dc::flops += 2 * double(np1) * double(mykl) * double(n12);
#if TIMER_PRINT > 1
        prof.end(67);
#endif
      }

      // 下側のdgemm
      if (np2 != 0 && n23 != 0) {
        const auto iq2 = np1;
        const auto jq2 = ctot[0 * lctot + pjcol];
        const auto iiu = ctot[0 * lctot + pjcol];
        const auto jju = 0;
        const auto iq = np1;
        const auto jq = 0;

#if TIMER_PRINT > 1
        prof.start(67);
#endif
        lapacke::gemm<Float>(CblasNoTrans, CblasNoTrans, np2, mykl, n23,
                             FS_const::ONE<Float>, &q2[jq2 * ldq2 + iq2], ldq2,
                             &u[jju * ldu + iiu], ldu, FS_const::ONE<Float>,
                             &q[jq * ldq + iq], ldq);
        eigen_dc::flops += 2 * double(np2) * double(mykl) * double(n23);
#if TIMER_PRINT > 1
        prof.end(67);
#endif
      }
#if TIMER_PRINT > 1
      prof.end(66);
#endif
    }
    // ブロッキングなし終了
  }

#if TIMER_PRINT
  prof.end(60);
#endif

#ifdef _DEBUGLOG
  if (FS_libs::FS_get_myrank() == 0) {
    printf("FS_pdlaed3 end. info= %d\n", info);
  }
#endif
  return info;
}
}  // namespace FS_pdlaed3

#endif
