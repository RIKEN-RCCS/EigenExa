#pragma once
#ifndef FS2EIGEN_PDLASRT_HPP
#define FS2EIGEN_PDLASRT_HPP
#include <mpi.h>

#include <algorithm>
#include <cstdio>
#include <memory>
#include <numeric>

#include "../eigen/eigen_libs0.hpp"
#include "FS_dividing.hpp"
#include "FS_libs.hpp"

namespace eigen_FS {
namespace FS2eigen {
using FS_dividing::bt_node;
using std::printf;

template <class Float>
class GpositionValue {
 public:
  int GRow;
  int GCol;
  Float MatrixValue;
};

template <class Float>
class CommBuf {
 public:
  int rank;
  int Ndata;
  MPI_Request req;
  MPI_Status sta;
  bool flag;
  GpositionValue<Float> *bufp;
};

class RANKLIST {
 public:
  int index;
  int *lid;
};

class NrankMaxsize {
 public:
  int nrank;
  int maxsize;
};

inline NrankMaxsize get_nrank_maxsize(int eigen_np,
                              const int comm_send_or_recv_info[]) {
  int send_or_recv_nrank = 0;
  int send_or_recv_maxsize = 1;
#pragma omp parallel for reduction(+ : send_or_recv_nrank) \
    reduction(max : send_or_recv_maxsize)
  for (auto i = 0; i < eigen_np; i++) {
    if (comm_send_or_recv_info[i] != 0) {
      send_or_recv_nrank += 1;
      if (send_or_recv_maxsize < comm_send_or_recv_info[i]) {
        send_or_recv_maxsize = comm_send_or_recv_info[i];
      }
    }
  }
  return NrankMaxsize{send_or_recv_nrank, send_or_recv_maxsize};
}

template <class Float>
void init_send(int np, const int comm_info[], CommBuf<Float> comm_buf[]) {
  int nrank = 0;
  for (auto i = 0; i < np; i++) {
    if (comm_info[i] != 0) {
      comm_buf[nrank].flag = true;
      comm_buf[nrank].rank = static_cast<int>(i);
      comm_buf[nrank].Ndata = (comm_info[i]);
      comm_buf[nrank].bufp = nullptr;
      nrank += 1;
    }
  }
}
template <class Float>
void init_recv(int np, const int comm_info[], CommBuf<Float> comm_buf[]) {
  int nrank = 0;
  for (auto i = 0; i < np; i++) {
    if (comm_info[i] != 0) {
      comm_buf[nrank].flag = true;
      comm_buf[nrank].rank = static_cast<int>(i);
      comm_buf[nrank].Ndata = comm_info[i];
      comm_buf[nrank].bufp = nullptr;
      nrank += 1;
    }
  }
}

template <class Float>
void send(CommBuf<Float> &comm_send_data, GpositionValue<Float> sendbuf[],
          MPI_Comm comm) {
  auto ncount = sizeof(GpositionValue<Float>) * comm_send_data.Ndata;
  auto srank = comm_send_data.rank;
  comm_send_data.flag = false;
  MPI_Send((void *)sendbuf, ncount, MPI_BYTE, srank, 1, comm);
}
template <class Float>
void irecv(CommBuf<Float> &comm_recv_data, GpositionValue<Float> recvbuf[],
           MPI_Comm comm) {
  auto ncount = sizeof(GpositionValue<Float>) * comm_recv_data.Ndata;
  auto rrank = comm_recv_data.rank;
  comm_recv_data.flag = false;
  MPI_Irecv((void *)recvbuf, ncount, MPI_BYTE, rrank, 1, comm,
            &comm_recv_data.req);
}
template <class Float>
int FS_nbroc_max(char comp, int n, const bt_node<Float> &subtree, int FS_nbroc,
                 int FS_myroc) {
  for (auto i = 0; i < FS_nbroc; i++) {
    auto lroc = i;
    auto groc = subtree.FS_index_L2G(comp, lroc, FS_myroc);
    if (n <= groc) {
      return i;
    }
  }
  return FS_nbroc;
}
template <class Float>
int get_FS_nbrow_max(int n, const bt_node<Float> &subtree, int FS_nbrow,
                     int FS_myrow) {
  return FS_nbroc_max('R', n, subtree, FS_nbrow, FS_myrow);
}
template <class Float>
int get_FS_nbcol_max(int n, const bt_node<Float> &subtree, int FS_nbcol,
                     int FS_mycol) {
  return FS_nbroc_max('C', n, subtree, FS_nbcol, FS_mycol);
}

inline int eigen_rank_xy2comm(char grid_major, int x_inod, int y_inod) {
  int nnod, x_nnod, y_nnod;
  eigen_libs0::eigen_get_procs(nnod, x_nnod, y_nnod);

  if (grid_major == 'R') {
    return y_inod + x_inod * y_nnod;
  } else {
    return x_inod + y_inod * x_nnod;
  }
}

/**
 * \brief 送信先が被らないように初回の通信相手を選択する
*/
template <class Float>
int select_first_communicater(int send_nrank, int eigen_np,int eigen_myrank, const CommBuf<Float> comm_send_data[]) {
  int i0 = 0;
  if (0 < send_nrank) {
    i0 = -1;

    for (auto k = 0; k < eigen_np; k++) {
      auto k1 = (k + eigen_myrank) % eigen_np;
      for (auto i = 0; i < send_nrank; i++) {
        if (k1 == comm_send_data[i].rank) {
          i0 = i;
          break;
        }
      }
      if (i0 != -1) {
        break;
      }
    }
  }
  return i0;
}

template <class Float>
int FS2eigen_pdlasrt(int n, Float d[], int ldq, Float q[],
                     const FS_dividing::bt_node<Float> &subtree, int ibuf[],
                     Float rbuf[], GpositionValue<Float> tbuf[], int indx[],
                     FS_prof::FS_prof &prof) {
  double prof_time[40];
  for (auto i = 0; i < 40; i++) {
    prof_time[i] = 0;
  }
#ifdef _DEBUGLOG
  if (FS_libs::FS_get_myrank() == 0) {
    printf("FS2eigen_PDLASRT start.");
  }
#endif
#if TIMER_PRINT
  prof.start(70);
#endif

  int eigen_np, eigen_nprow, eigen_npcol;
  int eigen_myrank, eigen_myrow, eigen_mycol;
  char eigen_grid_major;
  MPI_Comm eigen_comm, eigen_x_comm, eigen_y_comm;
  eigen_libs0::eigen_get_procs(eigen_np, eigen_nprow, eigen_npcol);
  eigen_libs0::eigen_get_id(eigen_myrank, eigen_myrow, eigen_mycol);
  eigen_libs0::eigen_get_grid_major(eigen_grid_major);
  eigen_libs0::eigen_get_comm(eigen_comm, eigen_x_comm, eigen_y_comm);

  int *comm_send_info = new int[eigen_np];
  std::unique_ptr<int[]> comm_recv_info(new int[eigen_np]);

#pragma omp parallel for
  for (auto i = 0; i < eigen_np; i++) {
    comm_send_info[i] = 0;
    comm_recv_info[i] = -1;
  }
  MPI_Bcast(d, n, FS_const::MPI_TYPE<Float>, 0, eigen_comm);


  std::iota(indx, indx + n, 0);
  std::sort(indx, &indx[n], [&d](int i1, int i2) { return d[i1] < d[i2]; });

#pragma omp parallel
  {
#pragma omp for
    for (auto i = 0; i < n; i++) {
      rbuf[i] = d[indx[i]];
    }
#pragma omp for
    for (auto i = 0; i < n; i++) {
      d[i] = rbuf[i];
    }
  }

  auto stime = MPI_Wtime();

  const auto FS_comm_member = FS_libs::is_comm_member();

  const auto FS_grid_info = FS_comm_member ? subtree.FS_grid_info()
                                           : FS_dividing::GridInfo{0, 0, 0, 0};
  const int FS_nblk = FS_comm_member ? subtree.FS_get_NBLK() : 0;
  const int FS_nb = FS_comm_member ? subtree.FS_get_NB() : 0;
  const int FS_nbrow = FS_comm_member ? (FS_nblk / FS_grid_info.nprow) * FS_nb : 0;
  const int FS_nbcol = FS_comm_member ? (FS_nblk / FS_grid_info.npcol) * FS_nb : 0;
  const int FS_nbrow_max =
      FS_comm_member ? get_FS_nbrow_max(n, subtree, FS_nbrow, FS_grid_info.myrow) : 0;
  // プロセス数で割り切れない場合に拡張した領域を除いた計算領域の総数
  const int FS_nbcol_max =
      FS_comm_member ? get_FS_nbcol_max(n, subtree, FS_nbcol, FS_grid_info.mycol) : 0;

  auto etime = MPI_Wtime();
  prof_time[0] = etime - stime;
  stime = etime;

  std::unique_ptr<int[]> lcol2gcol_index;
  std::unique_ptr<int[]> lrow2grow_index;

  if (FS_comm_member) {
    lcol2gcol_index.reset(new int[FS_nbcol]);
#pragma omp parallel for
    for (auto j = 0; j < FS_nbcol; j++) {
      lcol2gcol_index[j] = -1;
    }

    lrow2grow_index.reset(new int[FS_nbrow]);
#pragma omp parallel for
    for (auto i = 0; i < FS_nbrow; i++) {
      lrow2grow_index[i] = -1;
    }

#pragma omp parallel for
    for (auto lrow = 0; lrow < FS_nbrow_max; lrow++) {
      const auto grow = subtree.FS_index_L2G('R', lrow, FS_grid_info.myrow);
      lrow2grow_index[lrow] = grow;
    }

    // row <--> x
    // col <--> y
    // grow,
    // gcolはプロセス行で行列字数を割り切れない場合に対応するために拡張した字数の行列番号を返すので，
    // 拡張した範囲を省くために行列字数を超えたらexitする
#pragma omp parallel for reduction(+ : comm_send_info[0 : eigen_np]) \
    schedule(dynamic, 1)
    for (auto lcol = 0; lcol < FS_nbcol_max; lcol++) {
      // 固有値を並び替える前の列番号を取得
      auto gcol = subtree.FS_index_L2G('C', lcol, FS_grid_info.mycol);

      // 固有値を並び替えた後の列番号に変換
      for (auto k = 0; k < n; k++) {
        if (indx[k] == gcol) {
          gcol = k;
          lcol2gcol_index[lcol] = k;
          break;
        }
      }

      // グローバル情報から再分散後に担当するノードを求める
      const auto pcol = eigen_libs0::eigen_owner_node(gcol, eigen_npcol);
      for (auto lrow = 0; lrow < FS_nbrow_max; lrow++) {
        const auto grow = lrow2grow_index[lrow];
        const auto prow = eigen_libs0::eigen_owner_node(grow, eigen_nprow);
        const auto pn = eigen_rank_xy2comm(eigen_grid_major, prow, pcol);
        comm_send_info[pn] += 1;
      }
    }
  }

  etime = MPI_Wtime();
  prof_time[1] = etime - stime;
  stime = etime;

  MPI_Alltoall(comm_send_info, 1, MPI_INT, comm_recv_info.get(), 1, MPI_INT,
               eigen_comm);

  etime = MPI_Wtime();
  prof_time[2] = etime - stime;
  stime = etime;

  const auto send_nrank_maxsize = get_nrank_maxsize(eigen_np, comm_send_info);
  const auto send_nrank = send_nrank_maxsize.nrank;
  const auto send_maxsize = send_nrank_maxsize.maxsize;

  etime = MPI_Wtime();
  prof_time[3] = etime - stime;
  stime = etime;

  int pointer = 0;
  std::unique_ptr<CommBuf<Float>[]> comm_send_data;
  std::unique_ptr<RANKLIST[]> sendrank_list;
  GpositionValue<Float> *sendbuf;

  if (send_nrank != 0) {
    comm_send_data.reset(new CommBuf<Float>[send_nrank]);
    init_send(eigen_np, comm_send_info, comm_send_data.get());

    sendrank_list.reset(new RANKLIST[send_nrank]);

    for (auto k = 0; k < send_nrank; k++) {
      sendrank_list[k].lid = &ibuf[pointer];
      pointer += comm_send_data[k].Ndata;
      for (auto j = 0; j < comm_send_data[k].Ndata; j++) {
        sendrank_list[k].lid[j] = -1;
      }
    }
    pointer = (pointer + 4) / 4;

    sendbuf = &tbuf[pointer];

    pointer += send_maxsize;
#pragma omp parallel for
    for (auto i = 0; i < send_maxsize; i++) {
      sendbuf[i].GRow = -1;
      sendbuf[i].GCol = -1;
      sendbuf[i].MatrixValue = 0;
    }
  }

  etime = MPI_Wtime();
  prof_time[4] = etime - stime;
  stime = etime;

  const auto recv_nrank_maxsize = get_nrank_maxsize(eigen_np, comm_recv_info.get());
  const auto recv_nrank = recv_nrank_maxsize.nrank;

  etime = MPI_Wtime();
  prof_time[5] = etime - stime;
  stime = etime;

  // 受信バッファの設定
  if (!FS_comm_member) {
    pointer = 0;
  }
  std::unique_ptr<CommBuf<Float>[]> comm_recv_data;
  if (recv_nrank != 0) {
    comm_recv_data.reset(new CommBuf<Float>[recv_nrank]);
    init_recv(eigen_np, comm_recv_info.get(), comm_recv_data.get());
    for (auto k = 0; k < recv_nrank; k++) {
      if (comm_recv_data[k].rank + 1 != eigen_myrank) {
        comm_recv_data[k].bufp = &tbuf[pointer];
        pointer += comm_recv_data[k].Ndata;
        irecv(comm_recv_data[k], comm_recv_data[k].bufp, eigen_comm);
      }
    }
  }

  // 送信データのパック
  if (FS_comm_member) {
    for (auto k = 0; k < send_nrank; k++) {
      auto pn = comm_send_data[k].rank;
      comm_send_info[pn] = k;
      sendrank_list[k].index = 0;
    }

    for (auto lcol = 0; lcol < FS_nbcol_max; lcol++) {
      for (auto lrow = 0; lrow < FS_nbrow_max; lrow++) {
        auto gcol = lcol2gcol_index[lcol];
        auto pcol = eigen_libs0::eigen_owner_node(gcol, eigen_npcol);

        auto grow = lrow2grow_index[lrow];
        auto prow = eigen_libs0::eigen_owner_node(grow, eigen_nprow);

        auto pn = eigen_rank_xy2comm(eigen_grid_major, prow, pcol);

        const auto &info = comm_send_info[pn];

        sendrank_list[info].lid[sendrank_list[info].index] =
            lrow + lcol * FS_nbrow_max;
        sendrank_list[info].index += 1;
      }
    }

    etime = MPI_Wtime();
    prof_time[11] = etime - stime;
    stime = etime;

    const auto i0 = select_first_communicater(send_nrank, eigen_np, eigen_myrank,
                              comm_send_data.get());
    for (auto k0 = i0; k0 < send_nrank + i0; k0++) {
      const auto k = (k0 + 1) % send_nrank;
      if (comm_send_data[k].rank + 1 != eigen_myrank) {
        for (auto j = 0; j < comm_send_data[k].Ndata; j++) {
          const auto lcol = sendrank_list[k].lid[j] / FS_nbrow_max;
          const auto gcol = lcol2gcol_index[lcol];
          const auto lrow = sendrank_list[k].lid[j] % FS_nbrow_max;
          const auto grow = lrow2grow_index[lrow];

          sendbuf[j].GRow = grow;
          sendbuf[j].GCol = gcol;
          sendbuf[j].MatrixValue = q[lcol * ldq + lrow];
        }
        send(comm_send_data[k], sendbuf, eigen_comm);
      }
    }

    // Recv側のWait
    for (auto k0 = 0; k0 < recv_nrank; k0++) {
      if (comm_recv_data[k0].rank + 1 != eigen_myrank) {
        MPI_Wait(&comm_recv_data[k0].req, MPI_STATUS_IGNORE);
      }
    }

    // 送信データのアンパック
    // 送信先が自分自身であればデータコピー
    for (auto k0 = 0; k0 < send_nrank; k0++) {
      if (eigen_myrank == comm_send_data[k0].rank + 1) {
        for (auto i = 0; i < comm_send_data[k0].Ndata; i++) {
          const auto lcol = sendrank_list[k0].lid[i] / FS_nbrow_max;
          const auto gcol = lcol2gcol_index[lcol];
          const auto lrow = sendrank_list[k0].lid[i] % FS_nbrow_max;
          const auto grow = lrow2grow_index[lrow];

          sendbuf[i].GRow = grow;
          sendbuf[i].GCol = gcol;
          sendbuf[i].MatrixValue = q[lcol * ldq + lrow];
        }
        for (auto i = 0; i < comm_send_data[k0].Ndata; i++) {
          const auto gcol = sendbuf[i].GCol;
          const auto grow = sendbuf[i].GRow;
          const auto lcol = eigen_libs0::eigen_translate_g2l(gcol, eigen_npcol);
          const auto lrow = eigen_libs0::eigen_translate_g2l(grow, eigen_nprow);
          q[lcol * ldq + lrow] = sendbuf[i].MatrixValue;
        }
        break;
      }
    }

#pragma omp parallel for schedule(dynamic, 1)
    for (auto k = 0; k < recv_nrank; k++) {
      if (comm_recv_data[k].rank + 1 != eigen_myrank) {
        const auto Ndata = comm_recv_data[k].Ndata;
        for (auto i = 0; i < Ndata; i++) {
          const auto gcol = comm_recv_data[k].bufp[i].GCol;
          const auto grow = comm_recv_data[k].bufp[i].GRow;
          const auto lcol = eigen_libs0::eigen_translate_g2l(gcol, eigen_npcol);
          const auto lrow = eigen_libs0::eigen_translate_g2l(grow, eigen_nprow);
          q[lcol * ldq + lrow] = comm_recv_data[k].bufp[i].MatrixValue;
        }
      }
    }
  } else {
    if (recv_nrank != 0) {
      for (auto k = 0; k < recv_nrank; k++) {
        MPI_Wait(&comm_recv_data[k].req, MPI_STATUS_IGNORE);
        const auto Ndata = comm_recv_data[k].Ndata;
#pragma omp parallel for
        for (auto i = 0; i < Ndata; i++) {
          const auto gcol = comm_recv_data[k].bufp[i].GCol;
          const auto grow = comm_recv_data[k].bufp[i].GRow;
          const auto lcol = eigen_libs0::eigen_translate_g2l(gcol, eigen_npcol);
          const auto lrow = eigen_libs0::eigen_translate_g2l(grow, eigen_nprow);
          q[lcol * ldq + lrow] = comm_recv_data[k].bufp[i].MatrixValue;
        }
      }
    }
  }

  etime = MPI_Wtime();
  prof_time[6] = etime - stime;
  stime = etime;
  delete[] comm_send_info;
  etime = MPI_Wtime();
  prof_time[8] = etime - stime;
  stime = etime;

#if TIMER_PRINT
  prof.end(70);
#endif
#ifdef _DEBUGLOG
  if (FS_libs::FS_get_myrank() == 0) {
    printf("FS_pdlasrt end.");
  }
#endif
  return 0;
}

}  // namespace FS2eigen
}  // namespace eigen_FS

#endif
