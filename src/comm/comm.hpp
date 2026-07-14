#pragma once
#ifndef COMM_HPP
#define COMM_HPP
#define Kahan 0

#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include <type_traits>
#include <vector>
#include <algorithm>
#include "../eigen/eigen_devel.hpp"

namespace comm {
using eigen_devel::repro_reduce;
using eigen_devel::eigen_get_wtime;
using eigen_devel::time_reduce;
using eigen_devel::time_reduce_;
using eigen_devel::items_reduce;
using eigen_devel::counter_reduce_;
using eigen_devel::messages_reduce_;
using eigen_devel::TRD_inod;
using eigen_devel::ZERO;
using eigen_devel::x_inod;
using eigen_devel::x_nnod;
using eigen_devel::y_inod;
using eigen_devel::y_nnod;
using eigen_devel::z_inod;
using eigen_devel::z_nnod;
using eigen_devel::n_common;
using eigen_devel::x_COMM_WORLD;
using eigen_devel::y_COMM_WORLD;
using eigen_devel::z_COMM_WORLD;
using eigen_devel::w_COMM_WORLD;
using eigen_devel::p0_;
using eigen_devel::q0_;
using eigen_devel::items_bcast;
using eigen_devel::items_gather;
using eigen_devel::items_redist;
using eigen_devel::time_redist;
using eigen_devel::time_redist_;
using eigen_devel::counter_redist_;
using eigen_devel::messages_redist_;

constexpr int64_t BCAST_TAG = 100000;
constexpr int64_t BCAST_SEGMENT_SIZE = 4096;
constexpr int64_t BCASTW_TAG = 200000;

#define BCAST_MPI            (-1)
#define BCAST_MPI_SEGMENTED   (1)
#define BCAST_SEQUENTIAL      (2)
#define BCAST_BINOMIAL        (3)
#define BCAST_TRINOMIAL       (4)
#define BCAST_RELAY_ONEDIR    (5)
#define BCAST_RELAY_TWODIRS   (6)
#define BCAST_UNKNOWN         (0)

#define BCAST_ALGORITHM BCAST_MPI

#if BCAST_ALGORITHM == BCAST_UNKNOWN
#error "BCAST_ALGORITHM is undefined or unknown"
#endif

inline void barrier(MPI_Comm comm) {
  MPI_Barrier(comm);
}

template <typename T>
MPI_Datatype get_mpi_datatype();

template <>
inline MPI_Datatype get_mpi_datatype<double>() { return MPI_DOUBLE; }

template <>
inline MPI_Datatype get_mpi_datatype<float>() { return MPI_FLOAT; }

template <>
inline MPI_Datatype get_mpi_datatype<int>() { return MPI_INT; }

template <typename T>
inline void pack(const T* buf, int64_t n, T* buffer, int64_t& ptr) {
  std::copy(buf, buf + n, buffer + ptr);
  ptr += n;
}

template <typename T>
inline void pack1(T buf, T* buffer, int64_t& ptr) {
  buffer[ptr] = buf;
  ptr += 1;
}

template <typename T>
inline void unpack(T* buf, int64_t n, const T* buffer, int64_t& ptr) {
  std::copy(buffer + ptr, buffer + ptr + n, buf);
  ptr += n;
}

template <typename T>
inline void unpack1(T& buf, const T* buffer, int64_t& ptr) {
  buf = buffer[ptr];
  ptr += 1;
}

template <typename T>
inline void send(T* buf, int64_t n, int64_t idest, MPI_Comm icom) {
  MPI_Send(buf, n, get_mpi_datatype<T>(), idest - 1, 1, icom);
}

template <typename T>
inline void send_tagged(T* buf, int64_t n, int64_t idest, int64_t itag, MPI_Comm icom) {
  MPI_Send(buf, n, get_mpi_datatype<T>(), idest - 1, itag, icom);
}

template <typename T>
inline void isend(T* buf, int64_t n, int64_t idest, MPI_Request* ireq, MPI_Comm icom) {
  MPI_Isend(buf, n, get_mpi_datatype<T>(), idest - 1, 1, icom, ireq);
}

template <typename T>
inline void isend_tagged(T* buf, int64_t n, int64_t idest, int64_t itag, MPI_Request* ireq, MPI_Comm icom) {
  MPI_Isend(buf, n, get_mpi_datatype<T>(), idest - 1, itag, icom, ireq);
}

template <typename T>
inline void recv(T* buf, int64_t n, int64_t isrc, MPI_Comm icom) {
  MPI_Recv(buf, n, get_mpi_datatype<T>(), isrc - 1, 1, icom, MPI_STATUS_IGNORE);
}

template <typename T>
inline void irecv(T* buf, int64_t n, int64_t isrc, MPI_Request* ireq, MPI_Comm icom) {
  MPI_Irecv(buf, n, get_mpi_datatype<T>(), isrc - 1, 1, icom, ireq);
}

template <typename T>
inline void irecv_tagged(T* buf, int64_t n, int64_t isrc, int64_t itag, MPI_Request* ireq, MPI_Comm icom) {
  MPI_Irecv(buf, n, get_mpi_datatype<T>(), isrc - 1, itag, icom, ireq);
}

inline void wait(MPI_Request* ireq) {
  MPI_Wait(ireq, MPI_STATUS_IGNORE);
}

inline void waitall(int64_t n, MPI_Request* ireq) {
  MPI_Waitall(n, ireq, MPI_STATUSES_IGNORE);
}

template <typename T>
inline void bcast(T* buf, int64_t n, int iroot, int64_t col_id, MPI_Comm icom) {
  if (n < 1) return;

  int world_size = 0, my_rank = 0;
  MPI_Comm_size(icom, &world_size);
  if (world_size == 1) return;

  MPI_Comm_rank(icom, &my_rank);
  my_rank += 1;

#if BCAST_ALGORITHM == BCAST_MPI
  MPI_Bcast(buf, n, get_mpi_datatype<T>(), iroot - 1, icom);
#endif

#if BCAST_ALGORITHM == BCAST_MPI_SEGMENTED
  for (int64_t i = 0; i < n; i += BCAST_SEGMENT_SIZE) {
    int64_t j = std::min(n - i, BCAST_SEGMENT_SIZE);
    MPI_Bcast(buf + i, j, get_mpi_datatype<T>(), iroot - 1, icom);
  }
#endif

#if BCAST_ALGORITHM == BCAST_SEQUENTIAL
  if (my_rank == iroot) {
    std::vector<MPI_Request> ireq(world_size);
    int64_t j = 0;
    for (int64_t i = 1; i <= world_size; ++i) {
      if (i != iroot) {
        MPI_Isend(buf, n, get_mpi_datatype<T>(), i - 1, i, icom, &ireq[j++]);
      }
    }
    if (j > 0) MPI_Waitall(j, ireq.data(), MPI_STATUSES_IGNORE);
  } else {
    MPI_Recv(buf, n, get_mpi_datatype<T>(), iroot - 1, my_rank, icom, MPI_STATUS_IGNORE);
  }
#endif
#if BCAST_ALGORITHM == BCAST_BINOMIAL
  int64_t local_rank = (world_size + my_rank - iroot) % world_size;
  int64_t i = 1, Glog = 0;
  while (i < world_size) {
    Glog++;
    i *= 2;
  }

  std::vector<MPI_Request> ireq(Glog);
  i = 1;
  int64_t j = 0;
  while (i < world_size) {
    if (local_rank < i) {
      if (local_rank + i < world_size) {
        int64_t k = (my_rank + i - 1) % world_size + 1;
        int64_t tag = BCAST_TAG + k;
        MPI_Isend(buf, n, get_mpi_datatype<T>(), k - 1, tag, icom, &ireq[j++]);
      }
    } else if (local_rank < 2 * i) {
      int64_t k = (world_size + my_rank - i - 1) % world_size + 1;
      int64_t tag = BCAST_TAG + my_rank;
      MPI_Recv(buf, n, get_mpi_datatype<T>(), k - 1, tag, icom, MPI_STATUS_IGNORE);
    }
    i *= 2;
  }
  if (j > 0) MPI_Waitall(j, ireq.data(), MPI_STATUSES_IGNORE);
#endif

#if BCAST_ALGORITHM == BCAST_TRINOMIAL
  int64_t local_rank = (world_size + my_rank - iroot) % world_size;
  int64_t i = 1, Glog = 0;
  while (i < world_size) {
    Glog += 2;
    i *= 3;
  }

  std::vector<MPI_Request> ireq(Glog);
  i = 1;
  int64_t j = 0;
  while (i < world_size) {
    if (local_rank < i) {
      for (int64_t m = 1; m <= 2; ++m) {
        if (local_rank + m * i < world_size) {
          int64_t k = (my_rank + m * i - 1) % world_size + 1;
          int64_t tag = BCAST_TAG + k;
          MPI_Isend(buf, n, get_mpi_datatype<T>(), k - 1, tag, icom, &ireq[j++]);
        }
      }
    } else if (local_rank < 2 * i) {
      int64_t k = (world_size + my_rank - i - 1) % world_size + 1;
      int64_t tag = BCAST_TAG + my_rank;
      MPI_Recv(buf, n, get_mpi_datatype<T>(), k - 1, tag, icom, MPI_STATUS_IGNORE);
    } else if (local_rank < 3 * i) {
      int64_t k = (world_size + my_rank - 2 * i - 1) % world_size + 1;
      int64_t tag = BCAST_TAG + my_rank;
      MPI_Recv(buf, n, get_mpi_datatype<T>(), k - 1, tag, icom, MPI_STATUS_IGNORE);
    }
    i *= 3;
  }
  if (j > 0) MPI_Waitall(j, ireq.data(), MPI_STATUSES_IGNORE);
#endif

#if BCAST_ALGORITHM == BCAST_RELAY_ONEDIR
  int64_t tag = BCAST_TAG;
  int64_t rank_to   = (my_rank - 1 + 1) % world_size + 1;
  int64_t rank_from = (my_rank - 1 - 1 + world_size) % world_size + 1;

  if (my_rank != iroot) {
    MPI_Recv(buf, n, get_mpi_datatype<T>(), rank_from - 1, tag, icom, MPI_STATUS_IGNORE);
  }
  if (rank_to != iroot) {
    MPI_Send(buf, n, get_mpi_datatype<T>(), rank_to - 1, tag, icom);
  }
#endif

#if BCAST_ALGORITHM == BCAST_RELAY_TWODIRS
  int64_t tag = BCAST_TAG;
  int64_t Ghalf = world_size / 2;
  int64_t local_rank = (world_size + my_rank - iroot) % world_size;
  if (local_rank >= Ghalf) local_rank -= world_size;

  if (my_rank == iroot) {
    std::vector<MPI_Request> ireq(2);
    int64_t i = (my_rank - 1 + 1) % world_size + 1;
    int64_t j = (my_rank - 1 - 1 + world_size) % world_size + 1;
    if (i != j) {
      MPI_Isend(buf, n, get_mpi_datatype<T>(), i - 1, tag, icom, &ireq[0]);
      MPI_Isend(buf, n, get_mpi_datatype<T>(), j - 1, tag, icom, &ireq[1]);
      MPI_Waitall(2, ireq.data(), MPI_STATUSES_IGNORE);
    } else {
      MPI_Isend(buf, n, get_mpi_datatype<T>(), i - 1, tag, icom, &ireq[0]);
      MPI_Wait(&ireq[0], MPI_STATUS_IGNORE);
    }
  } else {
    int64_t rank_to, rank_from;
    if (local_rank > 0) {
      rank_to   = (my_rank - 1 + 1) % world_size + 1;
      rank_from = (my_rank - 1 - 1 + world_size) % world_size + 1;
    } else {
      rank_to   = (my_rank - 1 - 1 + world_size) % world_size + 1;
      rank_from = (my_rank - 1 + 1) % world_size + 1;
    }

    MPI_Recv(buf, n, get_mpi_datatype<T>(), rank_from - 1, tag, icom, MPI_STATUS_IGNORE);

    int64_t k = (world_size + rank_to - iroot) % world_size;
    if (k >= Ghalf) k -= world_size;
    if (k * local_rank > 0) {
      MPI_Send(buf, n, get_mpi_datatype<T>(), rank_to - 1, tag, icom);
    }
  }
#endif
}

template <typename T>
void bcastw(T* buf, int64_t n, int64_t iroot, int64_t lda, int64_t lpx,
                T* buffer, int64_t col_id, MPI_Comm icom) {
  if (lpx == 1) {
    bcast<T>(buf, n, iroot, col_id, icom);
    return;
  }
  int world_size = 0, my_rank = 0;
  MPI_Comm_size(icom, &world_size);
  MPI_Comm_rank(icom, &my_rank);
  my_rank = my_rank + 1;

  if (lpx <= 4) {
    for (int64_t i = 0; i < lpx; ++i) {
      int k = (iroot - 1 + i) % world_size + 1;
      bcast<T>(&buf[i * lda], n, k, col_id, icom);
    }
    return;
  }

#if TIMER_PRINT
  double timer = eigen_get_wtime();
#endif

  std::vector<int> rcounts(world_size, 0);
  std::vector<int> displs(world_size, 0);

  for (int64_t i = 1; i <= world_size; ++i) {
    int64_t k = (i - iroot + world_size) % world_size + 1;
    if (k <= lpx) {
      rcounts[i-1] = n;
      displs[i-1] = (k - 1) * lda;
    } else {
      rcounts[i-1] = 0;
      displs[i-1] = 0;
    }
  }

  int64_t k = (my_rank - iroot + world_size) % world_size + 1;
  int64_t j = (k <= lpx) ? n : 0;
  if (k > lpx) k = 1;

  std::copy(&buf[(k - 1) * lda], &buf[(k - 1) * lda] + n, buffer);

  MPI_Allgatherv(buffer, static_cast<int>(j), get_mpi_datatype<T>(),
                 buf, rcounts.data(), displs.data(),
                 get_mpi_datatype<T>(), icom);

#if TIMER_PRINT
  timer = eigen_get_wtime() - timer;
  time_bcast += timer;
  if (col_id >= 1 && col_id <= items_bcast) {
    time_bcast_[col_id]     += timer;
    counter_bcast_[col_id]  += n * lpx;
    messages_bcast_[col_id] += 1;
  }
#endif
}

inline void floor_log2(int64_t& n) {
  if (n <= 0) {
    // If n<=0, no log2(n) exists. So, assign n to dummy(-9999).
    n = -9999;
    return;
  }

  int64_t power_of_2 = 2;
  for (int64_t power = 0; power <= 29; ++power) {
    if (n < power_of_2) {
      n = power;
      return;
    }
    power_of_2 *= 2;
  }
  // 4 Byte integer type variable is definitely less than 2**31
  n = 30;
}

template <typename T>
inline void ALLREDUCE_sum(MPI_Comm comm, int64_t n, T* buff, T* buff0) {
  int ierr = 0, nprocess = 0, myrank = 0;
  int irank;

#ifdef Kahan
  std::vector<T> buff_c(n, 0.0);
#endif

  ierr = MPI_Comm_size(comm, &nprocess);
  ierr = MPI_Comm_rank(comm, &myrank);

  if (nprocess > 1) {
    if (myrank == 0) {
      std::vector<T> rbuff(n * nprocess, 0.0);
      std::vector<MPI_Request> req_recv(nprocess);

      for (int64_t i = 1; i < nprocess; ++i) {
        irank = i;
        MPI_Irecv(&rbuff[irank * n], n, get_mpi_datatype<T>(), irank, irank, comm, &req_recv[irank]);
      }

      for (int64_t i = 1; i < nprocess; ++i) {
        MPI_Wait(&req_recv[i], MPI_STATUS_IGNORE);
      }


      for (int64_t j = 1; j < nprocess; ++j) {
	irank = j;
#ifdef Kahan
        for (int64_t i = 0; i < n; ++i) {
          T yy = rbuff[i + irank * n] - buff_c[i];
          T tt = buff[i] + yy;
          buff_c[i] = (tt - buff[i]) - yy;
          buff[i] = tt;
        }
#else
        for (int64_t i = 0; i < n; ++i) {
          buff[i] += rbuff[i + irank * n];
        }
#endif
      }
    } else {  // myrank /= 0
      MPI_Request req_send;
      MPI_Isend(buff, n, get_mpi_datatype<T>(), 0, myrank, comm, &req_send);
      MPI_Wait(&req_send, MPI_STATUS_IGNORE);
    }

    ierr = MPI_Bcast(buff, n, get_mpi_datatype<T>(), 0, comm);
  }

  std::copy(buff, buff + n, buff0);
}

template <typename T>
inline void ALLREDUCE_binary_sum(MPI_Comm comm1, int64_t s, T* R, T* R0) {
  int irank, nprocess, ierr, ss_int4;
  ierr = MPI_Comm_size(comm1, &nprocess);
  ierr = MPI_Comm_rank(comm1, &irank);

  std::vector<T> R_c(s, 0.0), sbuff(R, R + s), rbuff(s);

  int64_t irank0 = 0, irank0_master = 0;

  int64_t istage = static_cast<int64_t>(nprocess);
  floor_log2(istage);

  int64_t irank_stage = irank;

  for (int64_t i = istage; i >= 0; --i) {
    if (nprocess < (irank0_master + (1 << i))) continue;

    if (irank_stage < (1 << i)) {
      // ----
      // irank0_master --> 集約された三角行列を持つランク
      // irank0_master - irank0 --> 集約した三角行列の送信先ランク
      // -----
      int64_t stage = 1, j = 0;
      while (stage < (1 << i)) {
        ++j;
        stage = 1 << j;
        if (irank % stage == 0) {
          if (irank == nprocess - 1) {
            int64_t nrk = irank0_master;
            stage = 1 << i;
	    // 送信先が自プロセスでない場合
            if (nrk != irank) {
              ss_int4 = static_cast<int>(s);
              MPI_Request req;
              MPI_Isend(sbuff.data(), ss_int4, get_mpi_datatype<T>(), nrk, nrk, comm1, &req);
              MPI_Wait(&req, MPI_STATUS_IGNORE);
            }
            break;
          } else {
            int64_t nrk = irank + stage / 2;
	    // 受信先がある場合
            if (nrk < nprocess) {
              ss_int4 = static_cast<int>(s);
              MPI_Request req;
              MPI_Irecv(rbuff.data(), ss_int4, get_mpi_datatype<T>(), nrk, irank, comm1, &req);
              MPI_Wait(&req, MPI_STATUS_IGNORE);

              for (int64_t k = 0; k < s; ++k) {
#ifdef Kahan
                double yy = rbuff[k] - R_c[k];
                double tt = R[k] + yy;
                R_c[k] = (tt - R[k]) - yy;
                R[k] = tt;
                sbuff[k] = R[k];
#else
                R[k] += rbuff[k];
                sbuff[k] = R[k];
#endif
              }
            }
          }
        } else {
          int64_t nrk = irank - stage / 2;
	  ss_int4 = static_cast<int>(s);
          MPI_Request req;
          MPI_Isend(sbuff.data(), ss_int4, get_mpi_datatype<T>(), nrk, nrk, comm1, &req);
          MPI_Wait(&req, MPI_STATUS_IGNORE);
          break;
        }
      }
      // -----
      // ここで1ステージ上のツリーに集約を行う

      if (irank0_master == irank && irank != 0) {
        if (irank + (1 << i) < nprocess) {
          int64_t nrk = irank + (1 << i);
	  ss_int4 = static_cast<int>(s);
          MPI_Request req;
          MPI_Irecv(rbuff.data(), ss_int4, get_mpi_datatype<T>(), nrk, irank, comm1, &req);
          MPI_Wait(&req, MPI_STATUS_IGNORE);

          for (int64_t k = 0; k < s; ++k) {
#ifdef Kahan
            double yy = rbuff[k] - R_c[k];
            double tt = R[k] + yy;
            R_c[k] = (tt - R[k]) - yy;
            R[k] = tt;
            sbuff[k] = R[k];
#else
            R[k] += rbuff[k];
            sbuff[k] = R[k];
#endif
          }
        }

        int64_t nrk = irank0_master - irank0;
        ss_int4 = static_cast<int>(s);
        MPI_Request req;
        MPI_Isend(sbuff.data(), ss_int4, get_mpi_datatype<T>(), nrk, nrk, comm1, &req);
        MPI_Wait(&req, MPI_STATUS_IGNORE);
      }

      if (irank == 0 && ((1 << i) < nprocess)) {
        int64_t nrk = irank + (1 << i);
        ss_int4 = static_cast<int>(s);
        MPI_Request req;
        MPI_Irecv(rbuff.data(), ss_int4, get_mpi_datatype<T>(), nrk, irank, comm1, &req);
        MPI_Wait(&req, MPI_STATUS_IGNORE);

        for (int64_t k = 0; k < s; ++k) {
#ifdef Kahan
          double yy = rbuff[k] - R_c[k];
          double tt = R[k] + yy;
          R_c[k] = (tt - R[k]) - yy;
          R[k] = tt;
          sbuff[k] = R[k];
#else
          R[k] += rbuff[k];
          sbuff[k] = R[k];
#endif
        }
      }

      break;
    } else {
      irank0 = 1 << i;
      irank0_master += irank0;
      irank_stage %= (1 << i);
    }
  }

  ss_int4 = static_cast<int>(s);
  ierr = MPI_Bcast(R, ss_int4, get_mpi_datatype<T>(), 0, comm1);
  std::copy(R, R + s, R0);
}

template <typename T>
inline void ALLREDUCE_binary_prod(MPI_Comm comm1, int64_t s, T* R, T* R0) {
  int irank, nprocess, ierr;
  ierr = MPI_Comm_size(comm1, &nprocess);
  ierr = MPI_Comm_rank(comm1, &irank);

  std::vector<T> sbuff(R, R + s), rbuff(s);
  int64_t istage = nprocess;
  floor_log2(istage);

  int64_t irank_stage = irank;
  int64_t irank0 = 0, irank0_master = 0;

  for (int64_t i = istage; i >= 0; --i) {
    if (nprocess < (irank0_master + (1 << i))) continue;

    if (irank_stage < (1 << i)) {
      int64_t stage = 1, j = 0;
      while (stage < (1 << i)) {
        ++j;
        stage = 1 << j;

        if (irank % stage == 0) {
          if (irank == nprocess - 1) {
            int64_t nrk = irank0_master;
            if (nrk != irank) {
              MPI_Request req;
              MPI_Isend(sbuff.data(), s, get_mpi_datatype<T>(), nrk, nrk, comm1, &req);
              MPI_Wait(&req, MPI_STATUS_IGNORE);
            }
            break;
          } else {
            int64_t nrk = irank + stage / 2;
            if (nrk < nprocess) {
              MPI_Request req;
              MPI_Irecv(rbuff.data(), s, get_mpi_datatype<T>(), nrk, irank, comm1, &req);
              MPI_Wait(&req, MPI_STATUS_IGNORE);
              for (int64_t k = 0; k < s; ++k) {
                R[k] *= rbuff[k];
                sbuff[k] = R[k];
              }
            }
          }
        } else {
          int64_t nrk = irank - stage / 2;
          MPI_Request req;
          MPI_Isend(sbuff.data(), s, get_mpi_datatype<T>(), nrk, nrk, comm1, &req);
          MPI_Wait(&req, MPI_STATUS_IGNORE);
          break;
        }
      }

      if (irank0_master == irank && irank != 0) {
        if (irank + (1 << i) < nprocess) {
          int64_t nrk = irank + (1 << i);
          MPI_Request req;
          MPI_Irecv(rbuff.data(), s, get_mpi_datatype<T>(), nrk, irank, comm1, &req);
          MPI_Wait(&req, MPI_STATUS_IGNORE);
          for (int64_t k = 0; k < s; ++k) {
            R[k] *= rbuff[k];
            sbuff[k] = R[k];
          }
        }

        int64_t nrk = irank0_master - irank0;
        MPI_Request req;
        MPI_Isend(sbuff.data(), s, get_mpi_datatype<T>(), nrk, nrk, comm1, &req);
        MPI_Wait(&req, MPI_STATUS_IGNORE);
      }

      if (irank == 0 && ((1 << i) < nprocess)) {
        int64_t nrk = irank + (1 << i);
        MPI_Request req;
        MPI_Irecv(rbuff.data(), s, get_mpi_datatype<T>(), nrk, irank, comm1, &req);
        MPI_Wait(&req, MPI_STATUS_IGNORE);
        for (int64_t k = 0; k < s; ++k) {
          R[k] *= rbuff[k];
          sbuff[k] = R[k];
        }
      }

      break;
    } else {
      irank0 = 1 << i;
      irank0_master += irank0;
      irank_stage %= (1 << i);
    }
  }

  ierr = MPI_Bcast(R, s, get_mpi_datatype<T>(), 0, comm1);
  std::copy(R, R + s, R0);
}

template <typename T>
inline void reduce(T* buf, T* wrk, int64_t n, int64_t col_id, MPI_Comm icom) {
  int64_t nnod = 0;

#if TIMER_PRINT
  double timer = eigen_get_wtime();
#endif

#if DEBUG
  if (TRD_inod == 1) {
    std::cout << "Reduction[" << n << "] :: " << std::boolalpha << repro_reduce << std::endl;
  }
#endif

  int innod = static_cast<int>(nnod);
  MPI_Comm_size(icom, &innod);
  nnod = static_cast<int64_t>(innod);

#if defined(__INTEL_LLVM_COMPILER)
  if (repro_reduce && n >= std::max(innod, 32)) {
#else
  if (repro_reduce) {
#endif
    MPI_Allreduce(buf, wrk, static_cast<int>(n), get_mpi_datatype<T>(), MPI_SUM, icom);
  } else if (n >= std::max(innod, 32)) {
    ALLREDUCE_binary_sum<T>(icom, n, buf, wrk);
  } else {
    ALLREDUCE_sum<T>(icom, n, buf, wrk);
  }

  std::copy(wrk, wrk + n, buf);

#if TIMER_PRINT
  timer = eigen_get_wtime() - timer;
  time_reduce += timer;
  if (col_id >= 1 && col_id <= items_reduce) {
    time_reduce_[col_id]     += timer;
    counter_reduce_[col_id]  += n;
    messages_reduce_[col_id] += 1;
  }
#endif
}

template <typename T>
inline void allgather(const T* buf, T* wrk, int64_t n, int64_t col_id, MPI_Comm icom) {
  int64_t ierr = 0;

#if TIMER_PRINT
  double timer = eigen_get_wtime();
#endif

#if __IBM_REGISTER_VARS && 0
  int64_t my_rank = 0, world_size = 0;
  MPI_Comm_size(icom, &world_size);
  MPI_Comm_rank(icom, &my_rank);

  if (my_rank == 0) {
    std::copy(buf, buf + n, wrk);
    for (int64_t i = 1; i < world_size; ++i) {
      MPI_Recv(wrk + i * n, n, get_mpi_datatype<T>(), i, 1, icom, MPI_STATUS_IGNORE);
    }
  } else {
    MPI_Send(buf, n, get_mpi_datatype<T>(), 0, 1, icom);
  }

  MPI_Bcast(wrk, world_size * n, get_mpi_datatype<T>(), 0, icom);
#else
  MPI_Allgather(buf, n, get_mpi_datatype<T>(),
                wrk, n, get_mpi_datatype<T>(), icom);
#endif

#if TIMER_PRINT
  timer = eigen_get_wtime() - timer;
  time_gather += timer;

  if (col_id > 0) {
    if (col_id >= 1 && col_id <= items_gather) {
      time_gather_[col_id]     += timer;
      counter_gather_[col_id]  += n;
      messages_gather_[col_id] += 1;
    }
  } else {
    int64_t cid = -col_id;
    if (cid >= 1 && cid <= items_gather + items_redist) {
      time_gather_[cid]     += timer;
      counter_gather_[cid]  += n;
      messages_gather_[cid] += 1;
    }
  }
#endif
}

inline void print_bcast_algorithm() {
#if BCAST_ALGORITHM == BCAST_MPI
  std::cout << "Bcast algorithm is MPI" << std::endl;
#endif
#if BCAST_ALGORITHM == BCAST_MPI_SEGMENTED
  std::cout << "Bcast algorithm is MPI_SEGMENTED" << std::endl;
#endif
#if BCAST_ALGORITHM == BCAST_SEQUENTIAL
  std::cout << "Bcast algorithm is SEQUENTIAL" << std::endl;
#endif
#if BCAST_ALGORITHM == BCAST_BINOMIAL
  std::cout << "Bcast algorithm is BINOMIAL" << std::endl;
#endif
#if BCAST_ALGORITHM == BCAST_TRINOMIAL
  std::cout << "Bcast algorithm is TRINOMIAL" << std::endl;
#endif
#if BCAST_ALGORITHM == BCAST_RELAY_ONEDIR
  std::cout << "Bcast algorithm is RELAY_ONE_DIR" << std::endl;
#endif
#if BCAST_ALGORITHM == BCAST_RELAY_TWODIRS
  std::cout << "Bcast algorithm is RELAY_TWO_DIRS" << std::endl;
#endif
}

inline void print_reduce_algorithm() {
  std::cout << "Reduce algorithm is MPI" << std::endl;
}

inline void print_gather_algorithm() {
  std::cout << "Gather algorithm is MPI" << std::endl;
}

template <typename T>
inline void datacast(T* u_y, const T* u_x, T* u_t, T* u_s, int64_t n, int64_t col_id) {
#if TIMER_PRINT
  double timer = eigen_get_wtime();
  double timer_excl1 = 0.0;
  double timer_excl2 = 0.0;
#endif

  if (x_nnod == 1 && y_nnod == 1) {
    std::copy(u_x, u_x + n, u_y);
    return;
  }

  int64_t n_x = (n - 1) / x_nnod + 1;
  int64_t n_y = (n - 1) / y_nnod + 1;

  if (x_nnod == 1) {
    for (int64_t i = 1; i <= n_y; ++i) {
      int64_t j = y_inod + y_nnod * (i - 1);
      u_y[i - 1] = u_x[j - 1];
    }
    return;
  }

  if (x_nnod == y_nnod) {
    if (x_inod == y_inod) {
      std::copy(u_x, u_x + n_y, u_y);
    }
    int64_t k = col_id + items_bcast;
    bcast<T>(u_y, n_y, static_cast<int>(y_inod), -k, x_COMM_WORLD);
    return;
  }

  if (p0_[x_inod - 1] > 0) {
    int64_t x_snod = x_nnod / n_common;
    int64_t y_snod = y_nnod / n_common;

    int64_t nx = (n_x - 1) / y_snod + 1;
    int64_t ny = (n_x - p0_[x_inod - 1]) / y_snod + 1;

    for (int64_t i = 1; i <= ny; ++i) {
      int64_t j = p0_[x_inod - 1] + y_snod * (i - 1);
      u_t[i - 1] = u_x[j - 1];
    }
    if (nx > ny) {
      std::fill(u_t + ny, u_t + nx, ZERO);
    }

#if TIMER_PRINT
    timer_excl2 = eigen_get_wtime();
#endif
    int64_t k = col_id + items_gather;
    allgather<T>(u_t, u_s, nx, -k, w_COMM_WORLD);
#if TIMER_PRINT
    timer_excl2 = eigen_get_wtime() - timer_excl2;
#endif

    for (int64_t ic = 0; ic < x_snod; ++ic) {
      int64_t his_rank = (x_inod - 1 + x_nnod + ic * n_common) % x_nnod + 1;
      if (n_x >= p0_[his_rank - 1] && n_y >= q0_[his_rank - 1]) {
        int64_t his_local = (his_rank - 1) / n_common + 1;
        int64_t ny2 = (n_x - p0_[his_rank - 1]) / y_snod + 1;
        int64_t nz2 = (n_y - q0_[his_rank - 1]) / x_snod + 1;
        for (int64_t i = 1; i <= std::min(ny2, nz2); ++i) {
          int64_t k2 = q0_[his_rank - 1] + x_snod * (i - 1);
          int64_t j = (his_local - 1) * nx + i;
          u_y[k2 - 1] = u_s[j - 1];
        }
      }
    }

    if (y_inod > (n - 1) % y_nnod + 1) {
      u_y[n_y - 1] = ZERO;
    }

#if TIMER_PRINT
    timer_excl1 = eigen_get_wtime();
#endif
    if (z_nnod > 1) {
      int64_t kb = col_id + items_bcast;
      bcast<T>(u_y, n_y, static_cast<int>(z_inod), -kb, z_COMM_WORLD);
    }
#if TIMER_PRINT
    timer_excl1 = eigen_get_wtime() - timer_excl1;
#endif

  } else {
    int64_t i = (y_inod - 1) % n_common;
    int64_t j = (x_inod - 1) % n_common;
    int64_t ic = (j - i + n_common * x_nnod * y_nnod) % n_common;
    int64_t his_rank = (x_inod - 1 + x_nnod * y_nnod - ic) % x_nnod + 1;
    int64_t his_local = (his_rank - 1) % n_common + 1;

#if TIMER_PRINT
    timer_excl1 = eigen_get_wtime();
#endif
    int64_t kb = col_id + items_bcast;
    bcast<T>(u_y, n_y, static_cast<int>(his_local), -kb, z_COMM_WORLD);
#if TIMER_PRINT
    timer_excl1 = eigen_get_wtime() - timer_excl1;
    timer_excl2 = 0.0;
#endif
  }

#if TIMER_PRINT
  timer = eigen_get_wtime() - timer;
  timer = timer - timer_excl1 - timer_excl2;
  time_redist += timer;
  if (col_id >= 1 && col_id <= items_redist) {
    time_redist_[col_id]     += timer;
    counter_redist_[col_id]  += n;
    messages_redist_[col_id] += 1;
  }
#endif
}

template <typename T>
inline void datacast2(T* ur_y, T* ui_y, const T* ur_x, const T* ui_x, T* u_t, T* u_s, int64_t n, int64_t col_id) {
#if TIMER_PRINT
  double timer = eigen_get_wtime();
  double timer_excl1 = 0.0;
  double timer_excl2 = 0.0;
#endif

  if (x_nnod == 1 && y_nnod == 1) {
    std::copy(ur_x, ur_x + n, ur_y);
    std::copy(ui_x, ui_x + n, ui_y);
    return;
  }

  int64_t n_x = (n - 1) / x_nnod + 1;
  int64_t n_y = (n - 1) / y_nnod + 1;

  if (x_nnod == 1) {
    for (int64_t i = 1; i <= n_y; ++i) {
      int64_t j = y_inod + y_nnod * (i - 1);
      ur_y[i - 1] = ur_x[j - 1];
      ui_y[i - 1] = ui_x[j - 1];
    }
    return;
  }

  if (x_nnod == y_nnod) {
    if (x_inod == y_inod) {
      std::copy(ur_x, ur_x + n_y, u_t);
      std::copy(ui_x, ui_x + n_y, u_t + n_y);
    }
    int64_t k = col_id + items_bcast;
    bcast<T>(u_t, 2 * n_y, static_cast<int>(y_inod), -k, x_COMM_WORLD);
    std::copy(u_t, u_t + n_y, ur_y);
    std::copy(u_t + n_y, u_t + 2 * n_y, ui_y);
    return;
  }

  if (p0_[x_inod - 1] > 0) {
    int64_t x_snod = x_nnod / n_common;
    int64_t y_snod = y_nnod / n_common;

    int64_t nx = (n_x - 1) / y_snod + 1;
    int64_t ny = (n_x - p0_[x_inod - 1]) / y_snod + 1;

    for (int64_t i = 1; i <= ny; ++i) {
      int64_t j = p0_[x_inod - 1] + y_snod * (i - 1);
      u_t[i - 1]      = ur_x[j - 1];
      u_t[nx + i - 1] = ui_x[j - 1];
    }
    if (nx > ny) {
      std::fill(u_t + ny,      u_t + nx,      ZERO);
      std::fill(u_t + nx + ny, u_t + nx + nx, ZERO);
    }

#if TIMER_PRINT
    timer_excl2 = eigen_get_wtime();
#endif
    int64_t k = col_id + items_gather;
    allgather<T>(u_t, u_s, 2 * nx, -k, w_COMM_WORLD);
#if TIMER_PRINT
    timer_excl2 = eigen_get_wtime() - timer_excl2;
#endif

    for (int64_t ic = 0; ic < x_snod; ++ic) {
      int64_t his_rank = (x_inod - 1 + x_nnod + ic * n_common) % x_nnod + 1;
      if (n_x >= p0_[his_rank - 1] && n_y >= q0_[his_rank - 1]) {
        int64_t his_local = (his_rank - 1) / n_common + 1;
        int64_t ny2 = (n_x - p0_[his_rank - 1]) / y_snod + 1;
        int64_t nz2 = (n_y - q0_[his_rank - 1]) / x_snod + 1;
        for (int64_t i = 1; i <= std::min(ny2, nz2); ++i) {
          int64_t k2 = q0_[his_rank - 1] + x_snod * (i - 1);
          int64_t j = 2 * (his_local - 1) * nx + i;
          ur_y[k2 - 1] = u_s[j - 1];
          ui_y[k2 - 1] = u_s[nx + j - 1];
        }
      }
    }

    if (y_inod > (n - 1) % y_nnod + 1) {
      ur_y[n_y - 1] = ZERO;
      ui_y[n_y - 1] = ZERO;
    }

    std::copy(ur_y, ur_y + n_y, u_s);
    std::copy(ui_y, ui_y + n_y, u_s + n_y);

#if TIMER_PRINT
    timer_excl1 = eigen_get_wtime();
#endif
    if (z_nnod > 1) {
      int64_t kb = col_id + items_bcast;
      bcast<T>(u_s, 2 * n_y, static_cast<int>(z_inod), -kb, z_COMM_WORLD);
    }
#if TIMER_PRINT
    timer_excl1 = eigen_get_wtime() - timer_excl1;
#endif

  } else {
    int64_t i = (y_inod - 1) % n_common;
    int64_t j = (x_inod - 1) % n_common;
    int64_t ic = (j - i + n_common * x_nnod * y_nnod) % n_common;
    int64_t his_rank = (x_inod - 1 + x_nnod * y_nnod - ic) % x_nnod + 1;
    int64_t his_local = (his_rank - 1) % n_common + 1;

#if TIMER_PRINT
    timer_excl1 = eigen_get_wtime();
#endif
    int64_t kb = col_id + items_bcast;
    bcast<T>(u_s, 2 * n_y, static_cast<int>(his_local), -kb, z_COMM_WORLD);
#if TIMER_PRINT
    timer_excl1 = eigen_get_wtime() - timer_excl1;
    timer_excl2 = 0.0;
#endif

    std::copy(u_s, u_s + n_y, ur_y);
    std::copy(u_s + n_y, u_s + 2 * n_y, ui_y);
  }

#if TIMER_PRINT
  timer = eigen_get_wtime() - timer;
  timer = timer - timer_excl1 - timer_excl2;
  time_redist += timer;
  if (col_id >= 1 && col_id <= items_redist) {
    time_redist_[col_id]     += timer;
    counter_redist_[col_id]  += 2 * n;
    messages_redist_[col_id] += 1;
  }
#endif
}

template <typename T>
inline void datacastx(int64_t nk, T* u_y, const T* u_x, int64_t ldv, T* u_t, T* u_s, int64_t n, int64_t col_id) {
#if TIMER_PRINT
  double timer = eigen_get_wtime();
  double timer_excl1 = 0.0;
  double timer_excl2 = 0.0;
#endif

  if (x_nnod == 1 && y_nnod == 1) {
    for (int64_t kk = 0; kk < nk; ++kk) {
      std::copy(u_x + kk * ldv, u_x + kk * ldv + n, u_y + kk * ldv);
    }
    return;
  }

  int64_t n_x = (n - 1) / x_nnod + 1;
  int64_t n_y = (n - 1) / y_nnod + 1;

  if (x_nnod == 1) {
    for (int64_t kk = 0; kk < nk; ++kk) {
      for (int64_t i = 1; i <= n_y; ++i) {
        int64_t j = y_inod + y_nnod * (i - 1);
        u_y[kk * ldv + i - 1] = u_x[kk * ldv + j - 1];
      }
    }
    return;
  }

  if (x_nnod == y_nnod) {
    if (x_inod == y_inod) {
      for (int64_t kk = 0; kk < nk; ++kk) {
        std::copy(u_x + kk * ldv, u_x + kk * ldv + n_y, u_t + kk * n_y);
      }
    }
    int64_t k = col_id + items_bcast;
    bcast<T>(u_t, n_y * nk, static_cast<int>(y_inod), -k, x_COMM_WORLD);
    for (int64_t kk = 0; kk < nk; ++kk) {
      std::copy(u_t + kk * n_y, u_t + kk * n_y + n_y, u_y + kk * ldv);
    }
    return;
  }

  if (p0_[x_inod - 1] > 0) {
    int64_t x_snod = x_nnod / n_common;
    int64_t y_snod = y_nnod / n_common;

    int64_t nx = (n_x - 1) / y_snod + 1;
    int64_t ny = (n_x - p0_[x_inod - 1]) / y_snod + 1;

    for (int64_t kk = 0; kk < nk; ++kk) {
      for (int64_t i = 1; i <= ny; ++i) {
        int64_t j = p0_[x_inod - 1] + y_snod * (i - 1);
        u_t[kk * nx + i - 1] = u_x[kk * ldv + j - 1];
      }
      if (nx > ny) {
        std::fill(u_t + kk * nx + ny, u_t + kk * nx + nx, ZERO);
      }
    }

#if TIMER_PRINT
    timer_excl2 = eigen_get_wtime();
#endif
    int64_t k = col_id + items_gather;
    allgather<T>(u_t, u_s, nk * nx, -k, w_COMM_WORLD);
#if TIMER_PRINT
    timer_excl2 = eigen_get_wtime() - timer_excl2;
#endif

    for (int64_t ic = 0; ic < x_snod; ++ic) {
      int64_t his_rank = (x_inod - 1 + x_nnod + ic * n_common) % x_nnod + 1;
      if (n_x >= p0_[his_rank - 1] && n_y >= q0_[his_rank - 1]) {
        int64_t his_local = (his_rank - 1) / n_common + 1;
        int64_t ny2 = (n_x - p0_[his_rank - 1]) / y_snod + 1;
        int64_t nz2 = (n_y - q0_[his_rank - 1]) / x_snod + 1;
        for (int64_t kk = 0; kk < nk; ++kk) {
          for (int64_t i = 1; i <= std::min(ny2, nz2); ++i) {
            int64_t k2 = q0_[his_rank - 1] + x_snod * (i - 1);
            int64_t j = (his_local - 1) * (nk * nx) + i + kk * nx;
            u_y[kk * ldv + k2 - 1] = u_s[j - 1];
          }
        }
      }
    }

    if (y_inod > (n - 1) % y_nnod + 1) {
      for (int64_t kk = 0; kk < nk; ++kk) {
        u_y[kk * ldv + n_y - 1] = ZERO;
      }
    }

    if (z_nnod > 1) {
      for (int64_t kk = 0; kk < nk; ++kk) {
        std::copy(u_y + kk * ldv, u_y + kk * ldv + n_y, u_s + kk * n_y);
      }
    }

#if TIMER_PRINT
    timer_excl1 = eigen_get_wtime();
#endif
    if (z_nnod > 1) {
      int64_t kb = col_id + items_bcast;
      bcast<T>(u_s, nk * n_y, static_cast<int>(z_inod), -kb, z_COMM_WORLD);
    }
#if TIMER_PRINT
    timer_excl1 = eigen_get_wtime() - timer_excl1;
#endif

  } else {
    int64_t i = (y_inod - 1) % n_common;
    int64_t j = (x_inod - 1) % n_common;
    int64_t ic = (j - i + n_common * x_nnod * y_nnod) % n_common;
    int64_t his_rank = (x_inod - 1 + x_nnod * y_nnod - ic) % x_nnod + 1;
    int64_t his_local = (his_rank - 1) % n_common + 1;

#if TIMER_PRINT
    timer_excl1 = eigen_get_wtime();
#endif
    int64_t kb = col_id + items_bcast;
    bcast<T>(u_s, nk * n_y, static_cast<int>(his_local), -kb, z_COMM_WORLD);
#if TIMER_PRINT
    timer_excl1 = eigen_get_wtime() - timer_excl1;
    timer_excl2 = 0.0;
#endif

    for (int64_t kk = 0; kk < nk; ++kk) {
      std::copy(u_s + kk * n_y, u_s + kk * n_y + n_y, u_y + kk * ldv);
    }
  }

#if TIMER_PRINT
  timer = eigen_get_wtime() - timer;
  timer = timer - timer_excl1 - timer_excl2;
  time_redist += timer;
  if (col_id >= 1 && col_id <= items_redist) {
    time_redist_[col_id]     += timer;
    counter_redist_[col_id]  += 2 * n;
    messages_redist_[col_id] += 1;
  }
#endif
}

}
#endif

