#pragma once
#ifndef FS_LIBS_HPP
#define FS_LIBS_HPP

#include <mpi.h>

#include "../eigen/eigen_libs0.hpp"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace FS_libs {
using eigen_libs0::eigen_free0;
using eigen_libs0::eigen_get_comm;
using eigen_libs0::eigen_get_id;
using eigen_libs0::eigen_get_procs;
using eigen_libs0::eigen_init0;
using std::log2;
using std::min;
using std::pow;
using std::printf;

extern MPI_Comm FS_COMM_WORLD;
extern int FS_MYRANK;
extern bool FS_COMM_MEMBER;
extern MPI_Group FS_GROUP;

class process_grid {
 public:
  int nnod, x_nnod, y_nnod;
  int inod, x_inod, y_inod;
};
extern process_grid FS_node;
class version_t {
 public:
  int Major_Version;
  int Minor_Version;
  int Patch_Level;
  char date[32];
  char vcode[32];
};
constexpr version_t FS_Version = {1, 1, 0, "Mar 31, 2019", "FS proto"};

extern char FS_GRID_major;

inline void FS_get_version(int &version, char date[32] = nullptr,
                    char vcode[32] = nullptr) {
  version = FS_Version.Major_Version * 100 + FS_Version.Minor_Version * 10 +
            FS_Version.Patch_Level;
  if (date != nullptr) {
    std::snprintf(date, 32, "%s\n", FS_Version.date);
  }
  if (vcode != nullptr) {
    std::snprintf(vcode, 32, "%s\n", FS_Version.vcode);
  }
}


inline void FS_show_version() {
  int id, x_id, y_id;
  eigen_get_id(id, x_id, y_id);
  const auto i = min(26, FS_Version.Patch_Level);
  const auto patchlevel = " abcdefghijklmnopqrstuvwxyz*"[i + 1];

  char version[256];
  snprintf(version, 256, "%d.%d%c", FS_Version.Major_Version,
           FS_Version.Minor_Version, patchlevel);

  if (id == 1) {
    printf("## FS version (%s) / (%s) / (%s)\n", version, FS_Version.date,
           FS_Version.vcode);
  }
}

inline void FS_init_cartesian(char GRID_major, int nnod, int inod) {
  auto x_nnod = int(sqrt(double(nnod)));
  int i = 1;
  const auto k = (nnod % i == 0) ? i : 1;

  while (true) {
    if (x_nnod <= k) {
      break;
    }
    if (x_nnod % k == 0 && nnod % x_nnod == 0) {
      break;
    }
    x_nnod -= 1;
  }

  const auto y_nnod = nnod / x_nnod;

  int x_inod, y_inod;
  if (GRID_major == 'R') {
    // row-major
    x_inod = (inod - 1) / y_nnod + 1;
    y_inod = (inod - 1) % y_nnod + 1;
  } else {
    // column-major
    x_inod = (inod - 1) % x_nnod + 1;
    y_inod = (inod - 1) / x_nnod + 1;
  }

  FS_node.nnod = nnod;
  FS_node.x_nnod = x_nnod;
  FS_node.y_nnod = y_nnod;
  FS_node.inod = inod;
  FS_node.x_inod = x_inod;
  FS_node.y_inod = y_inod;
}

inline void FS_init(MPI_Comm comm = MPI_COMM_WORLD, char order = 'C') {
  auto comm0 = comm;
  FS_GRID_major = order;
  if (FS_GRID_major == 'R' || FS_GRID_major == 'r') {
    FS_GRID_major = 'R';
  } else {
    FS_GRID_major = 'C';
  }
  eigen_init0(comm0, FS_GRID_major);
  MPI_Comm eigen_comm, eigen_x_comm, eigen_y_comm;
  eigen_get_comm(eigen_comm, eigen_x_comm, eigen_y_comm);

  // FS_COMM_WORLDの設定
  int nnod, x_nnod, y_nnod;
  int inod, x_inod, y_inod;
  eigen_get_procs(nnod, x_nnod, y_nnod);
  eigen_get_id(inod, x_inod, y_inod);

  const auto p = static_cast<int>(log2(nnod));
  int color = 0;
  if (inod <= pow(2, p)) {
    color = 0;
    FS_COMM_MEMBER = true;
  } else {
    color = 1;
    FS_COMM_MEMBER = false;
  }

  MPI_Comm_split(eigen_comm, color, inod, &FS_COMM_WORLD);

  if (FS_COMM_MEMBER) {
    MPI_Comm_rank(FS_COMM_WORLD, &FS_MYRANK);
    MPI_Comm_group(FS_COMM_WORLD, &FS_GROUP);

    MPI_Comm_size(FS_COMM_WORLD, &nnod);
    FS_init_cartesian(FS_GRID_major, nnod, FS_MYRANK + 1);
  } else {
    FS_MYRANK = -1;
    FS_node.nnod = -1;
    FS_node.x_nnod = -1;
    FS_node.y_nnod = -1;
    FS_node.inod = -1;
    FS_node.x_inod = -1;
    FS_node.y_inod = -1;
  }
}

inline void FS_free() {
  eigen_free0();
  MPI_Comm_free(&FS_COMM_WORLD);
}


inline void FS_get_procs(int *nnod, int *x_nnod, int *y_nnod) {
  *nnod = FS_node.nnod;
  *x_nnod = FS_node.x_nnod;
  *y_nnod = FS_node.y_nnod;
}

inline void FS_get_id(int *inod, int *x_inod, int *y_inod) {
  *inod = FS_node.inod;
  *x_inod = FS_node.x_inod;
  *y_inod = FS_node.y_inod;
}

inline void FS_get_matdims(int n, int &nx, int &ny) {
  int nnod, x_nnod, y_nnod;
  FS_get_procs(&nnod, &x_nnod, &y_nnod);
  int n1 = n / nnod;
  if (n % nnod != 0) {
    n1 += 1;
  }
  nx = n1 * (nnod / x_nnod);
  ny = n1 * (nnod / y_nnod);
}

inline void FS_WorkSize(int n, int &lwork, int &liwork) {
  int nnod, x_nnod, y_nnod;
  int np, nq;
  FS_get_procs(&nnod, &x_nnod, &y_nnod);
  FS_get_matdims(n, np, nq);

  lwork = 1 + 7 * n + 3 * np * nq + nq * nq;
  liwork = 1 + 8 * n + 2 * 4 * y_nnod;
}
inline void FS_get_grid_major(char *Major) {
  *Major = FS_GRID_major;
  return;
}
inline int FS_get_myrank() { return FS_MYRANK; }

inline MPI_Group FS_get_group() { return FS_GROUP; }

inline bool is_comm_member() { return FS_COMM_MEMBER; }

struct Nod {
  int nod;
  int x;
  int y;
};

inline Nod get_procs() {
  Nod nnod = {};
  FS_get_procs(&nnod.nod, &nnod.x, &nnod.y);
  return nnod;
}

inline Nod get_id() {
  Nod inod = {};
  FS_get_id(&inod.nod, &inod.x, &inod.y);
  return inod;
}

inline char get_grid_major() {
  char major = 0;
  FS_get_grid_major(&major);
  return major;
}

inline int get_myrank() { return FS_get_myrank(); }

}  // namespace FS_libs
#endif
