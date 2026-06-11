#pragma once
#ifndef FS_REDUCE_ZD_HPP
#define FS_REDUCE_ZD_HPP

#include <cstdio>

#include "FS_const.hpp"
#include "FS_dividing.hpp"
#include "FS_libs.hpp"
#include "FS_prof.hpp"
#include "MPI_Allreduce_group.hpp"

namespace eigen_FS {
using FS_dividing::bt_node;
using FS_libs::FS_COMM_WORLD;
using FS_prof::FS_prof;
using std::printf;
template <class Float>
void FS_reduce_zd(int n, const bt_node<Float> &subtree, Float work[], Float z[],
                  Float d[], FS_prof &prof) {
#ifdef _DEBUGLOG
  if (FS_libs::FS_get_myrank() == 0) {
    printf("FS_reduce_zd start.\n");
  }
#endif
#if TIMER_PRINT
  prof.start(45);
#endif
#if TIMER_PRINT
#ifdef PROF_DETAIL
  prof.start(46);
  MPI_Barrier(subtree.MERGE_COMM);
  prof.end(46);
  prof.start(47);
#endif
#endif

  // reduce
  MPI_Group_Allreduce(work, z, n * 2, FS_const::MPI_TYPE<Float>, MPI_SUM,
                      FS_COMM_WORLD, subtree.MERGE_GROUP_);

#pragma omp parallel for
  for (auto j = 0; j < n; j++) {
    d[j] = z[j + n];
  }

#if TIMER_PRINT
#ifdef PROF_DETAIL
  prof.end(47);
#endif
#endif

#if TIMER_PRINT
  prof.end(45);
#endif

#ifdef _DEBUGLOG
  if (FS_libs::FS_get_myrank() == 0) {
    printf("FS_reduce_zd end.\n");
  }
#endif
}

}  // namespace eigen_FS

#endif
