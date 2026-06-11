#pragma once
#ifndef FS_MERGE_D_HPP
#define FS_MERGE_D_HPP

#include <algorithm>

#include "../blas.hpp"
#include "FS_const.hpp"
#include "FS_dividing.hpp"

namespace eigen_FS {
using std::min;
template <class Float>
void FS_merge_d(int n, const Float d[],
                const FS_dividing::bt_node<Float> &subtree, Float d_out[]) {
  std::fill_n(d_out, n, FS_const::ZERO<Float>);

  const auto grid_info = subtree.FS_grid_info();
  const auto NB = subtree.FS_get_NB();

  if (subtree.direction_horizontal_ == 1) {
    const auto col = subtree.FS_info_G1L('C', 0).rocsrc;
    if (col == grid_info.mycol) {
      for (int i = 0; i < n; i += NB) {
        const auto row = subtree.FS_info_G1L('R', i).rocsrc;
        if (row == grid_info.myrow) {
          const auto NB1 = min(n, i + NB) - i;
          std::copy_n(&d[i], NB1, &d_out[i]);
        }
      }
    }
  } else {
    const auto row = subtree.FS_info_G1L('R', 0).rocsrc;
    if (row == grid_info.myrow) {
      for (int j = 0; j < n; j += NB) {
        const auto col = subtree.FS_info_G1L('C', j).rocsrc;
        if (col == grid_info.mycol) {
          const auto NB1 = min(n, j + NB) - j;
          std::copy_n(&d[j], NB1, &d_out[j]);
        }
      }
    }
  }
}
}  // namespace eigen_FS

#endif
