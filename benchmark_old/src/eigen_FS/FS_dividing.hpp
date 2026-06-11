#pragma once
#ifndef FS_DIVIDING
#define FS_DIVIDING
#include <mpi.h>

#include <algorithm>
#include <cstdio>
#include <memory>
#include <utility>

#include "FS_libs.hpp"
#include "FS_prof.hpp"

namespace FS_dividing {
using std::printf;

struct g1l {
  int l_index;
  int rocsrc;
};

struct GridIndex {
  int row;
  int col;
};

class GridInfo {
 public:
  int nprow;
  int npcol;
  int myrow;
  int mycol;
};

template <class Float>
class bt_node {
 public:
  int bt_id = 0;
  int layer_ = 0;
  bool direction_horizontal_ = true;
  int nstart_ = 0;
  int nend_ = 1;
  int nend_active_ = 1;
  int proc_istart_ = 0;  //  process start number of direction i
  int proc_iend_ = 1;    //  process end   number of direction i
  int proc_jstart_ = 0;  //  process start number of direction j
  int proc_jend_ = 1;    //  process end   number of direction j
  int block_start_ = 0;  //  merge block start number(refer to procs_i/j)
  int block_end_ = 1;    //  merge block end   number(refer to procs_i/j)
  std::unique_ptr<bt_node<Float>[]> sub_bt_node_;  //  sub tree node
  bt_node<Float> *parent_node_;  //  parent node
  int *procs_i_;                 //  process No. list of row
  int *procs_j_;                 //  process No. list of column
  int nnod_ = 0;                 //  nprocs of communicator
  int x_nnod_ = 0;               //  nprocs of X direction communicator
  int y_nnod_ = 0;               //  nprocs of Y direction communicator
  int inod_ = 0;                 //  inod in MERGE_COMM(1～)
  int x_inod_ = 0;               //  x_inod in MERGE_COMM_X(1～)
  int y_inod_ = 0;               //  y_inod in MERGE_COMM_Y(1～)
  int div_bit_ = -1;             //  bit stream of divided direction
  int div_nbit_ = 0;             //  number of dights of div_bit

  MPI_Group MERGE_GROUP_ = MPI_GROUP_NULL;        //  MERGE_COMM group
  MPI_Group MERGE_GROUP_X_ = MPI_GROUP_NULL;      //  MERGE_COMM_X group
  MPI_Group MERGE_GROUP_Y_ = MPI_GROUP_NULL;      //  MERGE_COMM_Y group
  std::unique_ptr<int[]> group_processranklist_;  //  list to convert from group
                                                  //  rank to communicator rank
  std::unique_ptr<int[]> group_X_processranklist_;
  std::unique_ptr<int[]> group_Y_processranklist_;

 public:
  int FS_dividing(int n, Float d[], const Float e[],
                  std::unique_ptr<bool[]> hint, FS_prof::FS_prof &prof);

  int FS_dividing_recursive(int n, Float d[], const Float e[],
                            std::unique_ptr<bool[]> &hint,
                            FS_prof::FS_prof &prof, int bt_id = 0);

  void dividing_setBitStream();

  inline bool FS_node_included() const {
    const FS_libs::Nod inod = FS_libs::get_id();
    if (inod.x < this->proc_istart_ || inod.x > this->proc_iend_) {
      return false;
    }
    if (inod.y < this->proc_jstart_ || inod.y > this->proc_jend_) {
      return false;
    }
    return true;
  }

  std::pair<const bt_node<Float> *, int> FS_dividing_getleaf(int info) const {
    const FS_libs::Nod inod = FS_libs::get_id();
    if (inod.x - 1 == this->proc_istart_ && inod.x == this->proc_iend_ &&
        inod.y - 1 == this->proc_jstart_ && inod.y == this->proc_jend_) {
      return std::make_pair(this, info);
    }

    if (this->sub_bt_node_ != nullptr) {
      for (int i = 0; i < 2; i++) {
        const bt_node<Float> &sub_bt_node = this->sub_bt_node_[i];
        const auto leaf_and_info = sub_bt_node.FS_dividing_getleaf(info);
        if (leaf_and_info.first != nullptr) {
          return leaf_and_info;
        }
        info = leaf_and_info.second;
      }
    }

    if (this->layer_ == 0) {
      return std::make_pair(nullptr, 9999);
    }
    return std::make_pair(nullptr, info);
  }

  void FS_create_merge_comm(FS_prof::FS_prof &prof);

  void FS_create_mergeXY_group();

  void FS_create_merge_comm_recursive();

  void FS_dividing_free();

  inline int FS_get_N() const { return this->nend_ - this->nstart_; }

  inline int FS_get_N_active() const {
    return this->nend_active_ - this->nstart_;
  }

  inline int FS_get_NBLK() const {
    return this->block_end_ - this->block_start_;
  }

  inline int FS_get_NB() const {
    return this->FS_get_N() / this->FS_get_NBLK();
  }

  GridIndex FS_get_QTOP() const;

  inline const GridInfo FS_grid_info() const {
    return GridInfo{.nprow = this->x_nnod_,
                    .npcol = this->y_nnod_,
                    .myrow = this->x_inod_,
                    .mycol = this->y_inod_};
  }

  g1l FS_info_G1L(char comp, int g_index) const;

  inline const GridIndex FS_info_G2L(int gr_index, int gc_index) const {
    const auto r = this->FS_info_G1L('R', gr_index);
    const auto i = this->FS_info_G1L('I', gc_index);
    return GridIndex{.row = r.l_index, .col = i.l_index};
  }

  inline int FS_index_G2L(char comp, int g_index) const {
    const auto g1l = this->FS_info_G1L(comp, g_index);
    return g1l.l_index;
  }

  int FS_index_L2G(char comp, int l_indx, int my_roc) const;

  void print_tree() const;

  void print_node() const;
};

void FS_create_hint(bool[]);

void bitprint(int kout, int title, int ibit, int nbit);

}  // namespace FS_dividing

namespace FS_dividing {
using std::abs;
using std::max;
using std::min;
using std::printf;

inline void FS_create_hint(bool hint[]) {
  FS_libs::Nod nnod = FS_libs::get_procs();

#if _TREEDIV == 1
  for (int layer = 0; nnod.x * nnod.y >= 1; layer++) {
    if (nnod.y >= nnod.x) {
      hint[layer] = false;
      nnod.y = nnod.y / 2;
    } else {
      hint[layer] = true;
      nnod.x = nnod.x / 2;
    }
  }

#elif _TREEDIV == 2
  for (int layer = 0; nnod.x * nnod.y >= 1; layer++) {
    if (nnod.x >= 2) {
      hint[layer] = true;
      nnod.x = nnod.x / 2;
    } else {
      hint[layer] = false;
      nnod.y = nnod.y / 2;
    }
  }
#elif _TREEDIV == 3
  for (int layer = 0; nnod.x * nnod.y >= 1; layer++) {
    if (nnod.y >= 2) {
      hint[layer] = false;
      nnod.y = nnod.y / 2;
    } else {
      hint[layer] = true;
      nnod.x = nnod.x / 2;
    }
  }
#else
  for (int layer = 0; nnod.x * nnod.y >= 1; layer++) {
    if (nnod.x >= nnod.y) {
      hint[layer] = true;
      nnod.x = nnod.x / 2;
    } else {
      hint[layer] = false;
      nnod.y = nnod.y / 2;
    }
  }
#endif

#ifdef _DEBUGLOG
  if (FS_libs::FS_get_myrank() == 0) {
    // TODO :: write(*,'(a,1000(1x,l))') "procdiv=",hint(1:layer)
    // Fortranの書式がわからないのでスキップする
  }
#endif
}

/**
 * \brief main routine of dividing tree
 */
template <class Float>
int bt_node<Float>::FS_dividing(int n, Float d[], const Float e[],
                                std::unique_ptr<bool[]> hint,
                                FS_prof::FS_prof &prof) {
#if TIMER_PRINT
  prof.start(21);
#endif
  const FS_libs::Nod nnod = FS_libs::get_procs();

  const auto Next = ((n % nnod.nod) == 0) ? n : ((n / nnod.nod + 1) * nnod.nod);

  this->layer_ = 0;  // root nodeは0, sub_node[0]と[1]は1
  this->direction_horizontal_ = hint[this->layer_];
  this->nstart_ = 0;  // このbt_nodeの担当する範囲であり、ループの開始地点
  this->nend_ = Next;
  this->nend_active_ = n;
  this->proc_istart_ = 0;
  this->proc_iend_ = nnod.x;
  this->proc_jstart_ = 0;
  this->proc_jend_ = nnod.y;
  this->block_start_ = 0;
  this->block_end_ = nnod.nod;
  this->procs_i_ = new int[nnod.nod];
  this->procs_j_ = new int[nnod.nod];
  std::fill_n(this->procs_i_, nnod.nod, -1);
  std::fill_n(this->procs_j_, nnod.nod, -1);
  this->parent_node_ = nullptr;
  this->bt_id = 0;
  const int info = this->FS_dividing_recursive(n, d, e, hint, prof);

#if TIMER_PRINT > 1
  prof.start(22);
#endif

  this->FS_create_merge_comm(prof);
#if TIMER_PRINT > 1
  prof.end(22);
#endif

#if TIMER_PRINT > 1
  prof.start(23);
  MPI_Barrier(FS_libs::FS_COMM_WORLD);
  prof.end(23);
#endif

#if TIMER_PRINT
  prof.end(21);
#endif

  return info;
}

/**
 * \brief create sub-tree recursive
 */
template <class Float>
int bt_node<Float>::FS_dividing_recursive(const int n, Float d[],
                                          const Float e[],
                                          std::unique_ptr<bool[]> &hint,
                                          FS_prof::FS_prof &prof, int bt_id) {
  const FS_libs::Nod nnod = FS_libs::get_procs();
  const auto x_lnod = this->proc_iend_ - this->proc_istart_;
  const auto y_lnod = this->proc_jend_ - this->proc_jstart_;
  const auto lnod = x_lnod * y_lnod;

  if (lnod == 0) {
    return -1;
  } else if (lnod == 1) {
    int i;
    for (i = bt_id; i < nnod.nod; i++) {
      if (this->procs_i_[i] < 0) {
        this->procs_i_[i] = this->proc_istart_;  // ここでのみ変更される
        this->procs_j_[i] = this->proc_jstart_;  // ここでのみ変更される

        this->block_start_ = i;
        this->block_end_ = i + 1;
        const auto nend = this->nend_;
        if (nend < n) {
          d[nend - 1] = d[nend - 1] - abs(e[nend - 1]);
          d[nend] = d[nend] - abs(e[nend - 1]);
        }
        break;
      }
    }
    this->bt_id = i;
    return 0;
  } else if ((hint[this->layer_] && x_lnod == 1) ||
             (!hint[this->layer_] && y_lnod == 1)) {
    return -2;
  } else {
    int info;
    this->sub_bt_node_.reset(new bt_node<Float>[2]);

    for (int i = 0; i < 2; i++) {
      const auto this_nstart = this->nstart_;
      const auto this_nend = this->nend_;
      const auto this_nstep = (this_nend - this_nstart) / 2;
      bt_node<Float> &subptr = this->sub_bt_node_[i];

      subptr.layer_ = this->layer_ + 1;
      subptr.direction_horizontal_ = hint[subptr.layer_];
      subptr.nstart_ = this_nstart + i * this_nstep;
      subptr.nend_ = this_nend - (1 - i) * this_nstep;
      subptr.nend_active_ = max(min(subptr.nend_, n), subptr.nstart_);
      subptr.sub_bt_node_.reset(nullptr);
      subptr.parent_node_ = this;
      subptr.procs_i_ = this->procs_i_;
      subptr.procs_j_ = this->procs_j_;

      subptr.proc_istart_ = this->proc_istart_;
      subptr.proc_iend_ = this->proc_iend_;
      subptr.proc_jstart_ = this->proc_jstart_;
      subptr.proc_jend_ = this->proc_jend_;
      if (hint[this->layer_]) {
        const auto proc_i_step = (this->proc_iend_ - this->proc_istart_) / 2;
        subptr.proc_istart_ += i * proc_i_step;
        subptr.proc_iend_ -= (1 - i) * proc_i_step;
      } else {
        const auto proc_j_step = (this->proc_jend_ - this->proc_jstart_) / 2;
        subptr.proc_jstart_ += i * proc_j_step;
        subptr.proc_jend_ -= (1 - i) * proc_j_step;
      }

      info = subptr.FS_dividing_recursive(n, d, e, hint, prof);

      if (info != 0) {
        return info;
      }
    }

    this->block_start_ = min(this->sub_bt_node_[0].block_start_,
                             this->sub_bt_node_[1].block_start_);
    this->block_end_ =
        max(this->sub_bt_node_[0].block_end_, this->sub_bt_node_[1].block_end_);
    this->dividing_setBitStream();
    return info;
  }
}

/**
 * \brief set bit stream of tree dividing direction of all child node to leaf
 */
template <class Float>
void bt_node<Float>::dividing_setBitStream() {
  if (this->sub_bt_node_ == nullptr) {
    return;
  }

  this->div_bit_ = 0;
  if (this->direction_horizontal_ == false) {
    this->div_bit_ |= (1 << 0);  // IBSET
  }

  this->div_nbit_ = 1;
  bt_node<Float> *node = &(this->sub_bt_node_[0]);
  while (node->sub_bt_node_ != nullptr) {
    this->div_bit_ <<= 1;  // ISHIFT
    if (node->direction_horizontal_ == false) {
      this->div_bit_ |= (1 << 0);  // IBSET
    }
    this->div_nbit_ += 1;
    node = &(node->sub_bt_node_[0]);
  }
}

/**
 * \brief create local merge X,Y group
 */
template <class Float>
void bt_node<Float>::FS_create_mergeXY_group() {
  const char order = FS_libs::get_grid_major();
  const FS_libs::Nod inod = FS_libs::get_id();
  const FS_libs::Nod nnod = FS_libs::get_procs();
  std::unique_ptr<int[]> ranklist_group(new int[nnod.nod]);

  {
    const auto proc_istart = this->proc_istart_;
    const auto proc_iend = this->proc_iend_;

    const auto ii_nrank = [&]() mutable {
      int ii_incl_nrank = 0;
      if (order == 'R') {
        const auto jj = inod.y - 1;
        for (auto ii = proc_istart; ii < proc_iend; ii++) {
          ranklist_group[ii_incl_nrank] = (ii)*nnod.y + (jj);
          ii_incl_nrank += 1;
        }
      } else {
        const auto jj = inod.y - 1;
        for (auto ii = proc_istart; ii < proc_iend; ii++) {
          ranklist_group[ii_incl_nrank] = (jj)*nnod.x + (ii);
          ii_incl_nrank += 1;
        }
      }
      return ii_incl_nrank;
    }();

    MPI_Group_incl(FS_libs::FS_get_group(), ii_nrank, ranklist_group.get(),
                   &this->MERGE_GROUP_X_);

    std::iota(ranklist_group.get(), ranklist_group.get() + ii_nrank, 0);

    this->group_X_processranklist_.reset(new int[ii_nrank]);

    MPI_Group_translate_ranks(this->MERGE_GROUP_X_, ii_nrank,
                              ranklist_group.get(), FS_libs::FS_get_group(),
                              this->group_X_processranklist_.get());
  }
  // j <--> y
  {
    const auto proc_jstart = this->proc_jstart_;
    const auto proc_jend = this->proc_jend_;
    const auto jj_nrank = [&]() mutable {
      int jj_incl_nrank = 0;
      if (order == 'R') {
        const auto ii = inod.x - 1;
        for (auto jj = proc_jstart; jj < proc_jend; jj++) {
          ranklist_group[jj_incl_nrank] = (ii)*nnod.y + (jj);
          jj_incl_nrank += 1;
        }
      } else {
        const auto ii = inod.x - 1;
        for (auto jj = proc_jstart; jj < proc_jend; jj++) {
          ranklist_group[jj_incl_nrank] = (jj)*nnod.x + (ii);
          jj_incl_nrank += 1;
        }
      }
      return jj_incl_nrank;
    }();

    MPI_Group_incl(FS_libs::FS_get_group(), jj_nrank, ranklist_group.get(),
                   &this->MERGE_GROUP_Y_);

    std::iota(ranklist_group.get(), ranklist_group.get() + jj_nrank, 0);

    this->group_Y_processranklist_.reset( new int[jj_nrank]);
    MPI_Group_translate_ranks(this->MERGE_GROUP_Y_, jj_nrank,
                              ranklist_group.get(), FS_libs::FS_get_group(),
                              this->group_Y_processranklist_.get());
  }
}

template <class Float>
void bt_node<Float>::FS_create_merge_comm(FS_prof::FS_prof &prof) {
  const FS_libs::Nod inod = FS_libs::get_id();
  const FS_libs::Nod nnod = FS_libs::get_procs();

  this->MERGE_GROUP_ = FS_libs::FS_get_group();

  if (this->MERGE_GROUP_ != MPI_GROUP_NULL) {
    this->inod_ = inod.nod - 1;
    this->x_inod_ = inod.x - 1;
    this->y_inod_ = inod.y - 1;
    this->nnod_ = nnod.nod;
    this->x_nnod_ = nnod.x;
    this->y_nnod_ = nnod.y;
#if TIMER_PRINT > 1
    prof.start(24);
#endif
    this->FS_create_mergeXY_group();
#if TIMER_PRINT > 1
    prof.end(24);
#endif
  }

#if TIMER_PRINT > 1
  prof.start(25);
#endif
  this->FS_create_merge_comm_recursive();
#if TIMER_PRINT > 1
  prof.end(25);
#endif
}

template <class Float>
void bt_node<Float>::FS_create_merge_comm_recursive() {
  const FS_libs::Nod nnod = FS_libs::get_procs();
  const char order = FS_libs::get_grid_major();

  if (this->sub_bt_node_ == nullptr) {
    return;
  }

  for (int n = 0; n < 2; n++) {
    bt_node<Float> &node = this->sub_bt_node_[n];

    const auto ni = node.proc_iend_ - node.proc_istart_;
    const auto nj = node.proc_jend_ - node.proc_jstart_;
    const auto ranklist_nrank = ni * nj;
    std::unique_ptr<int[]> ranklist(new int[ranklist_nrank]);
    std::unique_ptr<int[]> ranklist_group(new int[ranklist_nrank]);

    const auto nrank = [&]() mutable {
      int incl_nrank = 0;
      if (order == 'R') {
        for (auto ii = node.proc_istart_; ii < node.proc_iend_; ii++) {
          for (auto jj = node.proc_jstart_; jj < node.proc_jend_; jj++) {
            const auto i = ii - this->proc_istart_;
            const auto j = jj - this->proc_jstart_;
            ranklist[incl_nrank] = i * this->y_nnod_ + j;
            ranklist_group[incl_nrank] = (ii)*nnod.y + (jj);
            incl_nrank += 1;
          }
        }
      } else {
        for (auto jj = node.proc_jstart_; jj < node.proc_jend_; jj++) {
          for (auto ii = node.proc_istart_; ii < node.proc_iend_; ii++) {
            const auto i = ii - this->proc_istart_;
            const auto j = jj - this->proc_jstart_;
            ranklist[incl_nrank] = j * this->x_nnod_ + i;
            ranklist_group[incl_nrank] = (jj)*nnod.x + (ii);
            incl_nrank += 1;
          }
        }
      }
      return incl_nrank;
    }();

    if (nrank > 1) {
      MPI_Group_incl(FS_libs::FS_get_group(), nrank, ranklist_group.get(),
                     &node.MERGE_GROUP_);

      node.group_processranklist_.reset(new int[nrank]);
      std::iota(ranklist_group.get(), ranklist_group.get() + nrank, 0);

      MPI_Group_translate_ranks(node.MERGE_GROUP_, nrank, ranklist_group.get(),
                                FS_libs::FS_get_group(),
                                node.group_processranklist_.get());

      if (node.MERGE_GROUP_ != MPI_GROUP_NULL) {
        node.nnod_ = ni * nj;
        node.x_nnod_ = ni;
        node.y_nnod_ = nj;
        MPI_Group_rank(node.MERGE_GROUP_, &node.inod_);

        if (order == 'R') {
          const auto y_nnod = node.y_nnod_;
          node.x_inod_ = (node.inod_) / y_nnod;
          node.y_inod_ = ((node.inod_) % y_nnod);
        } else {
          const auto x_nnod = node.x_nnod_;
          node.x_inod_ = ((node.inod_) % x_nnod);
          node.y_inod_ = (node.inod_) / x_nnod;
        }

        node.FS_create_mergeXY_group();
      }
    }
  }

  for (int n = 0; n < 2; n++) {
    bt_node<Float> &node = this->sub_bt_node_[n];
    if (node.FS_node_included()) {
      node.FS_create_merge_comm_recursive();
    }
  }
}

template <class Float>
void bt_node<Float>::FS_dividing_free() {
  if (this->sub_bt_node_ != nullptr) {
    this->sub_bt_node_[0].FS_dividing_free();
    this->sub_bt_node_[1].FS_dividing_free();
  }

  if (this->MERGE_GROUP_X_ != MPI_GROUP_NULL) {
    MPI_Group_free(&this->MERGE_GROUP_X_);
  }
  if (this->MERGE_GROUP_Y_ != MPI_GROUP_NULL) {
    MPI_Group_free(&this->MERGE_GROUP_Y_);
  }
  if (this->MERGE_GROUP_ != FS_libs::FS_get_group() &&
      this->MERGE_GROUP_ != MPI_GROUP_NULL) {
    MPI_Group_free(&this->MERGE_GROUP_);
  }

  this->MERGE_GROUP_ = MPI_GROUP_NULL;
  this->MERGE_GROUP_X_ = MPI_GROUP_NULL;
  this->MERGE_GROUP_Y_ = MPI_GROUP_NULL;

  this->parent_node_ = nullptr;

  if (this->layer_ == 0) {
    if (this->procs_i_ != nullptr) {
      delete[] this->procs_i_;
    }
    if (this->procs_j_ != nullptr) {
      delete[] this->procs_j_;
    }
    this->procs_i_ = nullptr;
    this->procs_j_ = nullptr;
  }
}

template <class Float>
g1l bt_node<Float>::FS_info_G1L(char comp, int g_index) const {
  const auto NB = this->FS_get_NB();
  const auto IBLK = (g_index) / NB;

  int i_bit0 = 0;
  int i_bit1 = 0;
  if (this->div_nbit_ > 0) {
    for (auto i = this->div_nbit_; i > 0; i--) {
      const auto b = i - 1;
      if (!(this->div_bit_ & (1 << b))) {
        i_bit0 <<= 1;
        if (IBLK & (1 << b)) {
          i_bit0 |= 1;  // IBSET
        }
      } else {
        i_bit1 <<= 1;
        if (IBLK & (1 << b)) {
          i_bit1 |= 1;
        }
      }
    }
  }

  const auto rocsrc = (comp == 'R') ? i_bit0 : i_bit1;
  const auto LBLK = (comp == 'R') ? i_bit1 : i_bit0;
  // ローカルインデックス
  const auto l_index = LBLK * NB + (g_index % NB);
  return {l_index, rocsrc};
}

template <class Float>
GridIndex bt_node<Float>::FS_get_QTOP() const {
  const auto grid_info = this->FS_grid_info();
  const auto n = this->FS_get_N();
  const auto nb = this->FS_get_NB();

  const auto ii = [&]() {
    for (int i = 0; i < n; i += nb) {
      const auto row = this->FS_info_G1L('R', i);
      if (row.rocsrc == grid_info.myrow) {
        return i;
      }
    }
    return 0;
  }();

  const auto jj = [&]() {
    for (int j = 0; j < n; j += nb) {
      const auto col = this->FS_info_G1L('C', j);
      if (col.rocsrc == grid_info.mycol) {
        return j;
      }
    }
    return 0;
  }();

  const auto II = ii + this->nstart_;
  const auto JJ = jj + this->nstart_;

  const auto *root_node = this;

  while (root_node->parent_node_ != nullptr) {
    root_node = root_node->parent_node_;
  }

  return {root_node->FS_index_G2L('R', II), root_node->FS_index_G2L('C', JJ)};
}

template <class Float>
int bt_node<Float>::FS_index_L2G(char comp, int l_index, int my_roc) const {
  // 1ブロックの次数
  const auto nb = this->FS_get_NB();
  // ローカルインデックスが該当するブロック位置(0から)
  // gr_index, gc_indexはマージブロック内での相対インデックス
  const auto lblk = (l_index) / nb;

  // 行/列の場合分け
  const auto i_bit = (comp == 'R') ? std::make_pair(my_roc, lblk)
                                   : std::make_pair(lblk, my_roc);
  const auto &i_bit0 = i_bit.first;
  const auto &i_bit1 = i_bit.second;

  // 全体でのブロック位置を取得
  int iblk = 0;
  if (this->div_nbit_ > 0) {
    int ipnt0 = 0;
    int ipnt1 = 0;
    for (int i = 0; i < this->div_nbit_; i++) {
      if (!(this->div_bit_ & (1 << i))) {
        if (i_bit0 & (1 << ipnt0)) {
          iblk |= (1 << i);
        }
        ipnt0 += 1;
      } else {
        if (i_bit1 & (1 << ipnt1)) {
          iblk |= (1 << i);
        }
        ipnt1 += 1;
      }
    }
  }

  const auto g_index = iblk * nb + ((l_index) % nb);
  return g_index;
}

template <class Float>
void bt_node<Float>::print_tree() const {
  const auto nnod = FS_libs::get_procs();
  const auto inod = FS_libs::get_id();

  if (this->layer_ == 0) {
    printf("nnod, (x_nnod, y_nnod) = %d (%d, %d)\n", nnod.nod, nnod.x, nnod.y);
    printf("inod, (x_inod, y_inod) = %d (%d, %d)\n", inod.nod, inod.x, inod.y);
  }

  this->print_node();
  /*
    if(this->sub_bt_node_ != nullptr){
      for (int i = 0; i < 2; i++){
        this->sub_bt_node_[i].print_tree();
      }
    }
    */
}

template <class Float>
void bt_node<Float>::print_node() const {
  FS_libs::Nod nnod = FS_libs::get_procs();
  printf("******************************\n");
  printf("layer               = %d\n", this->layer_);
  printf("direction_horizontal= %d\n", this->direction_horizontal_);
  printf("nstart              = %d\n", this->nstart_);
  printf("nend                = %d\n", this->nend_);
  printf("nend_active         = %d\n", this->nend_active_);
  printf("proc_istart         = %d\n", this->proc_istart_);
  printf("proc_iend           = %d\n", this->proc_iend_);
  printf("proc_jstart         = %d\n", this->proc_jstart_);
  printf("proc_jend           = %d\n", this->proc_jend_);
  printf("block_start         = %d\n", this->block_start_);
  printf("block_end           = %d\n", this->block_end_);
  printf("parent_node = %p\n", this->parent_node_);
  printf("child_node  = %p\n", this->sub_bt_node_);
  printf("procs_i = ");
  for (int i = 0; i < nnod.nod; i++) {
    printf("%d ", this->procs_i_[i]);
  }
  printf("\nprocs_j = ");
  for (int j = 0; j < nnod.nod; j++) {
    printf("%d ", this->procs_j_[j]);
  }
  printf("\nmerge procs    = %d\n ", this->nnod_);
  printf("merge procs X  = %d\n ", this->x_nnod_);
  printf("merge procs Y  = %d\n ", this->y_nnod_);
  printf("merge rankid   = %d\n ", this->inod_);
  printf("merge rankid X = %d\n ", this->x_inod_);
  printf("merge rankid Y = %d\n ", this->y_inod_);
  printf("bit stream     = %d\n ", this->div_bit_);
  printf("#dights of bit = %d\n ", this->div_nbit_);
}
}  // namespace FS_dividing

#endif
