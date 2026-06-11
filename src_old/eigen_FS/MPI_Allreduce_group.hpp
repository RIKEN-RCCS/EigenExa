#pragma once
#ifndef MPI_ALLREDUCE_GROUP_HPP
#define MPI_ALLREDUCE_GROUP_HPP

#include <mpi.h>

#include <memory>
#include <cmath>
#include <memory>
#include <algorithm>
#include <numeric>

namespace eigen_FS {
namespace {
namespace MPI_Group_property {
template <typename Number>
class MPI_Group_type {
 public:
  Number *sbuf;
  Number *rbuf;
  size_t count;
  MPI_Datatype datatype;
  MPI_Op op;
  MPI_Comm comm;
  MPI_Group group;
  int group_rank;
  int group_size;
  int err;

  MPI_Group comm_group;
  int comm_group_rank;
  int comm_group_size;
  std::unique_ptr<int[]> comm_group_ranklist;
};
}  // namespace MPI_Group_property
namespace MPI_Allreduce_main {
using MPI_Group_property::MPI_Group_type;
using std::abs;
template <typename T>
const MPI_Datatype MPI_TYPE = 0;
template <>
const MPI_Datatype MPI_TYPE<int> = MPI_INT;
template <>
const MPI_Datatype MPI_TYPE<double> = MPI_DOUBLE;
template <>
const MPI_Datatype MPI_TYPE<float> = MPI_FLOAT;

template <typename Number>
void comm_op(size_t &n, size_t &head, Number sbuf[], Number rbuf[], MPI_Comm comm,
             MPI_Op op, int myrank, int pair) {
  const auto tag = 1;
  size_t head_r, head_s, count_r, count_s;
  if (myrank < pair) {
    head_r = 0;
    head_s = n / 2;
    count_r = n / 2;
    count_s = n - n / 2;
  } else {
    head_r = n / 2;
    head_s = 0;
    count_r = n - n / 2;
    count_s = n / 2;
  }
  MPI_Request req_s, req_r;
  MPI_Isend(&sbuf[head_s], count_s, MPI_TYPE<Number>, pair, tag, comm, &req_s);
  MPI_Irecv(rbuf, count_r, MPI_TYPE<Number>, pair, tag, comm, &req_r);

  MPI_Wait(&req_s, MPI_STATUS_IGNORE);
  MPI_Wait(&req_r, MPI_STATUS_IGNORE);
  if (op == MPI_SUM) {
    for (size_t i = 0; i < count_r; i++) {
      sbuf[head_r + i] += rbuf[i];
    }
  } else if (op == MPI_PROD) {
    for (size_t i = 0; i < count_r; i++) {
      sbuf[head_r + i] *= rbuf[i];
    }
  }

  n = count_r;
  head += head_r;
}
template <typename Number>
void comm_op_rev(size_t n, Number sbuf[], Number rbuf[], MPI_Comm comm, int myrank,
                 int pair) {
  auto tag = 1;
  size_t head_r, head_s, count_r, count_s;
  if (myrank < pair) {
    head_r = n / 2;
    head_s = 0;
    count_r = n - n / 2;
    count_s = n / 2;
  } else {
    head_r = 0;
    head_s = n / 2;
    count_r = n / 2;
    count_s = n - n / 2;
  }
  MPI_Request req_s, req_r;
  MPI_Isend(&sbuf[head_s], count_s, MPI_TYPE<Number>, pair, tag, comm, &req_s);
  MPI_Irecv(rbuf, count_r, MPI_TYPE<Number>, pair, tag, comm, &req_r);

  MPI_Wait(&req_s, MPI_STATUS_IGNORE);
  MPI_Wait(&req_r, MPI_STATUS_IGNORE);
  std::copy_n(rbuf, count_r, &sbuf[head_r]);
}

size_t get_level_size(size_t group_size) {
  size_t level = 1;
  size_t i = group_size;
  while (i != 1) {
    i /= 2;
    level += 1;
    if (i / 2 == 0) {
      return level;
    }
  }
  return level;
}

template <typename Number>
void Group_Allreduce(MPI_Group_type<Number> &mygroup) {
  const auto myrank = mygroup.group_rank;
  auto count = mygroup.count;

  auto level_size = get_level_size(mygroup.group_size);
  std::unique_ptr<ssize_t[]> count_level(new ssize_t[level_size]);
  std::unique_ptr<ssize_t[]> head_level(new ssize_t[level_size]);
  std::unique_ptr<ssize_t[]> pair_level(new ssize_t[level_size]);
  std::fill_n(count_level.get(), level_size, -1);
  std::fill_n(head_level.get(), level_size, -1);
  std::fill_n(pair_level.get(), level_size, -1);

  // 2べきのプロセスでのみ動作
  auto i = mygroup.group_size;
  auto step = mygroup.group_size / 2;
  if (step <= myrank) {
    step = -step;
  }
  int step0 = 0;
  size_t  level = 0;
  size_t head = 0;

  while (i != 1) {
    i /= 2;
    const auto pair = myrank + step;
    count_level[level] = count;
    head_level[level] = head;
    pair_level[level] = pair;

    comm_op<Number>(count, head, &mygroup.sbuf[head], mygroup.rbuf,
                    mygroup.comm, mygroup.op,
                    mygroup.comm_group_ranklist[myrank],
                    mygroup.comm_group_ranklist[pair]);
    level += 1;
    if (i / 2 == 0) {
      break;
    }
    if (step < 0) {
      step0 += i;
    }
    step = abs(step / 2);
    if ((step + step0) <= myrank) {
      step = -step;
    }
  }
  // 逆回転
  i = 2;
  while (i <= mygroup.group_size) {
    level -= 1;
    auto count = count_level[level];
    auto head = head_level[level];
    auto pair = pair_level[level];
    comm_op_rev<Number>(count, &mygroup.sbuf[head], mygroup.rbuf, mygroup.comm,
                         mygroup.comm_group_ranklist[myrank],
                        mygroup.comm_group_ranklist[pair]);
    i *= 2;
  }
  std::copy_n(mygroup.sbuf, mygroup.count, mygroup.rbuf);
}
}  // namespace MPI_Allreduce_main
namespace MPI_Allreduce_group {
using MPI_Allreduce_main::Group_Allreduce;
using MPI_Group_property::MPI_Group_type;
template <typename Number>
void set_group2comm_ranklist(MPI_Group_type<Number> &mygroup) {
  std::unique_ptr<int[]> group_ranks(new int[mygroup.group_size]);
  mygroup.err = MPI_Comm_group(mygroup.comm, &mygroup.comm_group);
  mygroup.err = MPI_Group_size(mygroup.comm_group, &mygroup.comm_group_size);
  mygroup.err = MPI_Group_rank(mygroup.comm_group, &mygroup.comm_group_rank);
  mygroup.comm_group_ranklist.reset( new int[mygroup.group_size]);

  std::iota(group_ranks.get(), group_ranks.get() + mygroup.group_size, 0);
  std::fill_n(mygroup.comm_group_ranklist.get(), mygroup.group_size, -1);

  mygroup.err = MPI_Group_translate_ranks(mygroup.group, mygroup.group_size,
                                          group_ranks.get(), mygroup.comm_group,
                                          mygroup.comm_group_ranklist.get());
}
template <typename Number>
void free_group(MPI_Group_type<Number> &mygroup, Number rbuf[], int count) {
  mygroup.err = MPI_Group_free(&mygroup.comm_group);
  if (count < mygroup.group_size) {
    std::copy_n(mygroup.rbuf, count, rbuf);
    delete[] mygroup.sbuf;
    delete[] mygroup.rbuf;
  }
}
template <typename Number>
void set_group(Number sbuf[], Number rbuf[], size_t count, MPI_Datatype datatype,
               MPI_Op op, MPI_Comm comm, MPI_Group group,
               MPI_Group_type<Number> &mygroup) {
  mygroup.count = count;
  mygroup.datatype = datatype;
  mygroup.op = op;
  mygroup.comm = comm;
  mygroup.group = group;
  mygroup.err = 0;
  MPI_Group_size(mygroup.group, &mygroup.group_size);
  MPI_Group_rank(mygroup.group, &mygroup.group_rank);
  set_group2comm_ranklist(mygroup);
  if (mygroup.count < static_cast<size_t>(mygroup.group_size)) {
    mygroup.sbuf = new Number[mygroup.group_size];
    mygroup.rbuf = new Number[mygroup.group_size];

    std::fill_n(mygroup.rbuf, mygroup.group_size, 0);
    std::copy_n(sbuf, mygroup.count, mygroup.sbuf);
    std::fill(&mygroup.sbuf[mygroup.count], &mygroup.sbuf[mygroup.group_size],
              0);
    mygroup.count = mygroup.group_size;
  } else {
    mygroup.sbuf = sbuf;
    mygroup.rbuf = rbuf;
  }
}
}  // namespace MPI_Allreduce_group
}  // namespace
using MPI_Allreduce_group::free_group;
using MPI_Allreduce_group::set_group;
using MPI_Allreduce_main::Group_Allreduce;
using MPI_Group_property::MPI_Group_type;
template <typename Number>
int MPI_Group_Allreduce(Number sbuf[], Number rbuf[], size_t count,
                        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                        MPI_Group group) {
  MPI_Group_type<Number> mygroup;
  set_group<Number>(sbuf, rbuf, count, datatype, op, comm, group, mygroup);
  Group_Allreduce<Number>(mygroup);
  free_group<Number>(mygroup, rbuf, count);
  return 0;
}
}  // namespace eigen_FS
#endif
