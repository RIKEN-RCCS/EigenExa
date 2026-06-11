#pragma once
#ifndef EIGEN_DEVEL_HPP
#define EIGEN_DEVEL_HPP

namespace eigen_devel {

extern "C" void eigen_abort();

extern "C" void eigen_timer_reset(int bcast, int reduce, int redist, int gather);

extern "C" double eigen_timer_print();




}  // namespace eigen_devel

#endif
