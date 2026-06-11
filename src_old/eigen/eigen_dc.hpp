#pragma once
#ifndef EIGEN_DC_HPP
#define EIGEN_DC_HPP
namespace eigen_dc {
extern "C" double flops, dgemm_time, dgemm_dummy[2];
extern "C" double p_times, p_timez, p_timer;
extern "C" double p_time0, p_time2, p_time3;
}  // namespace eigen_dc
#endif
