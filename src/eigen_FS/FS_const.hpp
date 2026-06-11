#pragma once
#ifndef FS_CONST_HPP
#define FS_CONST_HPP

#include <mpi.h>

namespace FS_const {
template <typename T>
const MPI_Datatype MPI_TYPE = 0;
template <>
const MPI_Datatype MPI_TYPE<int> = MPI_INT;
template <>
const MPI_Datatype MPI_TYPE<double> = MPI_DOUBLE;
template <>
const MPI_Datatype MPI_TYPE<float> = MPI_FLOAT;
template <class Float>
constexpr Float ZERO = 0.0e0;
template <class Float>
constexpr Float HALF = 0.5e0;
template <class Float>
constexpr Float ONE = 1.0e0;
template <class Float>
constexpr Float TWO = 2.0e0;
template <class Float>
constexpr Float THREE = 3.0e0;
template <class Float>
constexpr Float FOUR = 4.0e0;
template <class Float>
constexpr Float FIVE = 5.0e0;
template <class Float>
constexpr Float SIX = 6.0e0;
template <class Float>
constexpr Float SEVEN = 7.0e0;
template <class Float>
constexpr Float EIGHT = 8.0e0;
template <class Float>
constexpr Float NINE = 9.0e0;
template <class Float>
constexpr Float TEN = 1.0e1;

template <class Float>
constexpr Float MHALF = -0.5e0;
template <class Float>
constexpr Float MONE = -1.0e0;
template <class Float>
constexpr Float MTWO = -2.0e0;
}  // namespace FS_const

#endif
