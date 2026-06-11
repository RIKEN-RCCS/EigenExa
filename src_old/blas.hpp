#pragma once
#ifndef FS_BLAS_HPP
#define FS_BLAS_HPP

#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
#include <mkl.h>
#include <mkl_cblas.h>
#else
#include <cblas.h>
#if defined(__FUJITSU)
#include <plasma.h>
#include <core_dblas.h>
#include <core_sblas.h>
#include <lapacke.h>
#endif


extern "C" {
double dlaed4_(int *n, int *i, double d[], double z[], double delta[],
              double *rho, double *dlam, int *info);
float slaed4_(int *n, int *i, float d[], float z[], float delta[], float *rho,
             float *dlam, int *info);

double dlanst_(char *norn, int *n, const double *D, const double *E);
float slanst_(char *norn, int *n, const float *D, const float *E);

double dlamc3_(double *x, double *y);
float slamc3_(float *x, float *y);
}

#endif

namespace lapacke {
const CBLAS_LAYOUT FS_layout = CBLAS_LAYOUT::CblasColMajor;
const int FS_LAPACKE_LAYOUT = LAPACK_COL_MAJOR;

template <class Float>
inline int stedc(char compz, int n, Float *d, Float *e, Float *z, int ldz,
                 Float *work, int lwork, int *iwork, int liwork);

template <>
inline int stedc<double>(char compz, int n, double *d, double *e, double *z,
                         int ldz, double *work, int lwork, int *iwork,
                         int liwork) {
  return LAPACKE_dstedc_work(FS_LAPACKE_LAYOUT, compz, n, d, e, z, ldz, work, lwork,
                             iwork, liwork);
}

template <>
inline int stedc<float>(char compz, int n, float *d, float *e, float *z,
                        int ldz, float *work, int lwork, int *iwork,
                        int liwork) {
  return LAPACKE_sstedc_work(FS_LAPACKE_LAYOUT, compz, n, d, e, z, ldz, work, lwork,
                             iwork, liwork);
}

template <class Float>
inline int lascl(char type, int kl, int ku, Float cfrom, Float cto, int m,
                 int n, Float *a, int lda);

template <>
inline int lascl<double>(char type, int kl, int ku, double cfrom, double cto,
                         int m, int n, double *a, int lda) {
#if !defined(__FUJITSU) || defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
  return LAPACKE_dlascl_work(FS_LAPACKE_LAYOUT, type, kl, ku, cfrom, cto, m, n, a, lda);
#else
  if (type == 'G'){
    return CORE_dlascl(PlasmaGeneral , kl, ku, cfrom, cto, m, n, a, lda);
  }
  return -1;
#endif
}

template <>
inline int lascl<float>(char type, int kl, int ku, float cfrom, float cto,
                        int m, int n, float *a, int lda) {
#if !defined(__FUJITSU) || defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
  return LAPACKE_slascl_work(FS_LAPACKE_LAYOUT, type, kl, ku, cfrom, cto, m, n, a, lda);
#else
  if (type == 'G'){
    return CORE_slascl(PlasmaGeneral , kl, ku, cfrom, cto, m, n, a, lda);
  }
  return -1;
#endif
}

template <class Integer, class Float>
inline Float lanst(char norm, Integer n, const Float *D, const Float *E);

template <>
inline double lanst<int, double>(char norm, int n, const double *D,
                                 const double *E) {
  return dlanst_(&norm, &n, D, E);
}

template <>
inline float lanst<int, float>(char norm, int n, const float *D,
                               const float *E) {
  return slanst_(&norm, &n, D, E);
}

template <class Float>
inline Float lapy2(Float x, Float y);

template <>
inline double lapy2<double>(double x, double y) {
  return LAPACKE_dlapy2(x, y);
}

template <>
inline float lapy2<float>(float x, float y) {
  return LAPACKE_slapy2(x, y);
}

template <class Float>
inline void gemm(const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB,
                 const int M, const int N, const int K, const Float alpha,
                 const Float *A, const int lda, const Float *B, const int ldb,
                 const Float beta, Float *C, const int ldc);

template <>
inline void gemm<double>(const CBLAS_TRANSPOSE TransA,
                         const CBLAS_TRANSPOSE TransB, const int M, const int N,
                         const int K, const double alpha, const double *A,
                         const int lda, const double *B, const int ldb,
                         const double beta, double *C, const int ldc) {
  cblas_dgemm(FS_layout, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta,
              C, ldc);
}

template <>
inline void gemm<float>(const CBLAS_TRANSPOSE TransA,
                        const CBLAS_TRANSPOSE TransB, const int M, const int N,
                        const int K, const float alpha, const float *A,
                        const int lda, const float *B, const int ldb,
                        const float beta, float *C, const int ldc) {
  cblas_sgemm(FS_layout, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta,
              C, ldc);
}

template <class Float>
inline void copy(int n, const Float *X, int incX, Float *Y, int incY);

template <>
inline void copy<double>(int n, const double *X, int incX, double *Y,
                         int incY) {
  cblas_dcopy(n, X, incX, Y, incY);
}

template <>
inline void copy<float>(int n, const float *X, int incX, float *Y, int incY) {
  cblas_scopy(n, X, incX, Y, incY);
}


template <class Integer, class Float>
inline int iamax(Integer n, Float dx[], Integer incx);

template <>
inline int iamax<int, double>(int n, double dx[], int incx) {
  return cblas_idamax(n, dx, incx);
}

template <>
inline int iamax<int, float>(int n, float dx[], int incx) {
  return cblas_isamax(n, dx, incx);
}

template <class Float>
inline void scal(int n, Float alpha, Float x[], int incx);

template <>
inline void scal<double>(int n, double alpha, double x[], int incx) {
  cblas_dscal(n, alpha, x, incx);
}

template <>
inline void scal<float>(int n, float alpha, float x[], int incx) {
  cblas_sscal(n, alpha, x, incx);
}

template <class Float>
inline Float nrm2(int n, Float X[], int incX);

template <>
inline double nrm2(int n, double X[], int incX) {
  return cblas_dnrm2(n, X, incX);
}

template <>
inline float nrm2(int n, float X[], int incX) {
  return cblas_snrm2(n, X, incX);
}

template <class Float>
inline int laed4(int n, int i, Float d[], Float z[], Float delta[], Float rho,
                 Float &dlam);

template <>
inline int laed4(int n, int i, double d[], double z[], double delta[],
                 double rho, double &dlam) {
  int info;
  dlaed4_(&n, &i, d, z, delta, &rho, &dlam, &info);
  return info;
}

template <>
inline int laed4(int n, int i, float d[], float z[], float delta[], float rho,
                 float &dlam) {
  int info;
  slaed4_(&n, &i, d, z, delta, &rho, &dlam, &info);
  return info;
}

template <class Float>
inline Float lamc3(Float x, Float y);

template <>
inline double lamc3(double x, double y) {
  return dlamc3_(&x, &y);
}
template <>
inline float lamc3(float x, float y) {
  return slamc3_(&x, &y);
}

template <class Float>
inline void rot(int N, Float *X, int incX, Float *Y, int incY, Float c,
                Float s);

template <>
inline void rot(int N, double *X, int incX, double *Y, int incY, double c,
                double s) {
  cblas_drot(N, X, incX, Y, incY, c, s);
}

template <>
inline void rot(int N, float *X, int incX, float *Y, int incY, float c,
                float s) {
  cblas_srot(N, X, incX, Y, incY, c, s);
}
}  // namespace lapacke
#endif
