#include <memory>
#include "../eigen/eigen_trbak.hpp"
#include "trbakwy4.hpp"

namespace trbakwy4 {
  extern "C" void trbk_decide_overlap_level_Cpp(int i) {
    trbk_decide_overlap_level(static_cast<int64_t>(i));
  }

  extern "C" void eigen_common_trbakwy_dbl(
    int n,
    int nvec,
    const double* a, int nma0,
    double* z, int nmz0,
    double* beta,
    int m0,
    int iblk
  ) {
      eigen_common_trbakwy<double>(
	      n, nvec,
              a, nma0,
              z, nmz0,
              beta,
              m0,
              iblk);
  }

  extern "C" void eigen_trbakwy_body_dbl(
    int n, int nvec,
    const double* a, int nma,
    double* z, int nmz,
    double* beta,
    double* v1, double* v2, double* v3, int nm,
    int m,
    double* ss1, double* ss2, double* ss3,
    double* tt,
    int iblk,
    int nss
  ) {
	  eigen_trbakwy_body<double>(n, nvec,
			  a, nma,
			  z,nmz,
			  beta,
			  v1, v2, v3, nm,
			  m,
			  ss1, ss2, ss3,
			  tt,
			  iblk,
			  nss);
  }

  extern "C" void eigen_trbakwy_block_body_dbl(
    int local_nvec,
    double* z, int nmz,
    double* v, int nm, int m, int i,
    double* ss, double* tt, int nss, int iblk,
    double* dcom, double* dx, double* dy, double* dz) {
      eigen_trbakwy_block_body<double>(
      local_nvec, 
      z, nmz, 
      v, nm, m, i, 
      ss, tt, nss, iblk, 
      *dcom, *dx, *dy, *dz);
  }

  extern "C" void eigen_trbakwy_block_body1_dbl(
    double* z, int nmz,
    double* v, int nm, int m,
    double* ss, double* sm,
    int i_2, int i_3, int j_2, int j_3) {
      eigen_trbakwy_block_body1<double>(
      z, nmz, v, nm, m, ss, sm, i_2, i_3, j_2, j_3);
  }

  extern "C" void eigen_trbakwy_block_body2_dbl(
    double* z, int nmz,
    double* v, int nm, int m,
    double* ss, double* sm,
    int i_2, int i_3, int j_2, int j_3) {
      eigen_trbakwy_block_body2<double>(
      z, nmz, v, nm, m, ss, sm, i_2, i_3, j_2, j_3);
  }

  extern "C" void trbakwy_datacast_dbl(
    const int iloop_end, const int m, const int i,
    const double* a, const int nma,
    double* v, const int nm,
    double* ss,
#ifndef USE_BCASTW
    double* wk,
#endif
    const int iblk) {
      trbakwy_datacast<double>(iloop_end, m, i, a, nma, v, nm, ss,
#ifndef USE_BCASTW
      wk,
#endif
      iblk);
  }

  extern "C" void trbakwy_alloc_buffer(int size) {
      eigen_trbak::trbk_buf = std::make_unique<double[]>(size);
  }

}

