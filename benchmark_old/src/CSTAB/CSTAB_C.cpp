#include <complex>
#include "CSTAB.hpp"

namespace CSTAB {
  extern "C" void CSTAB_get_optdim_c(int n_min, int n_unroll, int delta_L1, int delta_L2, int* n_opt) {
    int64_t n_opt_64;
    get_optdim(static_cast<int64_t>(n_min), static_cast<int64_t>(n_unroll), static_cast<int64_t>(delta_L1), static_cast<int64_t>(delta_L2), n_opt_64);
    *n_opt = static_cast<int>(n_opt_64);
  }

  extern "C" void CSTAB_adjust_base_real(double* a, double* b, int64_t* offset) {
    adjust_base<double>(a, b, *offset);
  }

  extern "C" void CSTAB_adjust_page_real(double* a, double* b, int64_t* offset) {
    adjust_page<double>(a, b, *offset);
  }

  extern "C" void CSTAB_adjust_base_hermite(std::complex<double>* a, std::complex<double>* b, int64_t* offset) {
    adjust_base<std::complex<double>>(a, b, *offset);
  }

  extern "C" void CSTAB_adjust_page_hermite(std::complex<double>* a, std::complex<double>* b, int64_t* offset) {
    adjust_page<std::complex<double>>(a, b, *offset); 
  }

  extern "C" void CSTAB_round_offset(int* offset) {
    int64_t off64 = static_cast<int64_t>(*offset);
    round_offset(&off64);
    *offset = static_cast<int>(off64);
  }

  extern "C" void CSTAB_round_offset_h(int* offset) {
    int64_t off64 = static_cast<int64_t>(*offset);
    round_offset_h(&off64);
    *offset = static_cast<int>(off64);
  }

}
