#pragma once
#ifndef CSTAB_HPP
#define CSTAB_HPP

#include <cmath>
#include<stdio.h>

extern "C" {
  void get_delta_(void* a1, void* a2, int64_t* i);
}

namespace CSTAB {
  // Power 4
  // constexpr int64_t L1_SIZE   = 32 * 1024;
  // constexpr int64_t L1_WAY    = 2;
  // constexpr int64_t L2_SIZE   = 512 * 1024;
  // constexpr int64_t L2_WAY    = 8;
  
  // Pentium 4 (Northwood or prior)
  // constexpr int64_t L1_SIZE   = 8 * 1024;
  // constexpr int64_t L1_WAY    = 4;
  // constexpr int64_t L2_SIZE   = 512 * 1024;
  // constexpr int64_t L2_WAY    = 8;
  
  // Pentium 4 (Prescott)
  // constexpr int64_t L1_SIZE   = 16 * 1024;
  // constexpr int64_t L1_WAY    = 8;
  // constexpr int64_t L2_SIZE   = 1024 * 1024;
  // constexpr int64_t L2_WAY    = 8;
  
  // Celeron D (Prescott)
  // constexpr int64_t L1_SIZE   = 16 * 1024;
  // constexpr int64_t L1_WAY    = 8;
  // constexpr int64_t L2_SIZE   = 256 * 1024;
  // constexpr int64_t L2_WAY    = 4;
  
  // Itanium 2
  // constexpr int64_t L1_SIZE   = 16 * 1024;
  // constexpr int64_t L1_WAY    = 8;
  // constexpr int64_t L2_SIZE   = 256 * 1024;
  // constexpr int64_t L2_WAY    = 8;
  
  // Core2
  // constexpr int64_t L1_SIZE   = 32 * 1024;
  // constexpr int64_t L1_WAY    = 8;
  // constexpr int64_t L2_SIZE   = 2 * 1024 * 1024;
  // constexpr int64_t L2_WAY    = 16;
  
  // A64FX (active configuration)
  constexpr int64_t L1_SIZE       = 64 * 1024;
  constexpr int64_t L1_WAY        = 4;
  constexpr int64_t L1_LINE_SIZE  = 256;
  constexpr int64_t L2_SIZE       = 8 * 1024 * 1024;
  constexpr int64_t L2_WAY        = 16;
  constexpr int64_t L2_LINE_SIZE  = 256;
  
  // Derived cache parameters
  constexpr int64_t L1_LSIZE   = (L1_SIZE / L1_WAY) / 8;
  constexpr int64_t L1_WINDOW  = L1_LINE_SIZE / 8;
  constexpr int64_t L2_LSIZE   = (L2_SIZE / L2_WAY) / 8;
  constexpr int64_t L2_WINDOW  = 2 * L2_LINE_SIZE / 8;
  
  constexpr int64_t n_columns      = L2_LSIZE;
  constexpr int64_t PREFETCH_SIZE  = 512 / 8;
  
  // For IA32-Linux
  constexpr int64_t PAGE_SIZE   = 4096;
  constexpr int64_t PAGE_LSIZE  = PAGE_SIZE / 8;
  
  inline void get_optdim(int64_t n_min, int64_t n_unroll, int64_t delta_L1, int64_t delta_L2, int64_t& n_opt) {
    constexpr double ONE = 1.0;
    int64_t n_delta = 0;
    n_opt = n_min;
    while (true) {
      n_opt = ((n_opt - 1) / L1_WINDOW + 1);
      n_opt = (n_opt / 2) * 2 + 1;
      n_opt = n_opt * L1_WINDOW;

      n_delta = 0;

      for (int64_t i = 1; i <= ((n_unroll * 1.2 - ONE) / L1_WAY + 1); ++i) {
        int64_t k = (i * n_opt + L1_LSIZE / 2) % L1_LSIZE - L1_LSIZE / 2;
        if (std::abs(k) <= delta_L1 / 2) {
          n_delta = (delta_L1 / 2 - k - 1) / i + 1;
          break;
        }
      }

      if (n_delta == 0) {
        for (int64_t i = 1; i <= ((n_unroll * 1.2 - ONE) / L2_WAY + 1); ++i) {
          int64_t k = (i * n_opt + L2_LSIZE / 2) % L2_LSIZE - L2_LSIZE / 2;
          if (std::abs(k) <= delta_L2 / 2) {
            n_delta = (delta_L2 / 2 - k - 1) / i + 1;
            break;
          }
        }
      }
      if (n_delta == 0) break;
      n_opt += n_delta;
    }
  }

  template <typename T>
  inline void adjust_base(T* a, T* b, int64_t& offset) {
    get_delta_((void*)a, (void*)b, &offset);
    constexpr int64_t elem_size = sizeof(T);
    offset /= elem_size;
    int64_t L2_LSIZE_adj = L2_LSIZE / (elem_size / 8);
    if (offset > 0)
      offset = (L2_LSIZE_adj - (offset % L2_LSIZE_adj)) % L2_LSIZE_adj;
    else
      offset = (L2_LSIZE_adj + (-offset % L2_LSIZE_adj)) % L2_LSIZE_adj;
  }

  template <typename T>
  inline void adjust_page(T* a, T* b, int64_t& offset) {
    get_delta_((void*)a, (void*)b, &offset);
    constexpr int64_t elem_size = sizeof(T);
    offset /= elem_size;
    int64_t PAGE_LSIZE_adj = PAGE_LSIZE / (elem_size / 8);
    if (offset > 0)
      offset = (PAGE_LSIZE_adj - (offset % PAGE_LSIZE_adj)) % PAGE_LSIZE_adj;
    else
      offset = (PAGE_LSIZE_adj + (-offset % PAGE_LSIZE_adj)) % PAGE_LSIZE_adj;
  }

  // Round up the offset by L2 size (real)
  inline void round_offset(int64_t* offset) {
    *offset = *offset%L2_LSIZE;
  }

  // Round up the offset by L2 size (complex)
  inline void round_offset_h(int64_t* offset) {
    int64_t L2_LSIZE_h = L2_LSIZE / 2;
    *offset = *offset%L2_LSIZE_h;
  }

}
#endif
