#include <optional>
#include "eigen_libs0.hpp"
#include "eigen_devel.hpp"

std::vector<int> eigen_devel::p0_;
std::vector<int> eigen_devel::q0_;

int64_t version() {
  eigen_libs0::Eigen_Version;
}
namespace eigen_libs0 {
extern "C" void eigen_get_matdims0_C(
    int64_t n,
    int64_t* nx,
    int64_t* ny,
    const int64_t* m_forward,
    const int64_t* m_backward,
    const char* mode)
{
  int64_t mode_len = 1;
  std::optional<int64_t> mf = m_forward ? std::optional<int64_t>(*m_forward) : std::nullopt;
  std::optional<int64_t> mb = m_backward ? std::optional<int64_t>(*m_backward) : std::nullopt;
  std::optional<char> md = mode ? std::optional<char>(*mode) : std::nullopt;
}


extern "C" int64_t eigen_loop_start_C(int64_t istart, int64_t nnod, int64_t inod) {
  return eigen_loop_start(istart, nnod, inod);
}

extern "C" int64_t eigen_loop_start_XY_C(int64_t istart, const char* pdir_cstr, int64_t inod_val) {
  std::optional<int> inod = inod_val ? std::optional<int>(inod_val) : std::nullopt;
  std::string pdir(pdir_cstr);
  return eigen_loop_start(istart, pdir, inod);
}

extern "C" int64_t eigen_loop_end_C(int64_t iend, int64_t nnod, int64_t inod) {
  return eigen_loop_end(iend, nnod, inod);
}

extern "C" int64_t eigen_loop_end_XY_C(int64_t iend, const char* pdir_cstr, int64_t inod_val) {
  std::string pdir(pdir_cstr);
  std::optional<int> inod = inod_val ? std::optional<int>(inod_val) : std::nullopt;
  return eigen_loop_end(iend, pdir, inod);
}

extern "C" int64_t eigen_translate_l2g_C(int64_t ictr, int64_t nnod, int64_t inod) {
  return eigen_translate_l2g(ictr, nnod, inod);
}

extern "C" int64_t eigen_translate_l2g_XY_C(int64_t ictr, const char* pdir_cstr, int64_t inod_val) {
  std::string pdir(pdir_cstr);
  std::optional<int> inod = inod_val ? std::optional<int>(inod_val) : std::nullopt;
  return eigen_translate_l2g(ictr, pdir, inod);
}

extern "C" int64_t eigen_translate_g2l_C(int64_t ictr, int64_t nnod, int64_t inod) {
  return eigen_translate_g2l(ictr, nnod, inod);
}

extern "C" int64_t eigen_translate_g2l_XY_C(int64_t ictr, const char* pdir_cstr, int64_t inod_val) {
  std::string pdir(pdir_cstr);
  std::optional<int> inod = inod_val ? std::optional<int>(inod_val) : std::nullopt;
  return eigen_translate_g2l(ictr, pdir, inod);
}

extern "C" int64_t eigen_owner_node_C(int64_t ictr, int64_t nnod, int64_t inod) {
  int64_t ret = eigen_owner_node(ictr, nnod, inod);
  return ret;
}

extern "C" int64_t eigen_owner_node_XY_C(int64_t ictr, const char* pdir_cstr, int64_t inod_val) {
  std::string pdir(pdir_cstr);
  std::optional<int> inod = inod_val ? std::optional<int>(inod_val) : std::nullopt;
  int64_t ret = eigen_owner_node(ictr, pdir, inod);
  return ret;
}

}
