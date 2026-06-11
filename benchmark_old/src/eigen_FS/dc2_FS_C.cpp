#include <mpi.h>

#include "dc2_FS.hpp"
extern "C" {
void dc2_FS_f64(int n, int nvec, double d[], double e[], double z[], int ldz,
                int *info, double *ret) {
  eigen_FS::dc2_FS<double>(n, nvec, d, e, z, ldz, info, ret);
}

void dc2_FS_f32(int n, int nvec, float d[], float e[], float z[], int ldz,
                int *info, float *ret) {
  eigen_FS::dc2_FS<float>(n, nvec, d, e, z, ldz, info, ret);
}
}