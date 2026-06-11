#include <mpi.h>

#include "FS_libs.hpp"

namespace FS_libs {
MPI_Comm FS_COMM_WORLD = MPI_COMM_WORLD;
int FS_MYRANK = 0;
bool FS_COMM_MEMBER = false;
MPI_Group FS_GROUP = MPI_GROUP_NULL;

process_grid FS_node = {};
char FS_GRID_major = 'C';
extern "C" void FS_init_c(int comm, char order) {
  FS_init(MPI_Comm_f2c(comm), order);
}
extern "C" void FS_free_c() { FS_free(); }

extern "C" void FS_get_matdims_c(int n, int &nx, int &ny){
  FS_get_matdims(n, nx, ny);
}
extern "C" int FS_get_myrank_c() { return FS_MYRANK; }
}  // namespace FS_libs
