#pragma once


namespace eigen_libs {
struct comm{
  int comm;
  int x;
  int y;
};

struct id {
  int id;
  int x;
  int y;
};

struct procs {
  int procs;
  int x;
  int y;
};
extern "C" {
void eigen_get_procs(int* nnod, int* x_nnod, int* y_nnod);
void eigen_get_id(int* inod, int* x_inod, int* y_inod);
void eigen_get_comm(int* eigen_comm, int* eigen_x_comm, int* eigen_y_comm);
void eigen_get_grid_major(char *major);
int eigen_owner_node(int ictr, int nnod, int inod);
}

inline procs get_procs() { 
  eigen_libs::procs procs = {};
  eigen_get_procs(&procs.procs, &procs.x, &procs.y);
  return procs;
}

inline id get_id() { 
  eigen_libs::id id = {};
  eigen_get_id(&id.id, &id.x, &id.y);
  return id;
}

inline comm get_comm() { 
  eigen_libs::comm comm = {};
  eigen_get_comm(&comm.comm, &comm.x, &comm.y);
  return comm;
}

inline char get_grid_major() { 
  char major;
  eigen_get_grid_major(&major);
  return major;
}
}

