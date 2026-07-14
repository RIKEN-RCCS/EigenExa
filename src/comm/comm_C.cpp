
#include "comm.hpp"
#include <mpi.h>

namespace comm {

  extern "C" void comm_barrier(MPI_Fint* fcomm) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    barrier(comm);
  }

  extern "C" void comm_pack_dbl(double* buf, int* n, double* buffer, int* ptr) {
    int64_t p = static_cast<int64_t>(*ptr);
    pack<double>(buf, static_cast<int64_t>(*n), buffer, p);
    *ptr = static_cast<int>(p);
  }

  extern "C" void comm_pack1_dbl(double* buf, double* buffer, int* ptr) {
    int64_t p = static_cast<int64_t>(*ptr);
    pack1<double>(*buf, buffer, p);
    *ptr = static_cast<int>(p);
  }

  extern "C" void comm_unpack_dbl(double* buf, int* n, double* buffer, int* ptr) {
    int64_t p = static_cast<int64_t>(*ptr);
    unpack<double>(buf, static_cast<int64_t>(*n), buffer, p);
    *ptr = static_cast<int>(p);
  }

  extern "C" void comm_unpack1_dbl(double* buf, double* buffer, int* ptr) {
    int64_t p = static_cast<int64_t>(*ptr);
    unpack1<double>(*buf, buffer, p);
    *ptr = static_cast<int>(p);
  }

  extern "C" void comm_send_dbl(double* buf, int* n, int* idest, MPI_Fint* fcomm) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    send<double>(buf, static_cast<int64_t>(*n), static_cast<int64_t>(*idest), comm);
  }

  extern "C" void comm_send_dblt(double* buf, int* n, int* idest, int* itag, MPI_Fint* fcomm) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    send_tagged<double>(buf, static_cast<int64_t>(*n), static_cast<int64_t>(*idest), static_cast<int64_t>(*itag), comm);
  }

  extern "C" void comm_isend_dbl(double* buf, int* n, int* idest, MPI_Fint* ireq, MPI_Fint* fcomm) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    MPI_Request req;
    isend<double>(buf, static_cast<int64_t>(*n), static_cast<int64_t>(*idest), &req, comm);
    *ireq = MPI_Request_c2f(req);
  }

  extern "C" void comm_isend_dblt(double* buf, int* n, int* idest, int* itag, MPI_Fint* ireq, MPI_Fint* fcomm) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    MPI_Request req;
    isend_tagged<double>(buf, static_cast<int64_t>(*n), static_cast<int64_t>(*idest), static_cast<int64_t>(*itag), &req, comm);
    *ireq = MPI_Request_c2f(req);
  }

  extern "C" void comm_recv_dbl(double* buf, int* n, int* isrc, MPI_Fint* fcomm) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    recv<double>(buf, static_cast<int64_t>(*n), static_cast<int64_t>(*isrc), comm);
  }

  extern "C" void comm_irecv_dbl(double* buf, int* n, int* isrc, MPI_Fint* ireq, MPI_Fint* fcomm) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    MPI_Request req;
    irecv<double>(buf, static_cast<int64_t>(*n), static_cast<int64_t>(*isrc), &req, comm);
    *ireq = MPI_Request_c2f(req);
  }

  extern "C" void comm_irecv_dblt(double* buf, int* n, int* isrc, int* itag, MPI_Fint* ireq, MPI_Fint* fcomm) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    MPI_Request req;
    irecv_tagged<double>(buf, static_cast<int64_t>(*n), static_cast<int64_t>(*isrc), static_cast<int64_t>(*itag), &req, comm);
    *ireq = MPI_Request_c2f(req);
  }

  extern "C" void comm_wait_dbl(MPI_Fint* ireq) {
    MPI_Request req = MPI_Request_f2c(*ireq);
    wait(&req);
    *ireq = MPI_Request_c2f(req);
  }

  extern "C" void comm_waitall_dbl(int* n, MPI_Fint* ireq) {
    int64_t count = static_cast<int64_t>(*n);
    std::vector<MPI_Request> reqs(count);
    for (int64_t i = 0; i < count; ++i) reqs[i] = MPI_Request_f2c(ireq[i]);
    waitall(count, reqs.data());
    for (int64_t i = 0; i < count; ++i) ireq[i] = MPI_Request_c2f(reqs[i]);
  }

  extern "C" void comm_bcast_dbl(double* buf, int* n, int* root, int* col_id, MPI_Fint* fcomm) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    bcast<double>(buf, static_cast<int64_t>(*n), *root, static_cast<int64_t>(*col_id), comm);
  }

  extern "C" void comm_bcastw_dbl(double* buf, int* n, int* root, int* lda, int* lpx,
                     double* work, int* col_id, MPI_Fint* fcomm) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    bcastw<double>(buf, static_cast<int64_t>(*n), static_cast<int64_t>(*root), static_cast<int64_t>(*lda), static_cast<int64_t>(*lpx), work, static_cast<int64_t>(*col_id), comm);
  }
  extern "C" void comm_reduce_dbl(double* buf, double* work, int* n, int* col_id, MPI_Fint* fcomm) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    reduce<double>(buf, work, static_cast<int64_t>(*n), static_cast<int64_t>(*col_id), comm);
  }

  extern "C" void comm_allgather_dbl(double* sendbuf, double* recvbuf, int* n, int* col_id, MPI_Fint* fcomm) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    allgather<double>(sendbuf, recvbuf, static_cast<int64_t>(*n), static_cast<int64_t>(*col_id), comm);
  }

  extern "C" void comm_datacast_dbl(double* u_y, double* u_x, double* u_t, double* u_s, int* n, int* col_id) {
    datacast<double>(u_y, u_x, u_t, u_s, static_cast<int64_t>(*n), static_cast<int64_t>(*col_id));
  }

  extern "C" void comm_datacast_dbl2(double* ur_y, double* ui_y, double* ur_x, double* ui_x,
                     double* u_t, double* u_s, int* n, int* col_id) {
    datacast2<double>(ur_y, ui_y, ur_x, ui_x, u_t, u_s, static_cast<int64_t>(*n), static_cast<int64_t>(*col_id));
  }

  extern "C" void comm_datacast_dblx(int* nk, double* u_y, double* u_x, int* ldv,
                     double* u_t, double* u_s, int* n, int* col_id) {
    datacastx<double>(static_cast<int64_t>(*nk), u_y, u_x, static_cast<int64_t>(*ldv), u_t, u_s,
                       static_cast<int64_t>(*n), static_cast<int64_t>(*col_id));
  }

  extern "C" void comm_print_bcast_algorithm() {
    print_bcast_algorithm();
  }

  extern "C" void comm_print_reduce_algorithm() {
    print_reduce_algorithm();
  }

  extern "C" void comm_print_gather_algorithm() {
    print_gather_algorithm();
  }

  extern "C" void comm_allreduce_binary_sum_dbl(MPI_Fint* fcomm, int* s, double* R, double* R0) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    ALLREDUCE_binary_sum<double>(comm, static_cast<int64_t>(*s), R, R0);
  }

  extern "C" void comm_allreduce_binary_prod_dbl(MPI_Fint* fcomm, int* s, double* R, double* R0) {
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    ALLREDUCE_binary_prod<double>(comm, static_cast<int64_t>(*s), R, R0);
  }

}
