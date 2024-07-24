#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#pragma once

#ifdef __cplusplus
extern  "C" {
#endif

void eigen_init(MPI_Comm const comm, char const *order);
void eigen_free(void);

void eigen_s(int const n, int const nvec,
		double *a, int const lda,
		double *w,
		double *z, int const ldz,
		int const m_forward,
		int const m_backward,
		char const *mode);
void eigen_sx(int const n, int const nvec,
		double *a, int const lda,
		double *w,
		double *z, int const ldz,
		int const m_forward,
		int const m_backward,
		char const *mode);
void eigen_h(int const n, int const nvec,
		double *a, int const lda,
		double *w,
		double *z, int const ldz,
		int const m_forward,
		int const m_backward,
		char const *mode);

void eigen_get_version(int *version, char *date, char *vcode);
void eigen_get_procs(int *nnod, int *x_nnod, int *y_nnod);
void eigen_get_id(int *inod, int *x_inod, int *y_inod);
void eigen_get_comm(MPI_Comm *comm, MPI_Comm *x_comm, MPI_Comm *y_comm);

void eigen_get_matdims(int const n,
		int *nx, int *ny,
		int const m_forward,
		int const m_backward,
		char const *mode);

#ifdef __cplusplus
}
#endif

