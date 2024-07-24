#include "EigenExa.h"
#include "EigenExa.fh"

#ifdef __cplusplus
extern "C" {
#endif

void eigen_init(MPI_Comm const comm, char const *order)
{
  char order_F[2] = { 0, 0 };
  if ( order != NULL ) { order_F[0] = order[0]; }
  if ( order_F[0] == '\0' ) { order_F[0] = 'C'; }
  MPI_Fint comm_F = MPI_Comm_c2f(comm);
  F_eigen_init (&comm_F, order_F );
}

void eigen_free(void)
{
  F_eigen_free ();
}

void eigen_s(int const n, int const nvec,
	       	double *a, int const lda,
	       	double *w,
	       	double *z, int const ldz,
		int const m_forward,
		int const m_backward,
		char const *mode)
{
  int n_F = n;
  int nvec_F = nvec;
  int lda_F = lda;
  int ldz_F = ldz;
  int m_f_F = m_forward;
  int m_b_F = m_backward;
  char mode_F[2] = { 0, 0 };
  if ( mode != NULL ) { mode_F[0] = mode[0]; }
  if ( mode_F[0] == '\0' ) { mode_F[0] = 'A'; }
  F_eigen_s (&n_F, &nvec_F, a, &lda_F, w, z, &ldz_F, &m_f_F, &m_b_F, mode_F);
}

void eigen_sx(int const n, int const nvec,
	       	double *a, int const lda,
	       	double *w,
	       	double *z, int const ldz,
		int const m_forward,
		int const m_backward,
		char const *mode)
{
  int n_F = n;
  int nvec_F = nvec;
  int lda_F = lda;
  int ldz_F = ldz;
  int m_f_F = m_forward;
  int m_b_F = m_backward;
  char mode_F[2] = { 0, 0 };
  if ( mode != NULL ) { mode_F[0] = mode[0]; }
  if ( mode_F[0] == '\0' ) { mode_F[0] = 'A'; }
  F_eigen_sx (&n_F, &nvec_F, a, &lda_F, w, z, &ldz_F, &m_f_F, &m_b_F, mode_F);
}

void eigen_h(int const n, int const nvec,
	       	double *a, int const lda,
	       	double *w,
	       	double *z, int const ldz,
		int const m_forward,
		int const m_backward,
		char const *mode)
{
  int n_F = n;
  int nvec_F = nvec;
  int lda_F = lda;
  int ldz_F = ldz;
  int m_f_F = m_forward;
  int m_b_F = m_backward;
  char mode_F[2] = { 0, 0 };
  if ( mode != NULL ) { mode_F[0] = mode[0]; }
  if ( mode_F[0] == '\0' ) { mode_F[0] = 'A'; }
  F_eigen_h (&n_F, &nvec_F, a, &lda_F, w, z, &ldz_F, &m_f_F, &m_b_F, mode_F);
}

void eigen_get_version(int *version, char *date, char *vcode)
{
  char date_F[33];
  char vcode_F[33];
  F_eigen_get_version (version, date_F, vcode_F, 32, 32);
  for(int i=31;i>=0;i--) if( date_F[i]!=' ' ) { break; } else { date_F[i]=0; }
  for(int i=31;i>=0;i--) if( vcode_F[i]!=' ' ) { break; } else { vcode_F[i]=0; }
  strncpy(date, date_F, 32);
  strncpy(vcode, vcode_F, 32);
}

void eigen_get_procs(int *nnod, int *x_nnod, int *y_nnod)
{
  F_eigen_get_procs (nnod, x_nnod, y_nnod);
}

void eigen_get_id(int *inod, int *x_inod, int *y_inod)
{
  F_eigen_get_id (inod, x_inod, y_inod);
}

void eigen_get_comm(MPI_Comm *comm, MPI_Comm *x_comm, MPI_Comm *y_comm)
{
  MPI_Fint comm_F, comm_x, comm_y;
  F_eigen_get_id (&comm_F, &comm_x, &comm_y);
  *comm = MPI_Comm_f2c(comm_F);
  *x_comm = MPI_Comm_f2c(comm_x);
  *y_comm = MPI_Comm_f2c(comm_y);
}

void eigen_get_matdims(int const n,
                int *nx, int *ny,
                int const m_forward,
                int const m_backward,
                char const *mode)
{
  int n_F = n;
  int m_f_F = m_forward;
  int m_b_F = m_backward;
  char mode_F[2] = { 0, 0 };
  if ( mode != NULL ) { mode_F[0] = mode[0]; }
  if ( mode_F[0] == '\0' ) { mode_F[0] = 'O'; }
  F_eigen_get_matdims (&n_F, nx, ny, &m_f_F, &m_b_F, mode_F);
}

#undef	F_mod

#ifdef __cplusplus
}
#endif


