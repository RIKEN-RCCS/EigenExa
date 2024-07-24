#include <stdio.h>
#include <mpi.h>
#include "EigenExa.h"

void test(void)
{
  int n = 2; int nv = n;
  double a[n*n]; int lda = n;
  double w[n];
  double z[n*n]; int ldz = n;
  int mf = 1; int mb = 1;

  int nnod, xnnod, ynnod;
  int inod, xinod, yinod;

  eigen_get_procs(&nnod, &xnnod, &ynnod);
  eigen_get_id(&inod, &xinod, &yinod);

  if ( nnod == 1 ) {

  a[0] = -2;
  a[1] =  1;
  a[2] =  1;
  a[3] = -2;

  eigen_s(n, nv, a, lda, w, z, ldz, mf, mb, "A");

  printf("%le :: ", w[0]);
  printf("%le %le\n", z[0], z[1]);
  printf("%le :: ", w[1]);
  printf("%le %le\n", z[2], z[3]);

  }

  if ( nnod == 4 ) {

  if ( xinod==1 && yinod == 1 ) a[0] = -2;
  if ( xinod==2 && yinod == 1 ) a[0] =  1;
  if ( xinod==1 && yinod == 2 ) a[0] =  1;
  if ( xinod==2 && yinod == 2 ) a[0] = -2;

  eigen_s(n, nv, a, lda, w, z, ldz, mf, mb, "A");

  MPI_Barrier(MPI_COMM_WORLD);
  if(inod==1) printf("%le :: ", w[0]); fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  if ( xinod == 1 && yinod == 1 ) printf("%le ", z[0]); fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  if ( xinod == 2 && yinod == 1 ) printf("%le\n", z[0]); fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  if(inod==1) printf("%le :: ", w[1]); fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  if ( xinod == 1 && yinod == 2 ) printf("%le ", z[0]); fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  if ( xinod == 2 && yinod == 2 ) printf("%le\n", z[0]); fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  }
}

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  eigen_init(MPI_COMM_WORLD,"C");

  int version; char date[32]; char code[32];
  eigen_get_version(&version, date, code);
  printf("Version %06d, %s %s\n", version, date, code);

  test();

  eigen_free();

  return 0;
}

