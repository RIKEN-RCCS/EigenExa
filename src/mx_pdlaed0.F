!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: mx_pdlaed0_mod
!--------*---------*---------*---------*---------*---------*---------*-*

      MODULE MX_PDLAED0_MOD

      USE MX_PDLAED1_MOD
      USE COMM_MOD
      USE EIGEN_DEVEL_MOD
      USE EIGEN_DC_MOD
      USE MPI
!$    USE OMP_LIB
      IMPLICIT NONE
      PUBLIC

      CONTAINS

      SUBROUTINE MX_PDLAED0( N, D, E, Q, IQ, JQ, DESCQ,
     &     WORK, IWORK, INFO )
!     .. Scalar Arguments ..
      INTEGER            IQ, JQ, N
      INTEGER(8)         INFO
!     ..
!     .. Array Arguments ..
      INTEGER            DESCQ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), Q( * ), WORK( * )
!     ..
!     
!     Purpose
!     =======
!     
!     MX_PDLAED0 computes all eigenvalues and corresponding eigenvectors of a
!     symmetric tridiagonal matrix using the divide and conquer method.
!     
!     
!     Arguments
!     =========
!     
!     N       (global input) INTEGER
!     The order of the tridiagonal matrix T.  N >= 0.
!     
!     D       (global input/output) DOUBLE PRECISION array, dimension (N)
!     On entry, the diagonal elements of the tridiagonal matrix.
!     On exit, if INFO = 0, the eigenvalues in descending order.
!     
!     E       (global input/output) DOUBLE PRECISION array, dimension (N-1)
!     On entry, the subdiagonal elements of the tridiagonal matrix.
!     On exit, E has been destroyed.
!     
!     Q       (local output) DOUBLE PRECISION array,
!     global dimension (N, N),
!     local dimension ( LLD_Q, LOCc(JQ+N-1))
!     Q  contains the orthonormal eigenvectors of the symmetric
!     tridiagonal matrix.
!     On output, Q is distributed across the P processes in block
!     cyclic format.
!     
!     IQ      (global input) INTEGER
!     Q's global row index, which points to the beginning of the
!     submatrix which is to be operated on.
!     
!     JQ      (global input) INTEGER
!     Q's global column index, which points to the beginning of
!     the submatrix which is to be operated on.
!     
!     DESCQ   (global and local input) INTEGER array of dimension DLEN_.
!     The array descriptor for the distributed matrix Z.
!     
!     
!     WORK    (local workspace ) DOUBLE PRECISION array, dimension (LWORK)
!     LWORK = 6*N + 2*NP*NQ, with
!     NP = NUMROC( N, MB_Q, MYROW, IQROW, NPROW )
!     NQ = NUMROC( N, NB_Q, MYCOL, IQCOL, NPCOL )
!     IQROW = INDXG2P( IQ, NB_Q, MYROW, RSRC_Q, NPROW )
!     IQCOL = INDXG2P( JQ, MB_Q, MYCOL, CSRC_Q, NPCOL )
!     
!     IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
!     LIWORK = 2 + 7*N + 8*NPCOL
!     
!     INFO    (global output) INTEGER
!     = 0:  successful exit
!     < 0:  If the i-th argument is an array and the j-entry had
!     an illegal value, then INFO = -(i*100+j), if the i-th
!     argument is a scalar and had an illegal value, then
!     INFO = -i.
!     > 0:  The algorithm failed to compute the INFO/(N+1) th
!     eigenvalue while working on the submatrix lying in
!     global rows and columns mod(INFO,N+1).
!     
!     =====================================================================
!     
!     .. Parameters ..
!     
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     &     MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     &     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     &     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, ID, IDCOL, IDROW, IID, IINFO, IIQ, IM1, IM2,
     &     IQCOL, IQROW, J, JJD, JJQ, LDQ, MATSIZ,
     &     MYCOL, MYROW, N1, NB, NBL, NBL1, NPCOL, NPROW,
     &     SUBPBS, TSUBPBS
      INTEGER(8) :: IPQ
      INTEGER :: INFO_I4
!     ..
!     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGEBR2D, DGEBS2D, DGERV2D,
     &     DGESD2D, DSTEQR, INFOG2L, PXERBLA
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MIN
!     ..
      REAL(8) :: d1,d2
!     .. Executable Statements ..
!     
!     This is just to keep ftnchek and toolpack/1 happy
!      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
!     &     RSRC_.LT.0 )RETURN
!     
!     Test the input parameters.
!     
      CALL BLACS_GRIDINFO( DESCQ( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      INFO = 0
      IF( DESCQ( NB_ ).GT.N .OR. N.LT.2 )
     &     INFO = -1
      IF( INFO.NE.0 ) THEN
         INFO_I4 = INFO
         CALL PXERBLA( DESCQ( CTXT_ ), 'MX_PDLAED0', -INFO_I4 )
         INFO = INFO_I4
         RETURN
      END IF
!     
      NB = DESCQ( NB_ )
      LDQ = DESCQ( LLD_ )
      CALL INFOG2L( IQ, JQ, DESCQ, NPROW, NPCOL, MYROW, MYCOL, IIQ, JJQ,
     &     IQROW, IQCOL )
!     
!     Determine the size and placement of the submatrices, and save in
!     the leading elements of IWORK.
!     
      TSUBPBS = ( N-1 ) / NB + 1
      IWORK( 1 ) = TSUBPBS
      SUBPBS = 1
 10   CONTINUE
      IF( IWORK( SUBPBS ).GT.1 ) THEN
         DO 20 J = SUBPBS, 1, -1
            IWORK( 2*J ) = ( IWORK( J )+1 ) / 2
            IWORK( 2*J-1 ) = IWORK( J ) / 2
 20      CONTINUE
         SUBPBS = 2*SUBPBS
         GO TO 10
      END IF
      DO 30 J = 2, SUBPBS
         IWORK( J ) = IWORK( J ) + IWORK( J-1 )
 30   CONTINUE
!     
!     Divide the matrix into TSUBPBS submatrices of size at most NB
!     using rank-1 modifications (cuts).
!     
      DO 40 I = NB + 1, N, NB
         IM1 = I - 1
         D( IM1 ) = D( IM1 ) - ABS( E( IM1 ) )
         D( I ) = D( I ) - ABS( E( IM1 ) )
 40   CONTINUE
!     
!     Solve each submatrix eigenproblem at the bottom of the divide and
!     conquer tree. D is the same on each process.
!     
      d1 = MPI_WTIME()
      DO ID = 1, N, NB
         CALL INFOG2L( IQ-1+ID, JQ-1+ID, DESCQ, NPROW, NPCOL, MYROW,
     &        MYCOL, IID, JJD, IDROW, IDCOL )
         MATSIZ = MIN( NB, N-ID+1 )
         IF( MYROW.EQ.IDROW .AND. MYCOL.EQ.IDCOL ) THEN
            IPQ = IID + INT( JJD-1,8 )*INT(LDQ,8)
            INFO_I4 = INFO
            CALL DSTEQR( 'I', MATSIZ, D( ID ), E( ID ), Q( IPQ ), LDQ,
     &           WORK, INFO_I4 )
            INFO = INFO_I4
            IF( INFO.NE.0 ) THEN
               INFO_I4 = INFO
               CALL PXERBLA( DESCQ( CTXT_ ), 'DSTEQR', -INFO_I4 )
               INFO = INFO_I4
               RETURN
            END IF
         ELSE
            D(ID:ID+MATSIZ-1) = ZERO
         END IF
      END DO
!     
      call reduce_dbl( D, WORK, N, 1, TRD_COMM_WORLD )
      d2 = MPI_WTIME()
      p_time0 = p_time0 + (d2-d1)
!     
!     Successively merge eigensystems of adjacent submatrices
!     into eigensystem for the corresponding larger matrix.
!     
!     while ( SUBPBS > 1 )
!     
 60   CONTINUE
      IF( SUBPBS.GT.1 ) THEN
         IM2 = SUBPBS - 2
         DO 80 I = 0, IM2, 2
            IF( I.EQ.0 ) THEN
               NBL = IWORK( 2 )
               NBL1 = IWORK( 1 )
               IF( NBL1.EQ.0 )
     &              GO TO 70
               ID = 1
               MATSIZ = MIN( N, NBL*NB )
               N1 = NBL1*NB
            ELSE
               NBL = IWORK( I+2 ) - IWORK( I )
               NBL1 = NBL / 2
               IF( NBL1.EQ.0 )
     &              GO TO 70
               ID = IWORK( I )*NB + 1
               MATSIZ = MIN( NB*NBL, N-ID+1 )
               N1 = NBL1*NB
            END IF
!     
!     Merge lower order eigensystems (of size N1 and MATSIZ - N1)
!     into an eigensystem of size MATSIZ.
!     
            CALL MX_PDLAED1( MATSIZ, N1, D( ID ), ID, Q, IQ, JQ, DESCQ,
     &           E( ID+N1-1 ), WORK, IWORK( SUBPBS+1 ), INFO )
            IF( INFO.NE.0 ) THEN
               INFO = INT(INFO,8)*INT( N+1,8 ) + ID
            END IF
!     
 70         CONTINUE
            IWORK( I / 2+1 ) = IWORK( I+2 )
 80      CONTINUE
         SUBPBS = SUBPBS / 2
!     
         GO TO 60
      END IF
!     
!     end while
!     
 90   CONTINUE
      RETURN
!     
!     End of MX_PDLAED0
!     
      END SUBROUTINE MX_PDLAED0

      END MODULE MX_PDLAED0_MOD
