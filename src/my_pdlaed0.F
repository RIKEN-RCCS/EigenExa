!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: my_pdlaed0_mod
!--------*---------*---------*---------*---------*---------*---------*-*

      MODULE MY_PDLAED0_MOD

      USE MY_PDLAED1_MOD
      USE MY_PDLAEDZ_MOD
      USE LAPACK_EIGEN_MOD
      USE COMM_MOD
      USE EIGEN_DEVEL_MOD
      USE EIGEN_DC_MOD
      USE MPI
!$    USE OMP_LIB
      IMPLICIT NONE
      PUBLIC

      CONTAINS

      SUBROUTINE MY_PDLAED0(N, HBW, D, E, LDE, Q, IQ, JQ,
     &     DESCQ, WORK, IWORK, CP, CQ, INFO)
!     
      
!     This routine is a modified routine of PDLAED0
!     which is contained in the ScaLAPACK
!     Modification is done by T.Imamura
!     
!     .. Scalar Arguments ..
      INTEGER            IQ, JQ, N, CP, CQ, HBW, LDE
      INTEGER :: LDQ
      INTEGER(8) :: INFO
!     ..
!     .. Array Arguments ..
      INTEGER            DESCQ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), Q( * ), WORK( * )
      DOUBLE PRECISION   E( LDE, HBW )
!     ..
!     
!     Purpose
!     =======
!     
!     MY_PDLAED0 computes all eigenvalues and corresponding eigenvectors of a
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
!     E1      (global input/output) DOUBLE PRECISION array, dimension (N-1)
!     On entry, the subdiagonal elements of the tridiagonal matrix.
!     On exit, E has been destroyed.

!     (added by Phuong) 
!     E2          (global input) DOUBLE PRECISION array, dimention (N-2)
!     On entry, the second subdiagonal of the matrix
!     On exit, E has been destroyed(?)


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
     &     RSRC_ = 7, CSRC_ = 8, LLD_ = 9  )

!     ..
!     .. Local Scalars ..
      INTEGER            I, ID, IDCOL, IDROW, IID, IIQ, IM1, IM2,
     &     IQCOL, IQROW, J, K, JJD, JJQ, MATSIZ,
     &     MYCOL, MYROW, N1, NB, NBL, NBL1, NPCOL, NPROW,
     &     SUBPBS, TSUBPBS, ISTAGE
      integer :: ierr
      INTEGER(8) :: IPQ
      INTEGER :: INFO_I4

      DOUBLE PRECISION   RHO
!     .. Local Arrays ..

!     INTEGER          DESCZ( 50 )

      real(8) :: C(HBW,HBW), S(HBW), U(HBW,HBW), VVT(HBW, HBW),
     &     VT(HBW,HBW), UUT(HBW,HBW), V(HBW+HBW, HBW)
      real(8), allocatable ::  W(:)
!     ..
!     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGEBR2D, DGEBS2D, DGERV2D,
     &     DGESD2D, DSTEQR, INFOG2L, PXERBLA !, MAKELAMDAX
!     ..
!     .. External Functions ..
!     ..
!     DOUBLE PRECISION               DNRM2
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MIN, MAX

      real(8) :: t1,t2,d1,d2
      integer :: MY_RANK, IC, REQ(10)
!     ..
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
         CALL PXERBLA( DESCQ( CTXT_ ), 'MY_PDLAED0', -INFO_I4 )
         INFO = INFO_I4
         RETURN
      END IF

!     
!     
      allocate( W(N), stat=ierr )
      if ( ierr /= 0 ) then
         call eigen_abort( "Memory allocation fault.", 1 )
         return
      end if
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
         
         C(1:HBW, 1:HBW) = 0

         DO J = 1, HBW
            DO K=1, MIN(J, N-IM1)
               C(HBW+K-J, K) = E(IM1-J+K, J)
            END DO      
         END DO

!     SVD : C = U*S*VT
         INFO_I4 = INFO
         CALL DGESVD('A', 'A', HBW, HBW, C, HBW, S, U, HBW,  
     &        VT, HBW, WORK, 7*HBW, INFO_I4) ! WORK size should be much la
         INFO = INFO_I4
!     rger than 7HBW

         DO  J = 1, HBW
            S(J) = SQRT(S(J))
         END DO
         DO  J = 1, HBW
            DO  K = 1, HBW
               U(K, J)  = U(K, J)  * S(J)
               VT(K,J) =  VT(K,J)  * S(K)
            END DO
         END DO

         CALL DGEMM('N', 'T', HBW, HBW, HBW, ONE, U, HBW, U, 
     &        HBW, ZERO, UUT, HBW)

         CALL DGEMM('T', 'N', HBW, HBW, HBW, ONE, VT, HBW, VT, 
     &        HBW, ZERO, VVT, HBW)
         
         DO J = 1, HBW
            D(IM1-HBW+J) = D(IM1-HBW+J) - UUT(J, J)
            IF  ( IM1+J <= N ) THEN
               D(IM1+J)     = D(IM1+J)     - VVT(J, J)
            END IF
         END DO
         
         DO J = 1, HBW-1
            DO  K = HBW-J, 1, -1
               E(IM1-K-J+1, J) = E(IM1-K-J+1, J) 
     &              - UUT(HBW-K+1, 1-K+HBW-J) 
               IF  ( IM1+HBW-K+1 <= N ) THEN
                  E(IM1+HBW-K-J+1,J) = E(IM1+HBW-K-J+1,J)
     &                 - VVT(HBW-K+1, 1-K+HBW-J) 
               END IF
            END DO
         END DO    


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
            
            CALL LAPACK_EIGEN2(N, MATSIZ, HBW, ID, D( ID ), 
     &           E( ID, 1 ), LDE, Q( IPQ ), LDQ)
            
!     =            IF( INFO.NE.0 ) THEN
!     =               CALL PXERBLA( DESCQ( CTXT_ ), 'DSTEQR', -INFO )
!     =C-
!     =               deallocate( W )
!     =C-
!     =               RETURN
!     =            END IF

            WORK(ID:ID+MATSIZ-1) = D(ID:ID+MATSIZ-1)

         ELSE
            WORK(ID:ID+MATSIZ-1) = ZERO
         END IF
      END DO
      call reduce_dbl( WORK, D, N, 1, TRD_COMM_WORLD )
      D(1:N) = WORK(1:N)
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

!     CALL DESCINIT( DESCZ, MATSIZ, 1, MATSIZ, 1, 0, 0,
!     $                        DESCQ(CTXT_),  MATSIZ, INFO_I4)    

!     
!     Merge lower order eigensystems (of size N1 and MATSIZ - N1)
!     into an eigensystem of size MATSIZ.
!     
            IM1 =  ID+N1 - 1

            C(1:HBW, 1:HBW) = 0

            DO J = 1, HBW
               DO K=1, MIN(J, N-IM1)
                  C(HBW+K-J, K) = E(IM1-J+K, J)
               END DO
            END DO

            INFO_I4 = INFO
            CALL DGESVD('A', 'A', HBW, HBW, C, HBW, S, U, HBW, VT,  
     &           HBW, WORK, 7*HBW, INFO_I4)
            INFO = INFO_I4

            DO J = 1, HBW
               S(J)  = SQRT(S(J))
            END DO
            DO J = 1, HBW
               DO K = 1, HBW
                  U(K, J) = U(K,J)  * S(J)
                  VT(K,J) = VT(K,J) * S(K)
               END DO
            END DO

            DO K = 1, HBW
               DO J = 1, HBW
                  V(J,K) = U(J,K);   V(HBW+J, K) = VT(K,J)
               END DO
            END DO


            DO K=1, HBW

               ISTAGE = MIN(K, 2)
               d1 = MPI_WTIME()
!$OMP PARALLEL
               CALL MY_PDLAEDZ(ISTAGE,
     &              HBW, MATSIZ, N1, ID, Q, IQ, JQ, LDQ, 
     &              DESCQ, W, WORK (1 + N), V(1,K))
!$OMP END PARALLEL
               d2 = MPI_WTIME()
               p_timez = p_timez + (d2-d1)

               RHO = DNRM2(MATSIZ, W, 1)
               IF(RHO .EQ. 0) CYCLE
               CALL DSCAL( MATSIZ, 1D0/RHO, W, 1)
               RHO = RHO * RHO

               CALL MY_PDLAED1( MATSIZ, N1, D( ID ), ID, Q, IQ, JQ,  
     &              DESCQ, RHO, WORK, IWORK( SUBPBS+1 ), W, ISTAGE,
     &              INFO )

            END DO
            
            IF( INFO.NE.0 ) THEN
               INFO = INFO*( N+1 ) + ID
            END IF

 70         CONTINUE
            IWORK( I / 2+1 ) = IWORK( I+2 )
 80      CONTINUE
         SUBPBS = SUBPBS / 2

         GO TO 60

      END IF
!     
!     end while
!     
 90   CONTINUE

      deallocate( W )

      RETURN
!     
!     End of MY_PDLAED0
!     
      END SUBROUTINE MY_PDLAED0

      END MODULE MY_PDLAED0_MOD
