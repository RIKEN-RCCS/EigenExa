!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: my_pdlasrt_mod
!--------*---------*---------*---------*---------*---------*---------*-*

      MODULE MY_PDLASRT_MOD

      USE DLACPY_MOD
      USE COMM_MOD
      USE EIGEN_DEVEL_MOD
      USE EIGEN_DC_MOD
      USE MPI
!$    USE OMP_LIB
      IMPLICIT NONE
      PUBLIC

      CONTAINS

      SUBROUTINE MY_PDLASRT( ID, N, D, Q, IQ, JQ, DESCQ, WORK, LWORK, 
     &     IWORK, LIWORK, INFO )
!     
!     This routine is a modified routine of PDLASRT
!     which is contained in the ScaLAPACK
!     Modification is done by T.Imamura
!     
!     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            IQ, JQ, N
      INTEGER(8) :: INFO, LIWORK, LWORK
!     ..
!     .. Array Arguments ..
      INTEGER            DESCQ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), Q( * ), WORK( * )
!     ..
!     
!     Purpose
!     =======
!     
!     MY_PDLASRT Sort the numbers in D in increasing order and the
!     corresponding vectors in Q.
!     
!     Arguments
!     =========
!     
!     ID      (global input) CHARACTER*1
!     = 'I': sort D in increasing order;
!     = 'D': sort D in decreasing order. (NOT IMPLEMENTED YET)
!     
!     N       (global input) INTEGER
!     The number of columns to be operated on i.e the number of
!     columns of the distributed submatrix sub( Q ). N >= 0.
!     
!     D       (global input/output) DOUBLE PRECISION array, dimmension (N)
!     On exit, the number in D are sorted in increasing order.
!     
!     Q       (local input) DOUBLE PRECISION pointer into the local memory
!     to an array of dimension (LLD_Q, LOCc(JQ+N-1) ). This array
!     contains the local pieces of the distributed matrix sub( A )
!     to be copied from.
!     
!     IQ      (global input) INTEGER
!     The row index in the global array A indicating the first
!     row of sub( Q ).
!     
!     JQ      (global input) INTEGER
!     The column index in the global array A indicating the
!     first column of sub( Q ).
!     
!     DESCQ   (global and local input) INTEGER array of dimension DLEN_.
!     The array descriptor for the distributed matrix A.
!     
!     WORK    (local workspace/local output) DOUBLE PRECISION array,
!     dimension (LWORK)
!     LWORK   (local or global input) INTEGER
!     The dimension of the array WORK.
!     LWORK = MAX( N, NP * ( NB + NQ ))
!     where
!     NP = NUMROC( N, NB, MYROW, IAROW, NPROW ),
!     NQ = NUMROC( N, NB, MYCOL, DESCQ( CSRC_ ), NPCOL )
!     
!     IWORK   (local workspace/local output) INTEGER array,
!     dimension (LIWORK)
!     
!     LIWORK (local or global input) INTEGER
!     The dimension of the array IWORK.
!     LIWORK = N + 2*NB + 2*NPCOL
!     
!     INFO    (global output) INTEGER
!     = 0:  successful exit
!     < 0:  If the i-th argument is an array and the j-entry had
!     an illegal value, then INFO = -(i*100+j), if the i-th
!     argument is a scalar and had an illegal value, then
!     INFO = -i.
!     
!     =====================================================================
!     
!     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     &     MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     &     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     &     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
!     ..
!     .. Local Scalars ..
      INTEGER            CL, COL, DUMMY, I, ICTXT, IID, IIQ, INDCOL,
     &     INDX, INDXC, INDXG, IPQ2, J,
     &     JJQ, K, L, LDQ, LEND, LIWMIN, MYCOL,
     &     MYROW, NB, ND, NP, NPCOL, NPROW, NQ, PSQ, QCOL,
     &     QTOT, SBUF
      INTEGER :: II, IWK_PSQ_MYCOL
      INTEGER :: INFO_I4
      INTEGER(8) :: IPQ, IPW, IPWORK, LWMIN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2L, INDXG2P, NUMROC
      EXTERNAL           INDXG2L, INDXG2P, LSAME, NUMROC
!     ..
!     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PXERBLA, DCOPY,
     &     DGERV2D, DGESD2D, DLAPST
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
!     ..
      integer, allocatable :: RREQ(:)
      integer          :: IWK_(256)
      integer          :: ierr
!     .. Executable Statements ..
!     
!     This is just to keep ftnchek and toolpack/1 happy
!      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
!     &     RSRC_.LT.0 )RETURN
!     
      IF( N.EQ.0 )
     &     RETURN
!     
      CALL BLACS_GRIDINFO( DESCQ( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
!     
!     Test the input parameters
!     
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 600+CTXT_ )
      ELSE
         INFO_I4=INFO
         CALL CHK1MAT( N, 1, N, 1, IQ, JQ, DESCQ, 6, INFO_I4 )
         INFO=INFO_I4
         IF( INFO.EQ.0 ) THEN
            NB = DESCQ( NB_ )
            LDQ = DESCQ( LLD_ )
            NP = NUMROC( N, NB, MYROW, DESCQ( RSRC_ ), NPROW )
            NQ = NUMROC( N, NB, MYCOL, DESCQ( CSRC_ ), NPCOL )
            LWMIN = MAX( N, INT(NP,8)*INT( NB+NQ,8 ) )
            LIWMIN = N + 2*( NB+NPCOL )
            IF( .NOT.LSAME( ID, 'I' ) ) THEN
               INFO = -1
            ELSE IF( N.LT.0 ) THEN
               INFO = -2
            ELSE IF( LWORK.LT.LWMIN ) THEN
               INFO = -9
            ELSE IF( LIWORK.LT.LIWMIN ) THEN
               INFO = -11
            END IF
         END IF
      END IF
!     
      IF( INFO.NE.0 ) THEN
         INFO_I4 = INFO
         CALL PXERBLA( ICTXT, 'MY_PDLASRT', -INFO_I4 )
         INFO=INFO_I4
         RETURN
      END IF
!     
!     Set Pointers
!     
      INDXC = 1
      INDX = INDXC + N
      INDXG = INDX
      INDCOL = INDXG + NB
      QTOT = INDCOL + NB
      PSQ = QTOT + NPCOL
!     
      IID = 1
      IPQ2 = 1
      IPW = IPQ2 + INT(NP,8)*INT(NQ,8)
!     
      DUMMY = 0
      IIQ = INDXG2L( IQ, NB, DUMMY, DUMMY, NPROW )
!     
      allocate ( RREQ(NPCOL), stat=ierr )
      if ( ierr /= 0 ) then
         call eigen_abort( "Memory allocation fault.", 1 )
         return
      end if
!     
!     Sort the eigenvalues in D
!     
      INFO_I4=INFO
      CALL DLAPST( 'I', N, D, IWORK( INDX ), INFO_I4 )
      INFO = INFO_I4
!     
!$OMP PARALLEL DO
!$OMP+         PRIVATE(L)
      DO 10 L = 0, N - 1
         WORK( IID+L ) = D( IWORK( INDX+L ) )
         IWORK( INDXC-1+IWORK( INDX+L ) ) = IID + L
 10   CONTINUE
!$OMP END PARALLEL DO
!     CALL DCOPY( N, WORK, 1, D, 1 )
!$OMP PARALLEL DO
      DO I = 1, N
         D(I) = WORK(I)
      END DO
!$OMP END PARALLEL DO
!     
      ND = 0
 20   CONTINUE
      IF( ND.LT.N ) THEN
         LEND = MIN( NB, N-ND )
         J = JQ + ND
         QCOL = INDXG2P( J, NB, DUMMY, DESCQ( CSRC_ ), NPCOL )

         K = 0
         DO 30 L = 0, LEND - 1
            I = JQ - 1 + IWORK( INDXC+ND+L )
            CL = INDXG2P( I, NB, DUMMY, DESCQ( CSRC_ ), NPCOL )
            IWORK( INDCOL+L ) = CL
            IF( MYCOL.EQ.CL ) THEN
               IWORK( INDXG+K ) = IWORK( INDXC+ND+L )
               K = K + 1
            END IF
 30      CONTINUE
!     
         IF( MYCOL.EQ.QCOL ) THEN
            DO 40 CL = 0, NPCOL - 1
               IWORK( QTOT+CL ) = 0
 40         CONTINUE
!     
            DO 50 L = 0, LEND - 1
               IWORK( QTOT+IWORK( INDCOL+L ) ) =
     &              IWORK( QTOT+IWORK( INDCOL+L ) ) + 1
 50         CONTINUE
!     
            IWORK( PSQ ) = 1
            DO 60 CL = 1, NPCOL - 1
               IWORK( PSQ+CL ) = IWORK( PSQ+CL-1 ) + IWORK( QTOT+CL-1 )
 60         CONTINUE
!     
            DO L = 0, LEND - 1
               CL = IWORK( INDCOL+L )
               I = IWORK( PSQ+CL )
               IWK_(L+1) = I
               IWORK( PSQ+CL ) = I + 1
            END DO
!$OMP PARALLEL DO
!$OMP+         PRIVATE(CL,I,JJQ,IPQ,IPWORK)
            DO L = 0, LEND - 1
               CL = IWORK( INDCOL+L )
               I = JQ + ND + L
               JJQ = INDXG2L( I, NB, DUMMY, DUMMY, NPCOL )
               IPQ = IIQ + INT( JJQ-1,8 )*INT(LDQ,8)
               IPWORK = IPW + INT( IWK_(L+1)-1,8 )*INT(NP,8)
               CALL DCOPY( NP, Q( IPQ ), 1, WORK( IPWORK ), 1 )
            END DO
!$OMP END PARALLEL DO
            IWORK( PSQ ) = 1
!     
            DO 80 CL = 1, NPCOL - 1
               IWORK( PSQ+CL ) = IWORK( PSQ+CL-1 ) + IWORK( QTOT+CL-1 )
 80         CONTINUE
!     
            IWK_PSQ_MYCOL = IWORK( PSQ+MYCOL )
!$OMP PARALLEL DO
!$OMP+         PRIVATE(I,JJQ,IPQ,IPWORK)
            DO L = 0, K - 1
               I = IWORK( INDXG+L )
               JJQ = INDXG2L( I, NB, DUMMY, DUMMY, NPCOL )
               IPQ = IPQ2 + INT( JJQ-1,8 )*INT(NP,8)
               IPWORK = IPW + INT( IWK_PSQ_MYCOL+L-1,8 )*INT(NP,8)
               CALL DCOPY( NP, WORK( IPWORK ), 1, WORK( IPQ ), 1 )
            END DO
!$OMP END PARALLEL DO
            IWORK( PSQ+MYCOL ) = IWK_PSQ_MYCOL + K
!     
            II=0
            DO 100 CL = 1, NPCOL - 1
               COL = MOD( MYCOL+CL, NPCOL )
               SBUF = IWORK( QTOT+COL )
               IF( SBUF.NE.0 ) THEN
                  II = II + 1
                  IPWORK = IPW + INT( IWORK( PSQ+COL )-1,8 )*INT(NP,8)
                  call isend_dbl( WORK( IPWORK ), NP*SBUF,
     &                 COL+1, RREQ(II), y_COMM_WORLD )
               END IF
 100        CONTINUE
            IF ( II > 0 ) call waitall_dbl( II, RREQ )
!     
         ELSE
!     
            IF( K.NE.0 ) THEN
               call recv_dbl( WORK( IPW ), NP*K,
     &              QCOL+1, y_COMM_WORLD )
!$OMP PARALLEL DO
!$OMP+         PRIVATE(L,I,JJQ,IPQ,IPWORK)
               DO 110 L = 0, K - 1
                  I = JQ - 1 + IWORK( INDXG+L )
                  JJQ = INDXG2L( I, NB, DUMMY, DUMMY, NPCOL )
                  IPQ = 1 + INT( JJQ-1,8 )*INT(NP,8)
                  IPWORK = IPW + INT(L,8)*INT(NP,8)
                  CALL DCOPY( NP, WORK( IPWORK ), 1, WORK( IPQ ), 1 )
 110           CONTINUE
!$OMP END PARALLEL DO
            END IF
!     
         END IF
         ND = ND + NB
         GO TO 20
      END IF

      CALL DLACPY( 'Full', NP, NQ, WORK, NP, Q( IIQ ), LDQ )
      deallocate ( RREQ )
!     
!     End of MY_PDLASRT
!     
      END SUBROUTINE MY_PDLASRT

      END MODULE MY_PDLASRT_MOD
