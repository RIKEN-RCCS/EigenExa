!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: my_pdlaed3_mod
!--------*---------*---------*---------*---------*---------*---------*-*

      MODULE MY_PDLAED3_MOD

      USE EIGEN_DEVEL_MOD
      USE EIGEN_DC_MOD
      USE MPI
!$    USE OMP_LIB
      USE comm_mod
      IMPLICIT NONE
      PUBLIC

      CONTAINS

      SUBROUTINE MY_PDLAED3( ICTXT, K, N, NB, D, DROW, DCOL,
     &     RHO, DLAMDA,
     &     W, Z, U, LDU, BUF, INDX, INDCOL, INDROW,
     &     INDXR, INDXC, CTOT, NPCOL_, INFO )
!     
!     This routine is a modified routine of PDLAED3
!     which is contained in the ScaLAPACK
!     Modification is done by T.Imamura
!     
!     .. Scalar Arguments ..
      INTEGER            DCOL, DROW, ICTXT, K, LDU, N, NB, NPCOL_
      INTEGER(8) :: INFO
      DOUBLE PRECISION   RHO
!     ..
!     .. Array Arguments ..
      INTEGER            CTOT( 0: NPCOL_-1, 4 ), INDCOL( * ),
     &     INDROW( * ), INDX( * ), INDXC( * ), INDXR( * )
      DOUBLE PRECISION   BUF( 1:3*K ), D( * ), DLAMDA( * ), U( LDU, * ),
     &     W( * ), Z(1:K)
!     ..
!     
!     Purpose
!     =======
!     
!     MY_PDLAED3 finds the roots of the secular equation, as defined by the
!     values in D, W, and RHO, between 1 and K.  It makes the
!     appropriate calls to SLAED4
!     
!     This code makes very mild assumptions about floating point
!     arithmetic. It will work on machines with a guard digit in
!     add/subtract, or on those binary machines without guard digits
!     which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
!     It could conceivably fail on hexadecimal or decimal machines
!     without guard digits, but we know of none.
!     
!     Arguments
!     =========
!     
!     ICTXT  (global input) INTEGER
!     The BLACS context handle, indicating the global context of
!     the operation on the matrix. The context itself is global.
!     
!     K      (output) INTEGER
!     The number of non-deflated eigenvalues, and the order of the
!     related secular equation. 0 <= K <=N.
!     
!     N      (input) INTEGER
!     The dimension of the symmetric tridiagonal matrix.  N >= 0.
!     
!     NB      (global input) INTEGER
!     The blocking factor used to distribute the columns of the
!     matrix. NB >= 1.
!     
!     D      (input/output) DOUBLE PRECISION array, dimension (N)
!     On entry, D contains the eigenvalues of the two submatrices to
!     be combined.
!     On exit, D contains the trailing (N-K) updated eigenvalues
!     (those which were deflated) sorted into increasing order.
!     
!     DROW   (global input) INTEGER
!     The process row over which the first row of the matrix D is
!     distributed. 0 <= DROW < NPROW.
!     
!     DCOL   (global input) INTEGER
!     The process column over which the first column of the
!     matrix D is distributed. 0 <= DCOL < NPCOL.
!     
!     Q      (input/output) DOUBLE PRECISION array, dimension (LDQ, N)
!     On entry, Q contains the eigenvectors of two submatrices in
!     the two square blocks with corners at (1,1), (N1,N1)
!     and (N1+1, N1+1), (N,N).
!     On exit, Q contains the trailing (N-K) updated eigenvectors
!     (those which were deflated) in its last N-K columns.
!     
!     LDQ    (input) INTEGER
!     The leading dimension of the array Q.  LDQ >= max(1,NQ).
!     
!     RHO    (global input/output) DOUBLE PRECISION
!     On entry, the off-diagonal element associated with the rank-1
!     cut which originally split the two submatrices which are now
!     being recombined.
!     On exit, RHO has been modified to the value required by
!     MY_PDLAED3.
!     
!     DLAMDA (global output) DOUBLE PRECISION array, dimension (N)
!     A copy of the first K eigenvalues which will be used by
!     SLAED3 to form the secular equation.
!     
!     W      (global output) DOUBLE PRECISION array, dimension (N)
!     The first k values of the final deflation-altered z-vector
!     which will be passed to SLAED3.
!     
!     Z      (global input) DOUBLE PRECISION array, dimension (N)
!     On entry, Z contains the updating vector (the last
!     row of the first sub-eigenvector matrix and the first row of
!     the second sub-eigenvector matrix).
!     On exit, the contents of Z have been destroyed by the updating
!     process.
!     
!     U     (global output) DOUBLE PRECISION array
!     global dimension (N, N), local dimension (LDU, NQ).
!     Q  contains the orthonormal eigenvectors of the symmetric
!     tridiagonal matrix.
!     
!     LDU    (input) INTEGER
!     The leading dimension of the array U.
!     
!     QBUF   (workspace) DOUBLE PRECISION array, dimension 3*N
!     
!     
!     INDX   (workspace) INTEGER array, dimension (N)
!     The permutation used to sort the contents of DLAMDA into
!     ascending order.
!     
!     INDCOL (workspace) INTEGER array, dimension (N)
!     
!     
!     INDROW (workspace) INTEGER array, dimension (N)
!     
!     
!     INDXR (workspace) INTEGER array, dimension (N)
!     
!     
!     INDXC (workspace) INTEGER array, dimension (N)
!     
!     CTOT   (workspace) INTEGER array, dimension( NPCOL, 4)
!     
!     NPCOL   (global input) INTEGER
!     The total number of columns over which the distributed
!     submatrix is distributed.
!     
!     INFO   (output) INTEGER
!     = 0:  successful exit.
!     < 0:  if INFO = -i, the i-th argument had an illegal value.
!     > 0:  The algorithm failed to compute the ith eigenvalue.
!     
!     =====================================================================
!     
!     .. Parameters ..
!     DOUBLE PRECISION   ONE
!     PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            COL, GI, I, IINFO, IIU, IU, J, JJU, JU,
     &     KK, KL, KLC, KLR, MYCOL, MYKL, MYKLR, MYROW,
     &     NPCOL, NPROW, PDC, PDR, ROW
      DOUBLE PRECISION   AUX, TEMP
      INTEGER(8) :: IPD, IPD0
      INTEGER :: SINFO
!     ..
!     .. Shared Scalar ..
!$    INTEGER, SAVE :: KLR_SHARED
!     .. Local Arrays ..
      REAL(8), allocatable :: SZ(:), SBUF(:)
!     ..
!     .. External Functions ..
      INTEGER            INDXG2L
      DOUBLE PRECISION   DLAMC3 ! , DNRM2
      EXTERNAL           INDXG2L, DLAMC3 ! , DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DCOPY, DGEBR2D, DGEBS2D,
     &     DGERV2D, DGESD2D, DLAED4
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MOD, SIGN, SQRT
!     ..
      integer :: ERR_MPI
!     .. Executable Statements ..
!     
!     Test the input parameters.
!     
      INFO = 0
!     
!     Quick return if possible
!     
      IF( K.EQ.0 ) GOTO 190
!     
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
!     
      ROW = DROW
      COL = DCOL
!$OMP DO SCHEDULE(STATIC,1)
      DO I = 1, N, NB
         ROW = MOD(DROW+(I-1)/NB, NPROW)
         COL = MOD(DCOL+(I-1)/NB, NPCOL)
         DO J = 0, NB - 1
            IF( I+J.LE.N ) THEN
               INDROW( I+J ) = ROW
               INDCOL( I+J ) = COL
            END IF
         END DO
!     ROW = MOD( ROW+1, NPROW )
!     COL = MOD( COL+1, NPCOL )
      END DO
!$OMP END DO
      ROW = MOD(DROW+(NB-1)/NB+1, NPROW)
      COL = MOD(DCOL+(NB-1)/NB+1, NPCOL)
!     
      MYKL = CTOT( MYCOL, 1 ) + CTOT( MYCOL, 2 ) + CTOT( MYCOL, 3 )
      KLR = MYKL / NPROW
      IF( MYROW.EQ.DROW ) THEN
         MYKLR = KLR + MOD( MYKL, NPROW )
      ELSE
         MYKLR = KLR
      END IF
!     
      PDC = 1
      COL = DCOL
      DO
         IF( MYCOL == COL ) EXIT
         PDC = PDC + CTOT( COL, 1 ) + CTOT( COL, 2 ) + CTOT( COL, 3 )
         COL = MOD( COL+1, NPCOL )
      END DO
!     30 CONTINUE
!     IF( MYCOL.NE.COL ) THEN
!     PDC = PDC + CTOT( COL, 1 ) + CTOT( COL, 2 ) + CTOT( COL, 3 )
!     COL = MOD( COL+1, NPCOL )
!     GO TO 30
!     END IF
!     
      PDR = PDC
      KL = KLR + MOD( MYKL, NPROW )
      ROW = DROW
      DO
         IF( MYROW == ROW ) EXIT
         PDR = PDR + KL
         KL = KLR
         ROW = MOD( ROW+1, NPROW )
      END DO
!     40 CONTINUE
!     IF( MYROW.NE.ROW ) THEN
!     PDR = PDR + KL
!     KL = KLR
!     ROW = MOD( ROW+1, NPROW )
!     GO TO 40
!     END IF
!     
!---  
!     
!$OMP DO
      DO I = 1, K
         DLAMDA( I ) = DLAMC3( DLAMDA( I ), DLAMDA( I ) ) - DLAMDA( I )
      END DO
!$OMP END DO
!     
!---  
!     
      SINFO = INFO
      IF( MYKLR.GT.0 ) THEN
!     
         ALLOCATE ( SZ(1:K) )
         SZ( 1:K ) = ONE
         ALLOCATE ( SBUF(1:K) )

!$OMP DO SCHEDULE(STATIC,1)
         DO I = 1, MYKLR
            KK = PDR + I - 1
            CALL DLAED4( K, KK, DLAMDA, W, SBUF, RHO, AUX, IINFO )
            BUF(K+I) = AUX
            IF( IINFO.NE.0 ) THEN
               SINFO = KK
            END IF
!     
!     ..Compute part of z
!     
#if 0
            DO J = 1, KK - 1
               SZ( J ) = SZ( J )*( SBUF( J ) /
     &              ( DLAMDA( J )-DLAMDA( KK ) ) )
            END DO
            SZ( KK ) = SZ( KK )*SBUF( KK )
            DO J = KK + 1, K
               SZ( J ) = SZ( J )*( SBUF( J ) /
     &              ( DLAMDA( J )-DLAMDA( KK ) ) )
            END DO
#else
!OCL NOFP_RELAXED
!OCL NOFP_CONTRACT
!OCL NOEVAL
            DO J=1,K
               TEMP = DLAMDA( J )-DLAMDA( KK )
               IF ( J == KK ) TEMP = ONE
               SZ( J ) = SZ( J ) * ( SBUF( J ) / TEMP )
            ENDDO
#endif
         END DO
!$OMP END DO

         DEALLOCATE ( SBUF )

!$OMP MASTER
         Z(1:K) = SZ(1:K)
!
!
! count up the flops on the Loewner law's update
!
      flops = flops + DBLE(MYKLR) * DBLE(K*3)
!
!
!$OMP END MASTER
!$OMP BARRIER
!$       DO I = 1, OMP_GET_NUM_THREADS()-1
!$          IF ( OMP_GET_THREAD_NUM() == I ) THEN
!$             Z(1:K) = Z(1:K) * SZ(1:K)
!$          END IF
!$OMP BARRIER
!$       END DO
!$OMP BARRIER

         DEALLOCATE ( SZ )

      ELSE

!$OMP DO
         DO I=1, K
            Z( I ) = ONE
         END DO
!$OMP END DO

      END IF

!$OMP MASTER
      INFO = SINFO
!$OMP END MASTER
!$OMP BARRIER
!$    DO I = 1, OMP_GET_NUM_THREADS()-1
!$       IF ( OMP_GET_THREAD_NUM() == I ) THEN
!$          IF ( SINFO /= 0 ) INFO = SINFO
!$       END IF
!$OMP BARRIER
!$    END DO
!$OMP BARRIER

!$OMP MASTER

#ifdef ALLREDUCE_MPI
      call MPI_ALLREDUCE( Z( 1 ), BUF( 1 ),
     &     K, MPI_DOUBLE_PRECISION,
     &     MPI_PROD, TRD_COMM_WORLD, ERR_MPI )
#else
          call ALLREDUCE_binary_prod(
     &         TRD_COMM_WORLD,K, Z( 1 ), BUF( 1 ))
#endif

!$OMP END MASTER
!$OMP BARRIER
!$OMP DO
      DO I = 1, K
         Z( I ) = SIGN( SQRT( -BUF( I ) ), W( I ) )
      END DO
!$OMP END DO
!     
!---  
!     
!$OMP DO
      DO I = 1, MYKL
         BUF(2*K+I) = ZERO
      END DO
!$OMP END DO

      if ( myrow == drow ) then
         IPD = (2*K+1)
      else
         IPD = (2*K+1) + MOD( MYKL, NPROW )
     &        + INT(MOD( MYROW-DROW+NPROW, NPROW ),8) * INT(KLR,8)
      end if

      IF ( MYKLR > 0 ) THEN
!$OMP DO
         DO I = 1, MYKLR
            BUF(IPD+I-1) = BUF(K+I)
         END DO
!$OMP END DO
      END IF

!$OMP MASTER
      IF ( MYKL > 0 ) THEN

#ifdef ALLREDUCE_MPI
         call MPI_ALLREDUCE( BUF( 2*K+1 ), BUF( K+1 ),
     &        MYKL, MPI_DOUBLE_PRECISION,
     &        MPI_SUM, x_COMM_WORLD, ERR_MPI )
#else
          call ALLREDUCE_binary_sum(
     &         x_COMM_WORLD,MYKL, BUF( 2*K+1 ), BUF( K+1 ))

#endif
      END IF
!$OMP END MASTER
!$OMP BARRIER
!     
!---  
!     
!$OMP DO
      DO I = 1, K
         BUF(2*K+I) = ZERO
      END DO
!$OMP END DO

      IPD = (2*K+1)
      IPD0 = IPD ! just initialized
      DO I = 0, NPCOL - 1
         COL = MOD( DCOL+I, NPCOL )
         IF ( mycol == col ) IPD0 = IPD
         KL = CTOT( COL, 1 ) + CTOT( COL, 2 ) + CTOT( COL, 3 )
         IPD = IPD + KL
      END DO
      IPD = IPD0

      IF ( MYKL > 0 ) THEN
!$OMP DO
         DO I = 1, MYKL
            BUF(IPD+I-1) = BUF(K+I)
         END DO
!$OMP END DO
      END IF

!$OMP MASTER

#ifdef ALLREDUCE_MPI
      call MPI_ALLREDUCE( BUF( 2*K+1 ), BUF( K+1 ),
     &     K, MPI_DOUBLE_PRECISION,
     &     MPI_SUM, y_COMM_WORLD, ERR_MPI )
#else
          call ALLREDUCE_binary_sum(
     &     y_COMM_WORLD,K,BUF( 2*K+1 ), BUF( K+1 ))

#endif
!     
!     Copy of D at the good place
!     
      KLC = 0
      KLR = 0
      DO I = 1, K
         GI = INDX( I )
         D( GI ) = BUF( K+I )
         COL = INDCOL( GI )
         ROW = INDROW( GI )
         IF( COL.EQ.MYCOL ) THEN
            KLC = KLC + 1
            INDXC( KLC ) = I
         END IF
         IF( ROW.EQ.MYROW ) THEN
            KLR = KLR + 1
            INDXR( KLR ) = I
         END IF
      END DO
!     
!     Compute eigenvectors of the modified rank-1 modification.
!     
!     
!     Since the result of KLR on the master thread is accessed from the other
!     threads, it is broadcast via a shared variable KLR_SHARED hereafter.
!     
!$    KLR_SHARED = KLR
!$OMP END MASTER
!$OMP BARRIER
!$    KLR = KLR_SHARED
!$OMP BARRIER
!     
!---  
!     
      IF( MYKL.NE.0 ) THEN

         ALLOCATE( SBUF(1:K) )

!$OMP DO
         DO J = 1, MYKL
            KK = INDXC( J )
            JU = INDX( KK )
            JJU = INDXG2L( JU, NB, J, J, NPCOL )
            CALL DLAED4( K, KK, DLAMDA, W, SBUF, RHO, AUX, IINFO )
            IF( K == 1 .OR. K == 2 ) THEN
               DO I = 1, KLR
                  KK = INDXR( I )
                  IU = INDX( KK )
                  IIU = INDXG2L( IU, NB, J, J, NPROW )
                  U( IIU, JJU ) = SBUF( KK )
               END DO
               CYCLE
            END IF
!     
            SBUF( 1:K ) = Z( 1:K ) / SBUF( 1:K )
            TEMP = DNRM2( K, SBUF, 1 )
!     
            DO I = 1, KLR
               KK = INDXR( I )
               IU = INDX( KK )
               IIU = INDXG2L( IU, NB, J, J, NPROW )
               U( IIU, JJU ) = SBUF( KK ) / TEMP
            END DO
         END DO
!$OMP END DO

         DEALLOCATE( SBUF )

      END IF
!     
!---  
!     
 190  CONTINUE
!     
      RETURN
!     
!     End of MY_PDLAED3
!     
      END SUBROUTINE MY_PDLAED3

      END MODULE MY_PDLAED3_MOD
