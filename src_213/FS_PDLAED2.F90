!>
!> @file   FS_PDLAED2.F90
!> @brief  subroutine FS_PDLAED2
!!>
!
!
!
!> subroutine FS_PDLAED2
!>
!> @brief @n
!>  Purpose @n
!>  ======= @n
!>  FS_PDLAED2 sorts the two sets of eigenvalues together into a single  @n
!>  sorted set.  Then it tries to deflate the size of the problem.       @n
!>  There are two ways in which deflation can occur:  when two or more   @n
!>  eigenvalues are close together or if there is a tiny entry in the    @n
!>  Z vector.  For each such occurrence the order of the related secular
!>  equation problem is reduced by one.
!
!  Arguments
!  =========
!>
!> @param[out]    K        (output) INTEGER @n
!>                         The number of non-deflated eigenvalues, and the order of the @n
!>                         related secular equation. 0 <= K <=N.
!>
!> @param[in]     N        (input) INTEGER @n
!>                         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!>
!> @param[in]     N1       (input) INTEGER @n
!>                         The location of the last eigenvalue in the leading sub-matrix. @n
!>                         min(1,N) <= N1 <= N.
!>
!> @param[in,out] D        (input/output) DOUBLE PRECISION array, dimension (N)           @n
!>                         On entry, D contains the eigenvalues of the two submatrices to @n
!>                         be combined.                                                   @n
!>                         On exit, D contains the trailing (N-K) updated eigenvalues     @n
!>                         (those which were deflated) sorted into increasing order.
!>
!> @param[in,out] Q        (input/output) DOUBLE PRECISION array, dimension (LDQ, NQ)  @n
!>                         On entry, Q contains the eigenvectors of two submatrices in @n
!>                         the two square blocks with corners at (1,1), (N1,N1)        @n
!>                         and (N1+1, N1+1), (N,N).                                    @n
!>                         On exit, Q contains the trailing (N-K) updated eigenvectors @n
!>                         (those which were deflated) in its last N-K columns.
!>
!> @param[in]     LDQ      (local input) INTEGER @n
!>                         The leading dimension of the array Q.  LDQ >= max(1,NP).
!>
!> @param[in]     SUBTREE  (input) type(bt_node) @n
!>                         sub-tree information of merge block.
!>
!> @param[in,out] RHO      (global input/output) DOUBLE PRECISION                        @n
!>                         On entry, the off-diagonal element associated with the rank-1 @n
!>                         cut which originally split the two submatrices which are now  @n
!>                         being recombined.                                             @n
!>                         On exit, RHO has been modified to the value required by       @n
!>                         FS_PDLAED3.
!>
!> @param[in,out] Z        (global input) DOUBLE PRECISION array, dimension (N)           @n
!>                         On entry, Z contains the updating vector (the last             @n
!>                         row of the first sub-eigenvector matrix and the first row of   @n
!>                         the second sub-eigenvector matrix).                            @n
!>                         On exit, the contents of Z have been destroyed by the updating @n
!>                         process.
!>
!> @param[out]    W        (global output) DOUBLE PRECISION array, dimension (N)      @n
!>                         The first k values of the final deflation-altered z-vector @n
!>                         which will be passed to FS_PDLAED3.
!>
!> @param[out]    DLAMDA   (global output) DOUBLE PRECISION array, dimension (N)   @n
!>                         A copy of the first K eigenvalues which will be used by @n
!>                         FS_PDLAED3 to form the secular equation.
!>
!> @param[out]    Q2       (output) DOUBLE PRECISION array, dimension (LDQ2, NQ) @n
!>                         The eigen vectors which sorted by COLTYP
!>
!> @param[in]     LDQ2     (input) INTEGER @n
!>                         The leading dimension of the array Q2.
!>
!> @param[out]    INDX     (output) INTEGER array, dimension (N)                    @n
!>                         The permutation used to sort the contents of DLAMDA into @n
!>                         ascending order which will be passed to FS_PDLAED3.
!>
!> @param[out]    CTOT     (output) INTEGER array, dimension (NPCOL, 4) @n
!>                         The number of COLTYP of each process column  @n
!>                         which will be passed to FS_PDLAED3.
!>
!> @param         QBUF     (workspace) DOUBLE PRECISION array, dimension (N)
!>
!> @param         COLTYP   (workspace) INTEGER array, dimension (N)                   @n
!>                         During execution, a label which will indicate which of the @n
!>                         following types a column in the Q2 matrix is:              @n
!>                         1 : non-zero in the upper half only;                       @n
!>                         2 : dense;                                                 @n
!>                         3 : non-zero in the lower half only;                       @n
!>                         4 : deflated.
!>
!> @param         INDCOL   (workspace) INTEGER array, dimension (N)
!>
!> @param         INDXC    (workspace) INTEGER array, dimension (N)                       @n
!>                         The permutation used to arrange the columns of the deflated    @n
!>                         Q matrix into three groups:  the first group contains non-zero @n
!>                         elements only at and above N1, the second contains             @n
!>                         non-zero elements only below N1, and the third is dense.
!>
!> @param         INDXP    (workspace) INTEGER array, dimension (N)                      @n
!>                         The permutation used to place deflated values of D at the end @n
!>                         of the array.  INDXP(1:K) points to the nondeflated D-values  @n
!>                         and INDXP(K+1:N) points to the deflated eigenvalues.
!>
!> @param         PSM      (workspace) INTEGER array, dimension (NPCOL, 4)
!>
!> @param[out]    INFO     (global output) INTEGER @n
!>                         = 0: successful exit   @n
!>                         /=0: error exit
!>
!> @param[out]    prof     (global output) type(FS_prof) @n
!>                         profiling information of each subroutines.
!>
!> @note This routine is modified from ScaLAPACK PDLAED2.f
!>
subroutine FS_PDLAED2( K, N, N1, D, Q, LDQ, SUBTREE, RHO, Z, W, &
                       DLAMDA, Q2, LDQ2, INDX, CTOT, &
                       QBUF, COLTYP, INDCOL, INDXC, INDXP, PSM, &
                       INFO, prof )
      use FS_const_mod
      use FS_libs_mod
      use FS_prof_mod
      use FS_dividing_mod
      use comm_mod
      implicit none
!
!     .. Scalar Arguments ..
      integer          , intent(out)   :: K
      integer          , intent(in)    :: N, N1, LDQ, LDQ2
      type(bt_node)    , intent(in)    :: SUBTREE
      real(kind(0.0d0)), intent(inout) :: RHO
      integer          , intent(out)   :: INFO
      type(FS_prof)    , intent(inout) :: prof
!
!     .. Array Arguments ..
      real(kind(0.0d0)), intent(inout) :: D(N)
      real(kind(0.0d0)), intent(inout) :: Q(LDQ,*)
      real(kind(0.0d0)), intent(inout) :: Z(N)
      real(kind(0.0d0)), intent(out)   :: W(N)
      real(kind(0.0d0)), intent(out)   :: DLAMDA(N)
      real(kind(0.0d0)), intent(out)   :: Q2(LDQ2,*)
      integer          , intent(out)   :: INDX(N)
      integer          , intent(out)   :: CTOT(0:SUBTREE%y_nnod-1,4)
      ! work
      real(kind(0.0d0)) :: QBUF(N)
      integer :: COLTYP(N), INDCOL(N), INDXC(N), INDXP(N)
      integer :: PSM(0:SUBTREE%y_nnod-1,4)
!     ..
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
      integer :: I,J, JS, JJS, JJQ2
      integer :: COL, CT
      integer :: NBLK, NB, NP, NPROW, NPCOL, MYROW, MYCOL
      integer :: N2, N1P1
      integer :: K2,NJ,PJ,NJJ,PJJ,NJCOL,PJCOL
      integer :: IMAX, JMAX
      real(kind(0.0d0)) :: T, EPS, TOL, S, C, TAU
      integer :: NPA, ROW
!     ..
!     .. Local Arrays ..
      integer :: PTT(4)
      integer :: REQ(2)
!     ..
!     .. External Functions ..
      integer           :: IDAMAX
      real(kind(0.0d0)) :: DLAPY2, DLAMCH

      integer :: ierr

      external   IDAMAX, DLAPY2, DLAMCH
!     ..
!     .. External Subroutines ..
      external   DLAPST
!     ..
!     .. Executable Statements ..
!
      INFO = 0
#ifdef _DEBUGLOG
!      CALL MPI_Barrier(SUBTREE%MERGE_COMM,ierr)
      if( FS_MYRANK.eq.0 ) then
        write(*,'(a)') "FS_PDLAED2 start."
      endif
#endif
!
#if TIMER_PRINT
      CALL FS_prof_start(prof, 50)
#endif
!
!     Quick return if possible
!
      IF( N.EQ.0 ) THEN
        GO TO 220
      ENDIF
!
      CALL FS_GRIDINFO(SUBTREE, NPROW, NPCOL, MYROW, MYCOL)
      NBLK = FS_get_NBLK(SUBTREE)
      NB   = FS_get_NB(SUBTREE)
      NP   = (NBLK / NPROW) * NB
!
      N2   = N - N1
      N1P1 = N1 + 1
!
      DLAMDA(1:N) = 0.d0
!
      IF( RHO.LT.ZERO ) THEN
         CALL DSCAL( N2, MONE, Z( N1P1 ), 1 )
      END IF
!     
!     Normalize z so that norm(z) = 1.  Since z is the concatenation of
!     two normalized vectors, norm2(z) = sqrt(2).
!     
      T = ONE / SQRT( TWO )
      CALL DSCAL( N, T, Z, 1 )
!     
!     RHO = ABS( norm(z)**2 * RHO )
!     
      RHO = ABS( TWO*RHO )
!     
!     Calculate the allowable deflation tolerance
!     
      IMAX = IDAMAX( N, Z, 1 )
      JMAX = IDAMAX( N, D, 1 )
!     EPS = PDLAMCH( ICTXT, 'Epsilon' )
      EPS = DLAMCH( 'Epsilon' )
      TOL = EIGHT*EPS*MAX( ABS( D( JMAX ) ), ABS( Z( IMAX ) ) )
!     
!     If the rank-1 modifier is small enough, no more needs to be done
!     except to reorganize Q so that its columns correspond with the
!     elements in D.
!     
      IF( RHO*ABS( Z( IMAX ) ).LE.TOL ) THEN
         K = 0
         GO TO 220
      END IF
!     
!     If there are multiple eigenvalues then the problem deflates.  Here
!     the number of equal eigenvalues are found.  As each equal
!     eigenvalue is found, an elementary reflector is computed to rotate
!     the corresponding eigensubspace so that the corresponding
!     components of Z are zero in this new basis.
!     
!     
      CALL DLAPST( 'I', N, D, INDX, INFO )
!
!$OMP PARALLEL PRIVATE(I,J,COL)
!$OMP DO
      DO 10 I = 1, N1
         COLTYP( I ) = 1
 10   CONTINUE
!$OMP END DO
!$OMP DO
      DO 20 I = N1P1, N
         COLTYP( I ) = 3
 20   CONTINUE
!$OMP END DO
!$OMP DO
      DO 40 I = 1, N, NB
         CALL FS_INFOG1L('C', I, SUBTREE, J, COL)
         DO 30 J = 0, NB - 1
            IF( I+J.LE.N ) &
               INDCOL( I+J ) = COL
 30      CONTINUE
 40   CONTINUE
!$OMP END DO
!$OMP END PARALLEL
      NPA = 0
!$OMP PARALLEL DO PRIVATE(I,J,ROW) REDUCTION(+:NPA)
      DO I = 1, N, NB
         CALL FS_INFOG1L('R', I, SUBTREE, J, ROW)
         IF( ROW.EQ.MYROW ) THEN
            DO J = 0, NB - 1
               IF( I+J.LE.N ) THEN
                 NPA = NPA + 1
               ENDIF
           END DO
         ENDIF
      END DO
!$OMP END PARALLEL DO
!
      K = 0
      K2 = N + 1
      DO 50 J = 1, N
         NJ = INDX( J )
         IF( RHO*ABS( Z( NJ ) ).LE.TOL ) THEN
!     
!     Deflate due to small z component.
!     
            K2 = K2 - 1
            COLTYP( NJ ) = 4
            INDXP( K2 ) = NJ
            IF( J.EQ.N ) &
               GO TO 80
         ELSE
            PJ = NJ
            GO TO 60
         END IF
 50   CONTINUE
 60   CONTINUE
      J = J + 1
      IF( J.GT.N ) &
          GO TO 80
      NJ = INDX( J )
      IF( RHO*ABS( Z( NJ ) ).LE.TOL ) THEN
!     
!     Deflate due to small z component.
!     
         K2 = K2 - 1
         COLTYP( NJ ) = 4
         INDXP( K2 ) = NJ
      ELSE
!     
!     Check if eigenvalues are close enough to allow deflation.
!     
         S = Z( PJ )
         C = Z( NJ )
!     
!     Find sqrt(a**2+b**2) without overflow or
!     destructive underflow.
!     
         TAU = DLAPY2( C, S )
         T = D( NJ ) - D( PJ )
         C = C / TAU
         S = -S / TAU
         IF( ABS( T*C*S ).LE.TOL ) THEN
!     
!     Deflation is possible.
!     
            Z( NJ ) = TAU
            Z( PJ ) = ZERO
            IF( COLTYP( NJ ).NE.COLTYP( PJ ) ) &
               COLTYP( NJ ) = 2
            COLTYP( PJ ) = 4

!$OMP PARALLEL
!$OMP MASTER
            CALL FS_INFOG1L('C', NJ, SUBTREE, NJJ, NJCOL)
            CALL FS_INFOG1L('C', PJ, SUBTREE, PJJ, PJCOL)

            IF( INDCOL( PJ ).EQ.INDCOL( NJ ) .AND. MYCOL.EQ.NJCOL ) THEN
               CALL DROT( NPA, Q( 1, PJJ ), 1, Q( 1, NJJ ), 1, C, S )
            ELSE IF( MYCOL.EQ.PJCOL ) THEN
!               write(*,*)"@ 1 communication @"
#if 0
               CALL irecv_dbl( QBUF,      NPA, &
                    NJCOL+1, REQ(2), SUBTREE%MERGE_COMM_Y )
               CALL isend_dbl( Q(1, PJJ), NPA, &
                    NJCOL+1, REQ(1), SUBTREE%MERGE_COMM_Y )
#else
               CALL irecv_dbl( QBUF,      NPA, &
                    SUBTREE%group_Y_processranklist(NJCOL+1) + 1, &
                    REQ(2), FS_COMM_WORLD )
               CALL isend_dbl( Q(1, PJJ), NPA, &
                    SUBTREE%group_Y_processranklist(NJCOL+1) + 1, &
                    REQ(1), FS_COMM_WORLD )
#endif
               CALL waitall_dbl( 2, REQ )
               CALL DROT( NPA, Q( 1, PJJ ), 1, QBUF, 1, C, S )
            ELSE IF( MYCOL.EQ.NJCOL ) THEN
!               write(*,*)"@ 2 communication @"
#if 0
               CALL irecv_dbl( QBUF,      NPA, &
                    PJCOL+1, REQ(2), SUBTREE%MERGE_COMM_Y )
               CALL isend_dbl( Q(1, NJJ), NPA, &
                    PJCOL+1, REQ(1), SUBTREE%MERGE_COMM_Y )
#else

               CALL irecv_dbl( QBUF,      NPA, &
                    SUBTREE%group_Y_processranklist(PJCOL+1) + 1, &
                    REQ(2), FS_COMM_WORLD )
               CALL isend_dbl( Q(1, NJJ), NPA, &
                    SUBTREE%group_Y_processranklist(PJCOL+1) + 1, &
                    REQ(1), FS_COMM_WORLD )
#endif
               CALL waitall_dbl( 2, REQ )
               CALL DROT( NPA, QBUF, 1, Q( 1, NJJ ), 1, C, S )
            END IF
!$OMP END MASTER
!$OMP SINGLE
            T = D( PJ )*C**2 + D( NJ )*S**2
            D( NJ ) = D( PJ )*S**2 + D( NJ )*C**2
            D( PJ ) = T

            K2 = K2 - 1
            I = 1
 70         CONTINUE
            IF( K2+I.LE.N ) THEN
               IF( D( PJ ).LT.D( INDXP( K2+I ) ) ) THEN
                  INDXP( K2+I-1 ) = INDXP( K2+I )
                  INDXP( K2+I ) = PJ
                  I = I + 1
                  GO TO 70
               ELSE
                  INDXP( K2+I-1 ) = PJ
               END IF
            ELSE
               INDXP( K2+I-1 ) = PJ
            END IF
!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL 

            PJ = NJ
         ELSE
            K = K + 1
            DLAMDA( K ) = D( PJ )
            W( K ) = Z( PJ )
            INDXP( K ) = PJ
            PJ = NJ
         END IF
      END IF
      GO TO 60
 80   CONTINUE
!     
!     Record the last eigenvalue.
!     
      K = K + 1
      DLAMDA( K ) = D( PJ )
      W( K ) = Z( PJ )
      INDXP( K ) = PJ
!     
!     Count up the total number of the various types of columns, then
!     form a permutation which positions the four column types into
!     four uniform groups (although one or more of these groups may be
!     empty).
!     
!$OMP PARALLEL DO PRIVATE(I, J)
      DO 100 J = 1, 4
         DO 90 I = 0, NPCOL - 1
            CTOT( I, J ) = 0
 90      CONTINUE
         PTT( J ) = 0
 100  CONTINUE
!$OMP END PARALLEL DO
      DO 110 J = 1, N
         CT  = COLTYP( J )
         COL = INDCOL( J )
         CTOT( COL, CT ) = CTOT( COL, CT ) + 1
 110  CONTINUE
!     
!     PSM(*) = Position in SubMatrix (of types 1 through 4)
!     
!$OMP PARALLEL DO PRIVATE(COL)
      DO 120 COL = 0, NPCOL - 1
         PSM( COL, 1 ) = 1
         PSM( COL, 2 ) = 1 + CTOT( COL, 1 )
         PSM( COL, 3 ) = PSM( COL, 2 ) + CTOT( COL, 2 )
         PSM( COL, 4 ) = PSM( COL, 3 ) + CTOT( COL, 3 )
 120  CONTINUE
!$OMP END PARALLEL DO
!     
      PTT( 1 ) = 1
!$OMP PARALLEL DO PRIVATE(I,J,CT)
      DO I = 2, 4
         CT = 0
         DO J = 0, NPCOL - 1
            CT = CT + CTOT( J, I-1 )
         END DO
         PTT(I) = CT
      END DO
!$OMP END PARALLEL DO
      DO I = 2, 4
         PTT(I) = PTT( I-1 ) + PTT(I)
      END DO
!     
!     Fill out the INDXC array so that the permutation which it induces
!     will place all type-1 columns first, all type-2 columns next,
!     then all type-3's, and finally all type-4's.
!     
      DO 150 J = 1, N
         JS = INDXP( J )
         COL = INDCOL( JS )
         CT = COLTYP( JS )
!        I = INDXL2G( PSM( COL, CT ), NB, COL, DCOL, NPCOL )
         I = FS_INDXL2G( 'C', PSM( COL, CT ), COL, SUBTREE )
         INDX( J ) = I
         INDXC( PTT( CT ) ) = I
         PSM( COL, CT ) = PSM( COL, CT ) + 1
         PTT( CT ) = PTT( CT ) + 1
 150  CONTINUE
!
!$OMP PARALLEL DO PRIVATE(I,J,JS,JJS,JJQ2,COL)
      DO 160 J = 1, N
         JS = INDXP( J )
         COL = INDCOL( JS )
         IF( COL.EQ.MYCOL ) THEN
            JJS  = FS_INDXG2L( 'C', JS, SUBTREE )
            I    = INDX( J )
            JJQ2 = FS_INDXG2L( 'C', I, SUBTREE )
            CALL DCOPY( NPA, Q( 1, JJS ), 1, Q2( 1, JJQ2 ), 1 )
         END IF
 160  CONTINUE
!$OMP END PARALLEL DO
!     
!     The deflated eigenvalues and their corresponding vectors go back
!     into the last N - K slots of D and Q respectively.
!     
      CALL DCOPY( N, D, 1, Z, 1 )
!$OMP PARALLEL DO PRIVATE(I,J,JS)
      DO 170 J = K + 1, N
         JS = INDXP( J )
         I = INDX( J )
         D( I ) = Z( JS )
 170  CONTINUE
!$OMP END PARALLEL DO
!
 220  CONTINUE

#if TIMER_PRINT
      CALL FS_prof_end(prof, 50)
#endif
#ifdef _DEBUGLOG
!      CALL MPI_Barrier(SUBTREE%MERGE_COMM,ierr)
      if( FS_MYRANK.eq.0 ) then
        write(*,'(a,i0)') "FS_PDLAED2 end. INFO=", INFO
      endif
#endif
      RETURN
!
!     End of FS_PDLAED2
!     
end subroutine FS_PDLAED2

