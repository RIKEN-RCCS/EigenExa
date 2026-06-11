!>
!> @file   FS_PDLASRT.F90
!> @brief  subroutine FS_PDLASRT
!!>
!
!
!
!> subroutine FS_PDLASRT
!>
!> @brief @n
!>  Purpose @n
!>  ======= @n
!>  FS_PDLASRT Sort the numbers in D in increasing order and the @n
!>  corresponding vectors in Q.
!
!  Arguments
!  =========
!>
!> @param[in]     N        (global input) INTEGER @n
!>                         The number of columns to be operated on i.e the number of @n
!>                         columns of the distributed submatrix sub( Q ). N >= 0.
!>
!> @param[in,out] D        (global input/output) DOUBLE PRECISION array, dimension (N) @n
!>                         On exit, the number in D are sorted in increasing order.
!>
!> @param[in,out] Q        (input/output) DOUBLE PRECISION pointer into the local memory @n
!>                         to an array of dimension (LDQ, NQ). This array contains the   @n
!>                         local pieces of the distributed matrix sub( A ) to be copied  @n
!>                         from.
!>
!> @param[in]     LDQ      (local input) INTEGER @n
!>                         The leading dimension of the array Q.  LDQ >= max(1,NP).
!>
!> @param[in]     SUBTREE  (input) type(bt_node) @n
!>                         sub-tree information of merge block.
!>
!> @param         Q2       (workspace) DOUBLE PRECISION array, dimension (LDQ2, NQ)
!>
!> @param[in]     LDQ2     (input) INTEGER @n
!>                         The leading dimension of the array Q2. (=NP)
!>
!> @param         SENDQ    (workspace) DOUBLE PRECISION array, dimension (LDQ2*NQ)
!>
!> @param         RECVQ    (workspace) DOUBLE PRECISION array, dimension (LDQ2*NQ)
!>
!> @param         BUF      (workspace) DOUBLE PRECISION array, dimension (N)
!>
!> @param         INDROW   (workspace) INTEGER array, dimension (N)
!>
!> @param         INDCOL   (workspace) INTEGER array, dimension (N)
!>
!> @param         INDX     (workspace) INTEGER array, dimension (N)
!>
!> @param         INDRCV   (workspace) INTEGER array, dimension (N)
!>
!> @param[out]    INFO     (global output) INTEGER @n
!>                         = 0: successful exit    @n
!>                         /=0: error exit
!>
!> @param[out]    prof     (global output) type(FS_prof) @n
!>                         profiling information of each subroutines.
!>
!> @note This routine is modified from ScaLAPACK PDLASRT.f
!>

subroutine FS_PDLASRT( N, D, Q, LDQ, SUBTREE, Q2, LDQ2, SENDQ, RECVQ, &
                       BUF, INDROW, INDCOL, INDX, INDRCV, INFO, prof )
      use FS_libs_mod
      use FS_dividing_mod
      implicit none
!
!     .. Scalar Arguments ..
      integer          , intent(in)    :: N, LDQ, LDQ2
      type(bt_node)    , intent(in)    :: SUBTREE
      integer          , intent(out)   :: INFO
      type(FS_prof)    , intent(inout) :: prof
!
!     .. Array Arguments ..
      real(kind(0.0d0)), intent(inout) :: Q(LDQ,*)
      real(kind(0.0d0)), intent(inout) :: D(N)
      ! work
      real(kind(0.0d0)) :: Q2(LDQ2,*)
      real(kind(0.0d0)) :: SENDQ(*)
      real(kind(0.0d0)) :: RECVQ(*)
      real(kind(0.0d0)) :: BUF(*)
      integer :: INDROW(*), INDCOL(*), INDX(*), INDRCV(*)
!     ..
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
      integer :: I, J
      integer :: NBLK, NB, NP, NQ
      integer :: NPROW, NPCOL, MYROW, MYCOL
      integer :: ROW, COL, PJ, PJCOL, PJROW, GI
      integer :: NSEND, NRECV
      integer :: JCOL, JL, IROW, IL
      integer :: req, ierr, stat(MPI_STATUS_SIZE)
!     ..
!     .. Local Arrays ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
#ifdef _DEBUGLOG
      if( FS_MYRANK.eq.0 ) then
        write(*,'(a)') "FS_PDLASRT start."
      endif
#endif
#if TIMER_PRINT
      CALL FS_prof_start(prof, 70)
#endif
!
      CALL FS_GRIDINFO(SUBTREE, NPROW, NPCOL, MYROW, MYCOL)

      NBLK  = FS_get_NBLK(SUBTREE)
      NB    = FS_get_NB(SUBTREE)
      NP    = (NBLK / NPROW) * NB
      NQ    = (NBLK / NPCOL) * NB
!
!$OMP PARALLEL DO PRIVATE(I,J,COL,ROW) SCHEDULE(STATIC,1)
      DO I = 1, N, NB
         CALL FS_INFOG1L('R', I, SUBTREE, J, ROW)
         CALL FS_INFOG1L('C', I, SUBTREE, J, COL)
         DO J = 0, NB - 1
            IF( I+J.LE.N ) THEN
               INDROW( I+J ) = ROW
               INDCOL( I+J ) = COL
            ENDIF
         END DO
      END DO
!$OMP END PARALLEL DO
!
!     Sort the eigenvalues in D
!
      CALL DLAPST( 'I', N, D, INDX, INFO )
!
!$OMP PARALLEL PRIVATE(I)
!$OMP DO
      DO I = 1, N
        BUF( I ) = D( INDX(I) )
      END DO
!$OMP END DO
!$OMP DO
      DO I = 1, N
        D( I ) = BUF( I )
      END DO
!$OMP END DO
!$OMP END PARALLEL

      ! 列方向の入れ替え
      ! 固有値昇順に合わせたソートと非ブロックサイクリック化
      DO PJ = 0, NPCOL-1

        ! 入れ替え対象のプロセス列
        ! デッドロックしないように逆順に回す
        PJCOL = MOD( NPCOL-1-MYCOL-PJ+NPCOL, NPCOL )

        ! 
        NSEND = 0
        NRECV = 0
        DO J = 1, N

          ! 元のグローバルインデクス
          GI  = INDX(J)

          ! 元のグローバルインデクスを保持するプロセス列
          COL = INDCOL(GI)

          ! 入れ替え先インデクスJを保持するプロセス列
          JCOL = MOD(J-1, NPCOL)

          ! PJCOLに送信する列
          IF( COL.EQ.MYCOL .AND. JCOL.EQ.PJCOL ) THEN

            ! ローカルインデクス
            JL = FS_INDXG2L( 'C', GI, SUBTREE )

            ! 送信バッファに格納
            NSEND = NSEND + 1
            CALL DCOPY( NP, Q(1, JL), 1, SENDQ((NSEND-1)*NP+1), 1 )

          END IF

          ! PJCOLから受信する列
          IF( COL.EQ.PJCOL .AND. JCOL.EQ.MYCOL ) THEN

            ! 受信数
            NRECV = NRECV + 1

            ! 格納先のローカルインデクス
            INDRCV(NRECV) = (J+NPCOL-1)/NPCOL
          END IF

        END DO

        ! irecv
        IF( NRECV.GT.0 ) THEN
#if 0
          call MPI_IRECV(RECVQ, NP*NRECV, MPI_DOUBLE_PRECISION, PJCOL, &
                         1, SUBTREE%MERGE_COMM_Y, req, ierr)
#else
          call MPI_IRECV(RECVQ, NP*NRECV, MPI_DOUBLE_PRECISION, &
               SUBTREE%group_Y_processranklist(PJCOL+1), &
               1, FS_COMM_WORLD, req, ierr)

#endif
        END IF

        ! send
        IF( NSEND.GT.0 ) THEN
#if 0
          call MPI_SEND(SENDQ, NP*NSEND, MPI_DOUBLE_PRECISION, PJCOL, &
                        1, SUBTREE%MERGE_COMM_Y, ierr)
#else
          call MPI_SEND(SENDQ, NP*NSEND, MPI_DOUBLE_PRECISION, &
               SUBTREE%group_Y_processranklist(PJCOL+1), &
               1, FS_COMM_WORLD, ierr)
#endif
        END IF

        ! waitと展開
        IF( NRECV.GT.0 ) THEN
          call MPI_WAIT(req, stat, ierr)

!$OMP PARALLEL DO PRIVATE(J, JL)
          DO J = 1, NRECV
            JL = INDRCV(J)
            CALL DCOPY( NP, RECVQ((J-1)*NP+1), 1, Q2(1, JL), 1 )
          END DO
!$OMP END PARALLEL DO
        END IF

      END DO

      ! 行方向の入れ替え
      ! 非ブロックサイクリック化
      DO PJ = 0, NPROW-1

        ! 入れ替え対象のプロセス行
        ! デッドロックしないように逆順に回す
        PJROW = MOD( NPROW-1-MYROW-PJ+NPROW, NPROW )

        ! 
        NSEND = 0
        NRECV = 0
        DO I = 1, N

          ! 元のグローバルインデクスを保持するプロセス行
          ROW = INDROW(I)

          ! 入れ替え先インデクスIを保持するプロセス行
          IROW = MOD(I-1, NPROW)

          ! PJROWに送信する行
          IF( ROW.EQ.MYROW .AND. IROW.EQ.PJROW ) THEN

            ! ローカルインデクス
            IL = FS_INDXG2L( 'R', I, SUBTREE )

            ! 送信バッファに格納
            NSEND = NSEND + 1
            CALL DCOPY( NQ, Q2(IL, 1), LDQ2, SENDQ((NSEND-1)*NQ+1), 1 )

          END IF

          ! PJROWから受信する列
          IF( ROW.EQ.PJROW .AND. IROW.EQ.MYROW ) THEN

            ! 受信数
            NRECV = NRECV + 1

            ! 格納先のローカルインデクス
            INDRCV(NRECV) = (I+NPROW-1)/NPROW
          END IF

        END DO

        ! irecv
        IF( NRECV.GT.0 ) THEN
#if 0
          call MPI_IRECV(RECVQ, NRECV*NQ, MPI_DOUBLE_PRECISION, PJROW, &
                         1, SUBTREE%MERGE_COMM_X, req, ierr)
#else
          call MPI_IRECV(RECVQ, NRECV*NQ, MPI_DOUBLE_PRECISION, &
               SUBTREE%group_X_processranklist(PJROW+1), &
               1, FS_COMM_WORLD, req, ierr)

#endif
        END IF

        ! send
        IF( NSEND.GT.0 ) THEN
#if 0
          call MPI_SEND(SENDQ, NSEND*NQ, MPI_DOUBLE_PRECISION, PJROW, &
                        1, SUBTREE%MERGE_COMM_X, ierr)
#else
          call MPI_SEND(SENDQ, NSEND*NQ, MPI_DOUBLE_PRECISION,&
               SUBTREE%group_X_processranklist(PJROW+1), &
               1, FS_COMM_WORLD, ierr)
#endif
        END IF

        ! waitと展開
        IF( NRECV.GT.0 ) THEN
          call MPI_WAIT(req, stat, ierr)

!$OMP PARALLEL DO PRIVATE(I, IL)
          DO I = 1, NRECV
            IL = INDRCV(I)
            CALL DCOPY( NQ, RECVQ((I-1)*NQ+1), 1, Q(IL, 1), LDQ )
          END DO
!$OMP END PARALLEL DO
        END IF

      END DO
!
#if TIMER_PRINT
      CALL FS_prof_end(prof, 70)
#endif
#ifdef _DEBUGLOG
      if( FS_MYRANK.eq.0 ) then
        write(*,'(a)') "FS_PDLASRT end."
      endif
#endif
!
      RETURN
!
!     End of FS_PDLASRT
!
end subroutine FS_PDLASRT

