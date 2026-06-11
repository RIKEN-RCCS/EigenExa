!>
!> @file   FS_PDLAED3.F90
!> @brief  subroutine FS_PDLAED3
!!>
!
!
!
!> subroutine FS_PDLAED3
!>
!> @brief @n
!>  Purpose @n
!>  ======= @n
!>  FS_PDLAED3 finds the roots of the secular equation, as defined by the @n
!>  values in D, W, and RHO, between 1 and K.  It makes the               @n
!>  appropriate calls to DLAED4                                           @n
!>                                                                        @n
!>  The final stage consists of computing the updated eigenvectors        @n
!>  directly using the updated eigenvalues.  The eigenvectors for         @n
!>  the current problem are multiplied with the eigenvectors from         @n
!>  the overall problem.
!
!  Arguments
!  =========
!>
!> @param[in]     K        (output) INTEGER @n
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
!> @param[in,out] D        (input/output) DOUBLE PRECISION array, dimension (N)        @n
!>                         On entry, D contains the trailing (N-K) updated eigenvalues @n
!>                         (those which were deflated) sorted into increasing order.   @n
!>                         On exit, D contains the updated eigenvalues.
!>
!> @param[in]     RHO      (global output) DOUBLE PRECISION                        @n
!>                         The off-diagonal element associated with the rank-1 cut @n
!>                         cut which modified to the value for this subroutine.
!>
!> @param[in,out] DLAMDA   (global input) DOUBLE PRECISION array, dimension (K)   @n
!>                         A copy of the first K eigenvalues.
!>
!> @param[in]     W        (global input) DOUBLE PRECISION array, dimension (K)        @n
!>                         The first k values of the final deflation-altered z-vector. @n
!>
!> @param[out]    Q        (output) DOUBLE PRECISION array, dimension (LDQ, NQ)  @n
!>                         On exit, Q contains the updated eigenvectors.
!>
!> @param[in]     LDQ      (local input) INTEGER @n
!>                         The leading dimension of the array Q.  LDQ >= max(1,NP).
!>
!> @param[in]     SUBTREE  (input) type(bt_node) @n
!>                         sub-tree information of merge block.
!>
!> @param[in,out] Q2       (input/workspace) DOUBLE PRECISION array, dimension (LDQ2, NQ) @n
!>                         On entry, The eigen vectors which sorted by COLTYP             @n
!>                         On exit, Q2 has been destroyed.
!>
!> @param[in]     LDQ2     (input) INTEGER @n
!>                         The leading dimension of the array Q2.
!>
!> @param         U        (workspace) DOUBLE PRECISION array, dimension (LDU, NQ) @n
!>                         delta.
!>
!> @param[in]     LDU      (input) INTEGER @n
!>                         The leading dimension of the array U.
!>
!> @param[in]     SC       (input) INTEGER @n
!>                         blocking size.
!>
!> @param[in]     INDX     (input) INTEGER array, dimension (N)                     @n
!>                         The permutation used to sort the contents of DLAMDA into @n
!>                         ascending order.
!>
!> @param[in]     CTOT     (input) INTEGER array, dimension (NPCOL, 4) @n
!>                         The number of COLTYP of each process column.
!>
!> @param         Q2BUF1   (workspace) DOUBLE PRECISION array, dimension (LQ2BUF)
!>
!> @param         Q2BUF2   (workspace) DOUBLE PRECISION array, dimension (LQ2BUF)
!>
!> @param         LQ2BUF   (input) INTEGER @n
!>                         The leading dimension of the array LQ2BUF (=NP*NQ)
!>
!> @param         Z        (workspace) DOUBLE PRECISION array, dimension (K)
!>
!> @param         BUF      (workspace) DOUBLE PRECISION array, dimension (4*K)
!>
!> @param         INDROW   (workspace) INTEGER array, dimension (N)
!>
!> @param         INDCOL   (workspace) INTEGER array, dimension (N)
!>
!> @param         INDXC    (workspace) INTEGER array, dimension (N)
!>
!> @param         INDXR    (workspace) INTEGER array, dimension (N)
!>
!> @param         INDXCB   (workspace) INTEGER array, dimension (N)
!>
!> @param[out]    INFO     (global output) INTEGER @n
!>                         = 0: successful exit    @n
!>                         /=0: error exit
!>
!> @param[out]    prof     (global output) type(FS_prof) @n
!>                         profiling information of each subroutines.
!>
!> @note This routine is modified from ScaLAPACK PDLAED3.f
!>
subroutine FS_PDLAED3( K, N, N1, D, RHO, DLAMDA, W, Q, LDQ, SUBTREE, &
                       Q2, LDQ2, U, LDU, SC, INDX, CTOT, &
                       Q2BUF1, Q2BUF2, LQ2BUF, Z, BUF, INDROW, &
                       INDCOL, INDXC, INDXR, INDXCB, &
                       INFO, prof )
      USE EIGEN_DC_MOD
      use FS_const_mod
      use FS_libs_mod
      use FS_prof_mod
      use FS_dividing_mod
      use FS_MPI_Group
!$    use omp_lib
      implicit none
!
!     .. Scalar Arguments ..
      integer          , intent(in)    :: K
      integer          , intent(in)    :: N, N1, LDQ, LDQ2, LDU, LQ2BUF
      integer          , intent(in)    :: SC !blocking size of calculate U
      real(kind(0.0d0)), intent(inout) :: RHO
      type(bt_node)    , intent(in)    :: SUBTREE
      integer          , intent(out)   :: INFO
      type(FS_prof)    , intent(inout) :: prof
!
!     .. Array Arguments ..
      real(kind(0.0d0)), intent(out)   :: D(N)
      real(kind(0.0d0)), intent(inout) :: DLAMDA(K)
      real(kind(0.0d0)), intent(in)    :: W(K)
      real(kind(0.0d0)), intent(out)   :: Q(LDQ,*)
      real(kind(0.0d0)), intent(inout) :: Q2(LDQ2,*)
      real(kind(0.0d0)), intent(out)   :: U(LDU,*)
      integer          , intent(in)    :: INDX(K)
      integer          , intent(in)    :: CTOT(0:SUBTREE%y_nnod-1,4)
      ! work
      real(kind(0.0d0)), target :: Q2BUF1(LQ2BUF)
      real(kind(0.0d0)), target :: Q2BUF2(LQ2BUF)
      real(kind(0.0d0))         :: Z(K)
      real(kind(0.0d0)), target :: BUF(4*K)

      integer :: INDROW(N), INDCOL(N), INDXC(K), INDXR(K), INDXCB(K)
!     ..
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
      integer :: I, J, M, JJ, KK, ierr
      integer :: NBLK, NB, NP, NQ
      integer :: NPROW, NPCOL, MYROW, MYCOL
      integer :: MYKL, MYKLR, MYKLB, KL, KLR, COL, PDC, ROW, PDR
      integer :: SINFO, IINFO
      integer :: KLC, GI
      integer :: JU, JJU, IU, IIU
      integer :: PJ, PJCOL, NP1, NP2, N12, N23
      integer :: MINROW, MAXROW
      integer :: IQ, JQ, IQ2, JQ2, JS, PJCOLN
      integer :: DSTCOL, SRCCOL, NSEND, NRECV, req(2), stat(MPI_STATUS_SIZE)
      real(kind(0.0d0)) :: AUX, TEMP
      real(kind(0.0d0)), pointer :: SENDQ2(:), RECVQ2(:)
      real(kind(0.0d0)) :: SBUFD_MIN
      real(kind(0.0d0)) :: SBUFB_MIN
      integer :: NPA
#ifdef _MPITEST
      logical :: flag
#endif
      integer Kmax
!     ..
!     .. Local Arrays ..
      real(kind(0.0d0)), allocatable :: SBUF(:) !thread private
      real(kind(0.0d0)), allocatable :: SZ(:) !thread private
      real(kind(0.0d0)), pointer :: SBUFD(:) !minimum DELTA
      real(kind(0.0d0)), pointer :: SBUFB(:) !old eigenvalue at minimum DELTA
!     ..
!     .. External Functions ..
      real(kind(0.0d0)) :: DLAMC3, DNRM2
      external DLAMC3, DNRM2
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
      INFO = 0
#ifdef _DEBUGLOG
!      CALL MPI_Barrier(SUBTREE%MERGE_COMM,ierr)
      if( FS_MYRANK.eq.0 ) then
        write(*,'(a)') "FS_PDLAED3 start."
      endif
#endif
#if TIMER_PRINT
      CALL FS_prof_start(prof, 60)
#endif
!
!     Quick return if possible
      IF( K.EQ.0 ) THEN
        GO TO 190
      END IF
!
      CALL FS_GRIDINFO(SUBTREE, NPROW, NPCOL, MYROW, MYCOL)
      NBLK = FS_get_NBLK(SUBTREE)
      NB   = FS_get_NB(SUBTREE)
      NP   = (NBLK / NPROW) * NB
      NQ   = (NBLK / NPCOL) * NB
!
      ROW = 0
      COL = 0
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
      MYKL = CTOT( MYCOL, 1 ) + CTOT( MYCOL, 2 ) + CTOT( MYCOL, 3 )
      KLR = MYKL / NPROW
      IF( MYROW.EQ.0 ) THEN
         MYKLR = KLR + MOD( MYKL, NPROW )
      ELSE
         MYKLR = KLR
      END IF
!
      PDC = 1
      COL = 0
      DO
         IF( MYCOL == COL ) EXIT
         PDC = PDC + CTOT( COL, 1 ) + CTOT( COL, 2 ) + CTOT( COL, 3 )
         COL = MOD( COL+1, NPCOL )
      END DO
!
      PDR = PDC
      KL = KLR + MOD( MYKL, NPROW )
      ROW = 0
      DO
         IF( MYROW == ROW ) EXIT
         PDR = PDR + KL
         KL = KLR
         ROW = MOD( ROW+1, NPROW )
      END DO
!
!$OMP PARALLEL DO PRIVATE(I)
      DO I = 1, K
         DLAMDA( I ) = DLAMC3( DLAMDA( I ), DLAMDA( I ) ) - DLAMDA( I )
      END DO
!$OMP END PARALLEL DO
!     
!---  
!     
!$OMP PARALLEL DO PRIVATE(I)
      DO I = 1, 4*K
         BUF(I) = ZERO
      END DO
!$OMP END PARALLEL DO
!
#if TIMER_PRINT>1
      CALL FS_prof_start(prof, 61)
#endif
      SINFO = INFO
      IF( MYKLR.GT.0 ) THEN
!
!$OMP PARALLEL PRIVATE(I,J,KK,SBUF,AUX,IINFO,TEMP,SZ,SBUFD_MIN,SBUFB_MIN) REDUCTION(MAX:SINFO)
         allocate( SZ( 1:K ) )
         SZ( 1:K ) = ONE     
         allocate( SBUF( 1:K ) )
!$OMP DO SCHEDULE(STATIC,1)
         DO I = 1, MYKLR
            KK = PDR + I - 1
            CALL DLAED4( K, KK, DLAMDA, W, SBUF, RHO, AUX, IINFO )
            IF( K == 1 .OR. K == 2 ) THEN
               BUF(     (PDR-PDC)+I) = ZERO
               BUF(MYKL+(PDR-PDC)+I) = AUX
            ELSE
               ! minumum DELTA
               SBUFD_MIN =   SBUF(KK)
               SBUFB_MIN = DLAMDA(KK)
               if( KK.lt.K ) then
                 if( DABS(SBUF(KK+1)) .lt. DABS(SBUFD_MIN) ) then
                   SBUFD_MIN =   SBUF(KK+1)
                   SBUFB_MIN = DLAMDA(KK+1)
                 endif
               endif
               BUF(     (PDR-PDC)+I) = SBUFD_MIN
               BUF(MYKL+(PDR-PDC)+I) = SBUFB_MIN
            END IF
            IF( IINFO.NE.0 ) THEN
               SINFO = KK
            END IF
!     
!     ..Compute part of z
!     
!OCL NOFP_RELAXED
!OCL NOFP_CONTRACT
!OCL NOEVAL
            DO J=1,K
               TEMP = DLAMDA( J )-DLAMDA( KK )
               IF ( J == KK ) then
                 TEMP = ONE
               else
                 TEMP = TEMP
               endif
               SBUF( J ) = SBUF( J ) / TEMP
            ENDDO
            DO J=1,K
               SZ( J ) = SZ( J ) * SBUF( J )
            ENDDO

         END DO
!$OMP END DO
         deallocate( SBUF )
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
         deallocate( SZ )
!$OMP END PARALLEL
         INFO = SINFO

      ELSE

!$OMP PARALLEL DO PRIVATE(I)
         DO I=1, K
            Z( I ) = ONE
         END DO
!$OMP END PARALLEL DO

      END IF
#if TIMER_PRINT>1
      CALL FS_prof_end(prof, 61)
#endif
!
#if TIMER_PRINT>1
      CALL FS_prof_start(prof, 62)
#endif
!
!$OMP PARALLEL DO PRIVATE(I)
      DO I = 1, K
         BUF(2*K+I)=Z(I)
      END DO
!$OMP END PARALLEL DO
      call MPI_Group_Allreduce(&
           BUF(2*K+1:),&
           Z,&
           K,&
           MPI_DOUBLE_PRECISION,&
           MPI_PROD,&
           FS_COMM_WORLD,&
           SUBTREE%MERGE_GROUP,&
           ierr)

!$OMP PARALLEL DO PRIVATE(I)
      DO I = 1, K
         Z( I ) = SIGN( SQRT( -Z( I ) ), W( I ) )
      END DO
!$OMP END PARALLEL DO
!
      IF ( MYKL > 0 ) THEN
      call MPI_Group_Allreduce(&
           BUF(1:),&
           BUF(2*K+1:),&
           2*MYKL,&
           MPI_DOUBLE_PRECISION,&
           MPI_SUM,&
           FS_COMM_WORLD,&
           SUBTREE%MERGE_GROUP_X,&
           ierr)
      END IF
!
!---
!
!$OMP PARALLEL PRIVATE(I)
!$OMP DO
      DO I = 1, 2*MYKL
         BUF(I) = ZERO
      END DO
!$OMP END DO
!$OMP DO
      DO I = 1, MYKL
         BUF(  PDC+I-1) = BUF(2*K     +I)
         BUF(K+PDC+I-1) = BUF(2*K+MYKL+I)
      END DO
!$OMP END DO
!$OMP END PARALLEL
!
      call MPI_Group_Allreduce(&
           BUF(1:),&
           BUF(2*K+1:),&
           2*K,&
           MPI_DOUBLE_PRECISION,&
           MPI_SUM,&
           FS_COMM_WORLD,&
           SUBTREE%MERGE_GROUP_Y,&
           ierr)
!
!     Copy of D at the good place
!
      SBUFD => BUF( 2*K+1: )
      SBUFB => BUF( 3*K+1: )
      IF( K == 1 .OR. K == 2 ) THEN
         DO I = 1, K
            GI = INDX( I )
            D( GI ) = SBUFB( I )
         END DO
      ELSE
!$OMP PARALLEL DO PRIVATE(I,GI)
         DO I = 1, K
            GI = INDX( I )
            D( GI ) = SBUFB( I ) - SBUFD( I )
         END DO
!$OMP END PARALLEL DO
      END IF
!
#if TIMER_PRINT>1
      CALL FS_prof_end(prof, 62)
#endif
!
      ! 自身のCOLインデクスリストを作成
      KLC = 0
      DO I = 1, K
         GI = INDX( I )
         COL = INDCOL( GI )
         IF( COL.EQ.MYCOL ) THEN
            KLC = KLC + 1
            INDXC( KLC ) = I
         END IF
      END DO

#ifdef _BLOCKING_DGEMM
      ! ブロッキングのために列方向昇順の逆引きリストを作成
!$OMP PARALLEL DO PRIVATE(J, KK, JU, JJU)
      DO J = 1, MYKL
         KK = INDXC( J )
         JU = INDX( KK )
         JJU = FS_INDXG2L( 'C', JU, SUBTREE )
         INDXCB(JJU) = J
      END DO
!$OMP END PARALLEL DO
#endif
!
!     担当プロセス列の最大/最小インデクスを取得
!     N1, N2のどちらを保持しているかを判定するため
!     NP1 : 上側の担当次数
!     NP2 : 下側の担当次数
!
      MINROW = N
      MAXROW = 1
      NPA = 0
!$OMP PARALLEL DO PRIVATE(I) REDUCTION(MIN:MINROW) REDUCTION(MAX:MAXROW) REDUCTION(+:NPA)
      DO I = 1, N
         IF( INDROW(I).EQ.MYROW ) THEN
            MINROW = MIN(MINROW,I)
            MAXROW = MAX(MAXROW,I)
            NPA = NPA + 1 !アクティブなNP数(次数が割り切れないときに対応)
         ENDIF
      END DO
!$OMP END PARALLEL DO
      IF( MINROW.gt.N1 ) THEN
        ! 上側を持たない
        NP1 = 0
        NP2 = NPA
      ELSE IF( MAXROW.le.N1 ) THEN
        ! 下側を持たない
        NP1 = NPA
        NP2 = 0
      ELSE
        ! 両側を持つ
        NP1 = NP/2
        NP2 = NPA-NP/2
      ENDIF
!
!     Qの初期化
!     デフレーションされた列はQ2からコピー
!
#if TIMER_PRINT>1
      CALL FS_prof_start(prof, 63)
#endif
!$OMP PARALLEL PRIVATE(I, J)
!$OMP DO COLLAPSE(2)
      DO J = 1, MYKL
        DO I = 1, NP
          Q(I, J) = ZERO
        END DO
      END DO
!$OMP END DO
!$OMP DO COLLAPSE(2)
      DO J = MYKL+1, NQ
        DO I = 1, NP
          Q(I, J) = Q2(I, J)
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
#if TIMER_PRINT>1
      CALL FS_prof_end(prof, 63)
#endif

!
!     Compute eigenvectors of the modified rank-1 modification.
!
      DO PJ = 0, NPCOL-1

         ! 送受信バッファのポインタ
         IF( MOD(PJ,2).eq.0 ) THEN
            SENDQ2 => Q2BUF1
            RECVQ2 => Q2BUF2
         ELSE
            SENDQ2 => Q2BUF2
            RECVQ2 => Q2BUF1
         ENDIF

         ! 演算対象のプロセス列
         PJCOL = MOD( MYCOL+NPCOL+PJ, NPCOL )

         ! 処理する列数
         MYKL = CTOT( MYCOL, 1 ) + CTOT( MYCOL, 2 ) + CTOT( MYCOL, 3 )
         N12  = CTOT( PJCOL, 1 ) + CTOT( PJCOL, 2 )
         N23  = CTOT( PJCOL, 2 ) + CTOT( PJCOL, 3 )

        ! PJCOLのROWインデクスリストを作成
         KLR = 0
         DO I = 1, K
            GI = INDX( I )
            ROW = INDCOL( GI )
            IF( ROW.EQ.PJCOL ) THEN
               KLR = KLR + 1
               INDXR( KLR ) = I
            END IF
         END DO

         ! 前ループの送受信のwaitと展開
         ! wait -> RECVBUF -> Q2
         IF( PJ .ne. 0 .and. NPCOL .gt. 1 ) THEN

#if TIMER_PRINT>1
            CALL FS_prof_start(prof, 64)
#endif
            ! wait for irecv
            IF( NRECV.gt.0 ) THEN
              call MPI_WAIT(req(2), stat, ierr)
            ENDIF

            ! copy RECVBUF -> Q2
            ! 上側
            IF( NP1.gt.0 ) THEN
!$OMP PARALLEL DO PRIVATE(J,JQ2,JS)
               DO J = 1, N12
                  JQ2 = J
                  JS = (J-1) * NP1 + 1
                  CALL DCOPY( NP1, RECVQ2(JS), 1, Q2( 1, JQ2 ), 1 )
               END DO
!$OMP END PARALLEL DO
            END IF
            ! 下側
            IF( NP2.gt.0 ) THEN
!$OMP PARALLEL DO PRIVATE(J,JQ2,JS)
               DO J = 1, N23
                  JQ2 = J + CTOT( PJCOL, 1 )
                  JS = N12 * NP1 + (J-1) * NP2 + 1
                  CALL DCOPY( NP2, RECVQ2(JS), 1, Q2( NP1+1, JQ2 ), 1 )
               END DO
!$OMP END PARALLEL DO
            ENDIF

            ! wait for isend
            IF( NSEND.gt.0 ) THEN
              call MPI_WAIT(req(1), stat, ierr)
            ENDIF

#if TIMER_PRINT>1
            CALL FS_prof_end(prof, 64)
#endif
         ENDIF

         ! 送受信バッファのポインタの入れ替え
         IF( MOD(PJ,2).eq.0 ) THEN
            SENDQ2 => Q2BUF2
            RECVQ2 => Q2BUF1
         ELSE
            SENDQ2 => Q2BUF1
            RECVQ2 => Q2BUF2
         ENDIF

         ! 次ループ用の送受信
         ! Q2 -> SENDBUF -> isend
         IF( PJ .ne. NPCOL-1  .and. NPCOL .gt. 1 ) THEN

#if TIMER_PRINT>1
            CALL FS_prof_start(prof, 65)
#endif
            ! copy Q2 -> SENDBUF
            ! RECVBUFとSENDBUFを切り替えながら処理するためSENDBUFへの格納は初回のみ
            IF( PJ.eq.0 ) THEN
               ! 上側
               IF( NP1.gt.0 ) THEN
!$OMP PARALLEL DO PRIVATE(J,JQ2,JS)
                  DO J = 1, N12
                     JQ2 = J
                     JS = (J-1) * NP1 + 1
                     CALL DCOPY( NP1, Q2( 1, JQ2 ), 1, SENDQ2(JS), 1 )
                  END DO
!$OMP END PARALLEL DO
               END IF
               ! 下側
               IF( NP2.gt.0 ) THEN
!$OMP PARALLEL DO PRIVATE(J,JQ2,JS)
                  DO J = 1, N23
                     JQ2 = J + CTOT( PJCOL, 1 )
                     JS = N12 * NP1 + (J-1) * NP2 + 1
                     CALL DCOPY( NP2, Q2( NP1+1, JQ2 ), 1, SENDQ2(JS), 1 )
                  END DO
!$OMP END PARALLEL DO
               ENDIF
            ENDIF

            ! isend
            DSTCOL = MOD(MYCOL+NPCOL-1, NPCOL) !左側に送信
            NSEND  = NP1 * N12 + NP2 * N23
            IF( NSEND.gt.0 ) THEN
              DSTCOL=SUBTREE%group_Y_processranklist(DSTCOL+1)
              call MPI_ISEND(SENDQ2, NSEND, MPI_DOUBLE_PRECISION, &
                             DSTCOL, 1, FS_COMM_WORLD, req(1), &
                             ierr)

            ENDIF

            ! irecv
            SRCCOL = MOD(MYCOL+NPCOL+1, NPCOL) !右側から受信
            PJCOLN = MOD(MYCOL+NPCOL+PJ+1, NPCOL )
            NRECV  = NP1 * (CTOT(PJCOLN, 1) + CTOT(PJCOLN, 2)) &
                   + NP2 * (CTOT(PJCOLN, 2) + CTOT(PJCOLN, 3))
            IF( NRECV.gt.0 ) THEN
              SRCCOL=SUBTREE%group_Y_processranklist(SRCCOL+1)
              call MPI_IRECV(RECVQ2, NRECV, MPI_DOUBLE_PRECISION, SRCCOL, &
                             1, FS_COMM_WORLD, req(2), ierr)
            ENDIF

#ifdef _MPITEST
            IF( NSEND.gt.0 ) THEN
              call MPI_TEST(req(1), flag, stat, ierr)
            ENDIF
            IF( NRECV.gt.0 ) THEN
              call MPI_TEST(req(2), flag, stat, ierr)
            ENDIF
#endif

#if TIMER_PRINT>1
            CALL FS_prof_end(prof, 65)
#endif
         ENDIF

         ! deltaの計算と行列積
         IF( MYKL.NE.0 ) THEN

#if TIMER_PRINT>1
            CALL FS_prof_start(prof, 66)
#endif

            ! ブロッキングあり
#ifdef _BLOCKING_DGEMM

!$OMP PARALLEL PRIVATE(M,MYKLB,I,J,JJ,KK,IU,IIU,JU,JJU,SBUF,AUX,IINFO,TEMP)
            allocate( SBUF( 1:K ) )
            DO M = 1, MYKL, SC

               MYKLB = MIN(M+SC-1, MYKL) - M + 1
!$OMP DO SCHEDULE(DYNAMIC)
               DO JJ = 1, MYKLB
                  J = INDXCB(M + JJ - 1)
                  KK = INDXC( J )
                  JU = INDX( KK )
                  JJU = FS_INDXG2L( 'C', JU, SUBTREE )

                  IF( K == 1 .OR. K == 2 ) THEN
                     CALL DLAED4( K, KK, DLAMDA, W, SBUF, RHO, AUX, IINFO )
                  ELSE
                     SBUF( 1:K ) = DLAMDA( 1:K ) - SBUFB(KK)
                     CALL ADD_SCALAR(SBUF, K, SBUFD(KK))
                  END IF

                  IF( K == 1 .OR. K == 2 ) THEN
                     DO I = 1, KLR
                        KK = INDXR( I )
                        IU = INDX( KK )
                        IIU = FS_INDXG2L( 'C', IU, SUBTREE )
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
                     IIU = FS_INDXG2L( 'C', IU, SUBTREE )
                     U( IIU, JJU ) = SBUF( KK ) / TEMP
                  END DO
               END DO
!$OMP END DO
            deallocate( SBUF )
!$OMP MASTER

               ! Compute the updated eigenvectors.

               ! 上側のDGEMM
               IF( NP1.NE.0 .and. N12.NE.0 ) THEN

                  IIU = 1
                  JJU = M
                  IQ = 1
                  JQ = M

#if TIMER_PRINT>1
                  CALL FS_prof_start(prof, 67)
#endif
                  CALL DGEMM( 'N', 'N', NP1, MYKLB, N12, ONE, Q2, LDQ2, &
                              U(IIU, JJU), LDU, ONE, Q(IQ, JQ), LDQ )
                  flops = flops + 2*DBLE(NP1)*DBLE(MYKLB)*DBLE(N12)
#if TIMER_PRINT>1
                  CALL FS_prof_end(prof, 67)
#endif

               ENDIF

               ! 下側のDGEMM
               IF( NP2.NE.0 .and. N23.NE.0 ) THEN

                  IQ2 = NP1 + 1
                  JQ2 = CTOT(PJCOL,1) + 1
                  IIU = CTOT(PJCOL,1) + 1
                  JJU = M
                  IQ = NP1 + 1
                  JQ = M

#if TIMER_PRINT>1
                  CALL FS_prof_start(prof, 67)
#endif
                  CALL DGEMM( 'N', 'N', NP2, MYKLB, N23, ONE, Q2(IQ2, JQ2), &
                              LDQ2, U(IIU, JJU), LDU, ONE, Q(IQ, JQ), LDQ )
                  flops = flops + 2*DBLE(NP2)*DBLE(MYKLB)*DBLE(N23)
#if TIMER_PRINT>1
                  CALL FS_prof_end(prof, 67)
#endif


               ENDIF
!$OMP END MASTER

            END DO ! DO M

!$OMP END PARALLEL

#endif
            ! ブロッキングあり終了

            ! ブロッキングなし
#ifndef _BLOCKING_DGEMM

!$OMP PARALLEL PRIVATE(I,J,KK,IU,IIU,JU,JJU,SBUF,AUX,IINFO,TEMP)
            allocate( SBUF( 1:K ) )
!$OMP DO
            DO J = 1, MYKL
               KK = INDXC( J )
               JU = INDX( KK )
               JJU = FS_INDXG2L( 'C', JU, SUBTREE )

               IF( K == 1 .OR. K == 2 ) THEN
                  CALL DLAED4( K, KK, DLAMDA, W, SBUF, RHO, AUX, IINFO )
               ELSE
                  SBUF( 1:K ) = DLAMDA( 1:K ) - SBUFB(KK)
                  CALL ADD_SCALAR(SBUF, K, SBUFD(KK))
               END IF

               IF( K == 1 .OR. K == 2 ) THEN
                  DO I = 1, KLR
                     KK = INDXR( I )
                     IU = INDX( KK )
                     IIU = FS_INDXG2L( 'C', IU, SUBTREE )
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
                  IIU = FS_INDXG2L( 'C', IU, SUBTREE )
                  U( IIU, JJU ) = SBUF( KK ) / TEMP
               END DO
            END DO
!$OMP END DO
            deallocate( SBUF )
!$OMP END PARALLEL

            ! Compute the updated eigenvectors.

            ! 上側のDGEMM
            IF( NP1.NE.0 .and. N12.NE.0 ) THEN

#if TIMER_PRINT>1
              CALL FS_prof_start(prof, 67)
#endif
              CALL DGEMM( 'N', 'N', NP1, MYKL, N12, ONE, Q2, LDQ2, U, &
                          LDU, ONE, Q, LDQ )
              flops = flops + 2*DBLE(NP1)*DBLE(MYKL)*DBLE(N12)
#if TIMER_PRINT>1
              CALL FS_prof_end(prof, 67)
#endif

            ENDIF

            ! 下側のDGEMM
            IF( NP2.NE.0 .and. N23.NE.0 ) THEN

              IQ2 = NP1 + 1
              JQ2 = CTOT(PJCOL,1) + 1
              IIU = CTOT(PJCOL,1) + 1
              JJU = 1
              IQ = NP1 + 1
              JQ = 1

#if TIMER_PRINT>1
              CALL FS_prof_start(prof, 67)
#endif
              CALL DGEMM( 'N', 'N', NP2, MYKL, N23, ONE, Q2(IQ2, JQ2), &
                          LDQ2, U(IIU, JJU), LDU, ONE, Q(IQ, JQ), LDQ )
              flops = flops + 2*DBLE(NP2)*DBLE(MYKL)*DBLE(N23)
#if TIMER_PRINT>1
              CALL FS_prof_end(prof, 67)
#endif

           ENDIF

#endif
            ! ブロッキングなし終了

#if TIMER_PRINT>1
            CALL FS_prof_end(prof, 66)
#endif

         END IF !MYKL.NE.0

      END DO !PJ
!-
      SBUFD => null()
      SBUFB => null()
!-
 190  CONTINUE

#if TIMER_PRINT
      CALL FS_prof_end(prof, 60)
#endif
#ifdef _DEBUGLOG
!      CALL MPI_Barrier(SUBTREE%MERGE_COMM,ierr)
      if( FS_MYRANK.eq.0 ) then
        write(*,'(a,i0)') "FS_PDLAED3 end. INFO=", INFO
      endif
#endif
      RETURN
!
!     End of FS_PDLAED3
!     

      CONTAINS

!> subroutine ADD_SCALAR
!>
!> @brief @n
!>  Purpose @n
!>  ======= @n
!>  Add scalar value to vector array.
!
!  Arguments
!  =========
!>
!> @param[inout]  V        (input/output) DOUBLE PRECISION array @n
!>                          vector array.
!>
!> @param[in]     NV       (input) INTEGER @n
!>                          vector array size.
!>
!> @param[in]     S        (input) DOUBLE PRECISION @n
!>                          scalar value.
!>
subroutine ADD_SCALAR(V, NV, S)
      implicit none
      integer :: NV
      real(kind(0.0d0)) :: V(NV), S

      V(:) = V(:) + S
end subroutine ADD_SCALAR

end subroutine FS_PDLAED3

