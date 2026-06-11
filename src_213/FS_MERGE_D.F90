!>
!> @file   FS_MERGE_D.F90
!> @brief  subroutine FS_MERGE_D
!!>
!
!
!
!> subroutine FS_MERGE_D
!>
!> @brief  @n
!> Purpose @n
!> ======= @n
!> gather D of sub-matrix.
!
!  Arguments
!  =========
!>
!> @param[in]     N        (global input) INTEGER @n
!>                         The order of the tridiagonal matrix T.  N >= 0.
!>
!> @param[in]     D        (local input) DOUBLE PRECISION array, dimension(N) @n
!>                         The diagonal elements of the tridiagonal matrix.
!>
!> @param[in]     SUBTREE  (input) type(bt_node) @n
!>                         sub-tree information of merge block.
!>
!> @param[out]    DOUT     (local output) DOUBLE PRECISION array, dimension (N)
!>                         generated D before MPI_ALLREDUCE.
!>
subroutine FS_MERGE_D(N, D, SUBTREE, DOUT)
      use FS_const_mod
      use FS_dividing_mod
      implicit none
!
      integer          , intent(in)  :: N
      type(bt_node)    , intent(in)  :: SUBTREE
      real(kind(0.0d0)), intent(in)  :: D(N)
      real(kind(0.0d0)), intent(out) :: DOUT(N)

      integer :: NPROW, NPCOL, MYROW, MYCOL
      integer :: NB, NB1, J, JJ, COL, I, II, ROW
      integer :: ierr

      ! reduce用バッファのゼロクリア
      DOUT(1:N) = ZERO

      ! プロセス情報取得
      CALL FS_GRIDINFO(SUBTREE, NPROW, NPCOL, MYROW, MYCOL)
      NB = FS_get_NB(SUBTREE)

      ! 縦分割/横分割
      IF( SUBTREE%direction_horizontal ) THEN

        ! 横分割のとき

        ! 先頭列を含むプロセス列
        CALL FS_INFOG1L('C', 1, SUBTREE, JJ, COL)

        ! Dをコピー
        IF( COL.EQ.MYCOL ) THEN
          DO I=1, N, NB

            ! Iを含むプロセス行
            CALL FS_INFOG1L('R', I, SUBTREE, II, ROW)
            IF( ROW.EQ.MYROW ) THEN
              NB1 = MIN(N, I+NB-1) - I + 1
              CALL DCOPY( NB1, D(I), 1, DOUT(I), 1 )
            ENDIF

          END DO
        END IF

      ELSE

        ! 横分割のとき

        ! 先頭行を含むプロセス行
        CALL FS_INFOG1L('R', 1, SUBTREE, II, ROW)

        ! Dをコピー
        IF( ROW.EQ.MYROW ) THEN
          DO J=1, N, NB

            ! Jを含むプロセス列
            CALL FS_INFOG1L('C', J, SUBTREE, JJ, COL)
            IF( COL.EQ.MYCOL ) THEN
              NB1 = MIN(N, J+NB-1) - J + 1
              CALL DCOPY( NB1, D(J), 1, DOUT(J), 1 )
            ENDIF

          END DO
        END IF

      END IF

!      ! reduce
!      CALL MPI_Allreduce( WORK, D, N, MPI_DOUBLE_PRECISION, MPI_SUM, &
!                          SUBTREE%MERGE_COMM, ierr)

      return
end subroutine FS_MERGE_D

