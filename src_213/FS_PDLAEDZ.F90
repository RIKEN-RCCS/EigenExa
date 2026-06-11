!>
!> @file   FS_PDLAEDZ.F90
!> @brief  subroutine FS_PDLAEDZ
!!>
!
!
!
!> subroutine FS_PDLAEDZ
!>
!> @brief  @n
!> Purpose @n
!> ======= @n
!> FS_PDLAEDZ Form the z-vector which consists of the last row of Q_1 @n
!> and the first row of Q_2.
!
!  Arguments
!  =========
!>
!> @param[in]     N        (global input) INTEGER @n
!>                         The order of the tridiagonal matrix T.  N >= 0.
!>
!> @param[in]     N1       (input) INTEGER @n
!>                         The location of the last eigenvalue in the leading sub-matrix. @n
!>                         min(1,N) <= N1 <= N.
!>
!> @param[in]     Q        (local output) DOUBLE PRECISION array,                    @n
!>                         global dimension (N, N),                                  @n
!>                         local dimension (LDQ, NQ)                                 @n
!>                         Q  contains the orthonormal eigenvectors of the symmetric @n
!>                         tridiagonal matrix.
!>
!> @param[in]     LDQ      (local input) INTEGER @n
!>                         The leading dimension of the array Q.  LDQ >= max(1,NP).
!>
!> @param[in]     SUBTREE  (input) type(bt_node) @n
!>                         sub-tree information of merge block.
!>
!> @param[out]    Z        (input) DOUBLE PRECISION array, dimension (N)                   @n
!>                         The updating vector before MPI_ALLREDUCE (the last row of the   @n
!>                         first sub-eigenvector matrix and the first row of the second    @n
!>                         sub-eigenvector matrix).
!>
!> @param[out]    prof     (global output) type(FS_prof) @n
!>                         profiling information of each subroutines.
!>
!> @note This routine is modified from ScaLAPACK PDLAEDZ.f
!>
subroutine FS_PDLAEDZ(N, N1, Q, LDQ, SUBTREE, Z, prof)
      use FS_const_mod
      use FS_libs_mod
      use FS_prof_mod
      use FS_dividing_mod
      implicit none
!
!     .. Scalar Arguments ..
      integer          , intent(in)    :: N, N1, LDQ
      type(bt_node)    , intent(in)    :: SUBTREE
      type(FS_prof)    , intent(inout) :: prof
!
!     .. Array Arguments ..
      real(kind(0.0d0)), intent(in)    :: Q(LDQ,*)
      real(kind(0.0d0)), intent(out)   :: Z(N)
!     ..
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
      integer :: NB, NBLK, NPROW, NPCOL, MYROW, MYCOL, NB1
      integer :: IZ1, JZ1, IZ1ROW, JZ1COL
      integer :: IZ2, JZ2, IZ2ROW, JZ2COL
      integer :: J, ierr
!     ..
!     .. Local Arrays ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
#ifdef _DEBUGLOG
      if( FS_MYRANK.eq.0 ) then
        write(*,'(a)') "FS_PDLAEDZ start."
      endif
#endif
#if TIMER_PRINT
      CALL FS_prof_start(prof, 40)
#endif
!
      ! プロセス情報
      CALL FS_GRIDINFO(SUBTREE, NPROW, NPCOL, MYROW, MYCOL)
      NBLK = FS_get_NBLK(SUBTREE)
      NB   = FS_get_NB(SUBTREE)
!
      DO J = 1, N
         Z(J) = ZERO
      END DO
!
      ! Z1、Z2を含むプロセス行を取得
      CALL FS_INFOG1L('R', N1  , SUBTREE, IZ1, IZ1ROW)
      CALL FS_INFOG1L('R', N1+1, SUBTREE, IZ2, IZ2ROW)
!     
!     Form z1 which consist of the last row of Q1
!     
      IF( IZ1ROW.EQ.MYROW ) THEN

        ! z1を含むプロセス列
!$OMP PARALLEL DO PRIVATE(J,JZ1,JZ1COL,NB1)
        DO J=1,N1,NB
          CALL FS_INFOG1L('C', J, SUBTREE, JZ1, JZ1COL)
          IF( JZ1COL.EQ.MYCOL ) THEN
            NB1 = MIN(N1, J+NB-1) - J + 1
            CALL DCOPY( NB1, Q(IZ1,JZ1), LDQ, Z(J), 1 )
          END IF
        END DO
!$OMP END PARALLEL DO

      END IF
!     
!     Form z2 which consist of the first row of Q2
!     
      IF( IZ2ROW.EQ.MYROW ) THEN

        ! z2を含むプロセス列
!$OMP PARALLEL DO PRIVATE(J,JZ2,JZ2COL,NB1)
        DO J=N1+1,N,NB
          CALL FS_INFOG1L('C', J, SUBTREE, JZ2, JZ2COL)
          IF( JZ2COL.EQ.MYCOL ) THEN
            NB1 = MIN(N, J+NB-1) - J + 1
            CALL DCOPY( NB1, Q(IZ2,JZ2), LDQ, Z(J), 1 )
          END IF
        END DO
!$OMP END PARALLEL DO

      END IF
!
! ver.1.1.0 move to FS_REDUCE_ZD
!     reduce
!     CALL MPI_Allreduce( WORK, Z, N, MPI_DOUBLE_PRECISION, MPI_SUM, SUBTREE%MERGE_COMM, ierr)

 10   CONTINUE
!
#if TIMER_PRINT
      CALL FS_prof_end(prof, 40)
#endif
#ifdef _DEBUGLOG
      if( FS_MYRANK.eq.0 ) then
        write(*,'(a)') "FS_PDLAEDZ end."
      endif
#endif
!
      RETURN
!
!     End of FS_PDLAEDZ
!
end subroutine FS_PDLAEDZ

