!#define PROF_DETAIL
!>
!> @file   FS_REDUCE_ZD.F90
!> @brief  subroutine FS_REDUCE_ZD
!!>
!
!
!
!> subroutine FS_REDUCE_ZD
!>
!> @brief  @n
!> Purpose @n
!> ======= @n
!> MPI_ALLREDUCE Z and D
!
!  Arguments
!  =========
!>
!> @param[in]     N        (global input) INTEGER @n
!>                         The order of the tridiagonal matrix T.  N >= 0.
!>
!> @param[in]     SUBTREE  (input) type(bt_node) @n
!>                         sub-tree information of merge block.
!>
!> @param[in]     WORK     (input) DOUBLE PRECISION array, dimension (N,2)         @n
!>                         WORK(:,1) is the updating vector before MPI_ALLREDUCE.  @n
!>                         WORK(:,2) is the generated D before MPI_ALLREDUCE.
!>
!> @param[out]    Z        (local output) DOUBLE PRECISION array, dimension (N)                   @n
!>                         The updating vector (the last row of the first sub-eigenvector  @n
!>                         matrix and the first row of the second sub-eigenvector matrix).
!>
!> @param[out]    D        (local output) DOUBLE PRECISION array, dimension (N)
!>                         generated D.
!>
!> @param[out]    prof     (global output) type(FS_prof) @n
!>                         profiling information of each subroutines.
!>
subroutine FS_REDUCE_ZD( N, SUBTREE, WORK, Z, D, prof )
      use FS_const_mod
      use FS_libs_mod
      use FS_prof_mod
      use FS_dividing_mod
      use FS_MPI_Group
      implicit none
!
!     .. Scalar Arguments ..
      integer          , intent(in)    :: N
      type(bt_node)    , intent(in)    :: SUBTREE
      type(FS_prof)    , intent(inout) :: prof
!
!     .. Array Arguments ..
      ! real(kind(0.0d0)), intent(in)    :: WORK(N,2)
      ! real(kind(0.0d0)), intent(out)   :: Z(*)
      real(kind(0.0d0)),target   :: WORK(N,2)
      real(kind(0.0d0)),target   :: Z(*)
      real(kind(0.0d0)), intent(out)   :: D(*)
!     ..
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
      integer :: NB, NBLK, NPROW, NPCOL, MYROW, MYCOL, NB1
      integer :: IZ1, JZ1, IZ1ROW, JZ1COL
      integer :: IZ2, JZ2, IZ2ROW, JZ2COL
      integer :: J, ierr
      real(8) :: tmp

!     ..
!     .. Local Arrays ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
#ifdef _DEBUGLOG
      CALL MPI_Barrier(SUBTREE%MERGE_COMM,ierr)
      if( FS_MYRANK.eq.0 ) then
        write(*,'(a)') "FS_REDUCE_ZD start."
      endif
#endif
#if TIMER_PRINT
      CALL FS_prof_start(prof, 45)
#endif
!
#if TIMER_PRINT
#ifdef PROF_DETAIL
      CALL FS_prof_start(prof, 46)
      CALL MPI_Barrier(SUBTREE%MERGE_COMM,ierr)
      CALL FS_prof_end(prof, 46)
      CALL FS_prof_start(prof, 47)
#endif
#endif
!
!     reduce
#if 0
      CALL MPI_Allreduce( WORK, Z, N*2, MPI_DOUBLE_PRECISION, MPI_SUM, SUBTREE%MERGE_COMM, ierr)
#else
      call MPI_Group_Allreduce(&
           WORK(1:,1),&
           Z,&
           N*2,&
           MPI_DOUBLE_PRECISION,&
           MPI_SUM,&
           FS_COMM_WORLD,&
           SUBTREE%MERGE_GROUP,&
           ierr)
#endif          
!
!     copy D
!      D(1:N) = Z(N+1:N*2)

!$OMP PARALLEL DO PRIVATE(J)
      do J=1,N
         D(J) = Z(J+N)
      enddo
!$OMP END PARALLEL DO 

!
#if TIMER_PRINT
#ifdef PROF_DETAIL
      CALL FS_prof_end(prof, 47)
#endif
#endif
!
#if TIMER_PRINT
      CALL FS_prof_end(prof, 45)
#endif
#ifdef _DEBUGLOG
      CALL MPI_Barrier(SUBTREE%MERGE_COMM,ierr)
      if( FS_MYRANK.eq.0 ) then
        write(*,'(a)') "FS_REDUCE_ZD end."
      endif
#endif
!
      RETURN
!
!     End of FS_REDUCE_ZD
!
end subroutine FS_REDUCE_ZD

