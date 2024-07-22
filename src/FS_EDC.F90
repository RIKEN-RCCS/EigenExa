!>
!> @file   FS_EDC.F90
!> @brief  subroutine FS_EDC
!!>
!
!
!
MODULE FS_EDC_MOD
  use FS_const_mod
  use FS_libs_mod
  use FS_prof_mod

  use eigen_libs_mod
  use, intrinsic :: ieee_arithmetic

  implicit none
CONTAINS
!
!> subroutine FS_EDC
!> @brief @n
!>   Purpose @n
!>   ======= @n
!>   FS_EDC computes all eigenvalues and eigenvectors of a             @n
!>   symmetric tridiagonal matrix in parallel, using the divide and    @n
!>   conquer algorithm.
!
!  Arguments
!  =========
!>
!> @param[in]     N      (global input) INTEGER @n
!>                       The order of the tridiagonal matrix T.  N >= 0.
!>
!> @param[in,out] D      (global input/output) DOUBLE PRECISION array, dimension (N) @n
!>                       On entry, the diagonal elements of the tridiagonal matrix.  @n
!>                       On exit, if INFO = 0, the eigenvalues in descending order.
!>
!> @param[in,out] E      (global input/output) DOUBLE PRECISION array, dimension (N-1) @n
!>                       On entry, the subdiagonal elements of the tridiagonal matrix. @n
!>                       On exit, E has been destroyed.
!>
!> @param[in,out] Q      (local output) DOUBLE PRECISION array, local dimension (LDQ, *)   @n
!>                       Q contains the orthonormal eigenvectors of the symmetric tridiagonal matrix.   @n
!>                       On output, Q is distributed across the P processes in non block cyclic format. @n
!>
!> @param[in]     LDQ    (local input) INTEGER @n
!>                       leading dimension of array Q.
!>
!> @param         WORK   (local workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!>
!> @param[in]     LWORK  (local input/output) INTEGER, the dimension of the array WORK. @n
!>                       LWORK = 1 + 6*N + 3*NP*NQ + NQ*NQ @n
!>                       LWORK can be obtained from subroutine FS_WorkSize.
!>
!> @param         IWORK  (local workspace/output) INTEGER array, dimension (LIWORK)
!>
!> @param[in]     LIWORK (input) INTEGER                   @n
!>                       The dimension of the array IWORK. @n
!>                       LIWORK = 1 + 8*N + 8*NPCOL        @n
!>                       LIWORK can be obtained from subroutine FS_WorkSize.
!>
!> @param[out]    info   (global output) INTEGER @n
!>                       = 0: successful exit   @n
!>                       /=0: error exit
!>
!> @param[out]    prof   (global output) type(FS_prof) @n
!>                       profiling information of each subroutines.
!>
!> @note This routine is modified from ScaLAPACK PDSTEDC.f
!>
subroutine FS_EDC(N, D, E, Q, LDQ, WORK, LWORK, IWORK, LIWORK, INFO, prof)
      implicit none

      integer          , intent(in)    :: N, LDQ
      integer(8)       , intent(in)    :: LWORK, LIWORK
      real(kind(0.0d0)), intent(inout) :: D(N)
      real(kind(0.0d0)), intent(inout) :: E(N-1)
      real(kind(0.0d0)), intent(out)   :: Q(LDQ,*)
      integer(8)       , intent(out)   :: INFO
      type(FS_prof)    , intent(inout), optional :: prof
      ! work
      real(kind(0.0d0)) :: WORK(LWORK)
      integer           :: IWORK(LIWORK)

!     .. Parameters ..

!     .. Local Scalars ..
      integer          :: nnod, x_nnod, y_nnod
      integer          :: inod, x_inod, y_inod

      integer          :: FS_nnod, FS_x_nnod, FS_y_nnod
      integer          :: FS_inod, FS_x_inod, FS_y_inod

      integer          :: INFO_I4, LWORK_I4, LIWORK_I4

      DOUBLE PRECISION :: ORGNRM
      integer :: I
      type(FS_prof)    :: prof_tmp

!     .. External Functions ..
      DOUBLE PRECISION :: DLANST
!
!     .. Executable Statements ..
!
      INFO = 0
#ifdef _DEBUGLOG
      if( FS_MYRANK.eq.0 ) then
        write(*,'(a)') "FS_EDC start."
      endif
#endif
!
#if TIMER_PRINT
      if( PRESENT(prof) ) then
        prof_tmp = prof
      else
        call FS_prof_init(prof_tmp)
      endif
      call FS_prof_start(prof_tmp, 10)
#endif
!
      call FS_get_procs(FS_nnod,FS_x_nnod,FS_y_nnod)
      call FS_get_id(FS_inod,FS_x_inod,FS_y_inod)
      call eigen_get_procs(nnod,x_nnod,y_nnod)
      call eigen_get_id(inod,x_inod,y_inod)

!     
!     Quick return
!     
      IF( N.EQ.0 ) THEN
        GO TO 10
      END IF
      IF( N.EQ.1 ) THEN
        IF( (x_inod.EQ.1) .AND. (y_inod.EQ.1) ) THEN
          Q(1,1) = ONE
        END IF
        GO TO 10
      END IF
!     
!     If P=NPROW*NPCOL=1, solve the problem with DSTEDC.
!     
      IF( x_nnod*y_nnod.EQ.1 ) THEN
#if TIMER_PRINT
        CALL FS_prof_start(prof_tmp, 11)
#endif
        LWORK_I4 = LWORK
        LIWORK_I4 = LIWORK
        CALL DSTEDC( 'I', N, D, E, Q, LDQ, WORK, LWORK_I4, IWORK, LIWORK_I4, INFO_I4 )
        INFO = INFO_I4
#if TIMER_PRINT
        CALL FS_prof_end(prof_tmp, 11)
#endif
        GO TO 10
      END IF

!     
!     Scale matrix to allowable range, if necessary.
!     
      ORGNRM = DLANST( 'M', N, D, E )
      if ( ieee_is_nan( ORGNRM ) ) ORGNRM = ZERO
      IF( ORGNRM.NE.ZERO ) THEN
        CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO_I4 )
        IF( N-1 >= 1 ) THEN
          CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N-1, 1, E, N-1, INFO_I4 )
        END IF
        INFO = INFO_I4
      END IF
!
      call FS_PDLAED0(N, D, E, Q, LDQ, WORK, LWORK, IWORK, LIWORK, INFO_I4, prof_tmp)
      INFO = INFO_I4
      IF( INFO.NE.0 ) THEN
        GO TO 10
      END IF
!     
!     Scale back.
!     
      IF( ORGNRM.NE.ZERO ) THEN
        CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO_I4 )
        INFO = INFO_I4
      END IF
!     
 10   CONTINUE
!
#ifdef _DEBUGLOG
      if( FS_MYRANK.eq.0 ) then
        write(*,'(a,i0)') "FS_EDC end. INFO=", INFO
      endif
#endif
#if TIMER_PRINT
      call FS_prof_end(prof_tmp, 10)
      if( PRESENT(prof) ) then
        prof = prof_tmp
      else
        call FS_prof_finalize(prof_tmp)
      endif
#endif
!
      RETURN
end subroutine FS_EDC

END MODULE FS_EDC_MOD
