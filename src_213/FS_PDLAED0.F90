!>
!> @file   FS_PDLAED0.F90
!> @brief  subroutine FS_PDLAED0
!!>
!
!
!
!> subroutine FS_PDLAED0
!> @brief  @n
!> Purpose @n
!> ======= @n
!> FS_PDLAED0 computes all eigenvalues and corresponding eigenvectors of a @n
!> symmetric tridiagonal matrix using the divide and conquer method.
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
!> @param[in,out] Q      (local output) DOUBLE PRECISION array,                    @n
!>                       global dimension (N, N),                                  @n
!>                       local dimension (LDQ, NQ)                                 @n
!>                       Q contains the orthonormal eigenvectors of the symmetric  @n
!>                       tridiagonal matrix.                                       @n
!>                       On output, Q is distributed across the P processes in non @n
!>                       block cyclic format.
!>
!> @param[in]     LDQ    (local input) INTEGER @n
!>                       The leading dimension of the array Q.  LDQ >= max(1,NP).
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
!> @param[out]    INFO   (global output) INTEGER @n
!>                       = 0: successful exit   @n
!>                       /=0: error exit
!>
!> @param[out]    prof   (global output) type(FS_prof) @n
!>                       profiling information of each subroutines.
!>
!> @note This routine is modified from ScaLAPACK PDLAED0.f
!>
subroutine FS_PDLAED0(N, D, E, Q, LDQ, WORK, LWORK, IWORK, LIWORK, INFO, prof)
      use FS_const_mod
      use FS_libs_mod
      use FS_prof_mod
      use FS_dividing_mod

      use eigen_libs_mod 
      implicit none
!
!     .. Scalar Arguments ..
      integer          , intent(in)    :: N, LDQ, LWORK, LIWORK
      integer          , intent(out)   :: INFO
      type(FS_prof)    , intent(inout) :: prof
!
!     .. Array Arguments ..
      real(kind(0.0d0)), intent(inout) :: D(N)
      real(kind(0.0d0)), intent(inout) :: E(N-1)
      real(kind(0.0d0)), intent(out)   :: Q(LDQ,*)
      ! work
!      real(kind(0.0d0)) :: WORK(LWORK)
      real(kind(0.0d0)) :: WORK(LWORK)
      integer           :: IWORK(LIWORK)
!     ..
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
      integer :: nnod, x_nnod, y_nnod
      integer :: inod, x_inod, y_inod
      integer :: MATSIZ, ID, N0, N1, NB
      integer :: IPQ, JPQ, IDROW, IDCOL
      real(kind(0.0d0)) :: RHO
      integer :: NBLK, NP, NQ, LDQ2
      integer :: IPQ2, ISENDQ, IRECVQ, IBUF
      integer :: INDROW, INDCOL, INDX, INDRCV
      type(bt_node) :: ROOT_NODE
      type(bt_node), pointer :: NODE => null()
      type(bt_node), pointer :: PARENT_NODE => null()
      integer :: I, J
      integer ierr
      integer eigen_comm, eigen_x_comm, eigen_y_comm
      type(FS_prof) :: prof_layer

      integer eigen_NP, eigen_NPROW, eigen_NPCOL
!     ..
!     .. Local Arrays ..
      logical, allocatable :: HINT(:)
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
!
      INFO = 0
      if(FS_COMM_MEMBER)then
! #ifdef _DEBUGLOG
!       if( FS_MYRANK.eq.0 ) then
!          write(*,*)"+++++++++++++++++++++"
!          DO I=1,N-1
!             write(*,*)"D,E,",I,D(I),E(I)
!          ENDDO
!          write(*,*)"D,E,",N,D(N)
!          write(*,*)"+++++++++++++++++++++"
!         write(*,'(a)') "FS_PDLAED0 start."
!       endif
! #endif

#if TIMER_PRINT
      CALL FS_prof_start(prof, 20)
#endif
      CALL FS_get_procs(nnod,x_nnod,y_nnod)
      CALL FS_get_id(inod,x_inod,y_inod)
!
!     .. divide tree
!
      ALLOCATE( HINT(nnod) )
      CALL FS_create_hint(HINT)
      CALL FS_dividing(N, D, E, ROOT_NODE, HINT, INFO, prof)

      DEALLOCATE( HINT )
      IF( INFO.NE.0 ) THEN
        GO TO 10
      END IF
!
!     .. clear Q
!
      NBLK = FS_get_NBLK(ROOT_NODE)
      NB   = FS_get_NB(ROOT_NODE)
      NP   = (NBLK / ROOT_NODE%x_nnod) * NB
      NQ   = (NBLK / ROOT_NODE%y_nnod) * NB
!$OMP PARALLEL DO PRIVATE(I,J) COLLAPSE(2)
      DO J = 1, NQ
      DO I = 1, NP
        Q(I, J) = ZERO
      END DO
      END DO
!$OMP END PARALLEL DO
!
!     .. get leaf node for own process
!
      NODE => null()
      CALL FS_dividing_getleaf(ROOT_NODE, NODE, INFO)
      IF( INFO.NE.0 ) THEN
        GO TO 10
      END IF
!
!     Solve each submatrix eigenproblem at the bottom of the divide and
!     conquer tree. D is the same on each process.
!
#if TIMER_PRINT
      CALL FS_prof_start(prof, 28)
#endif
      ID = NODE%nstart
      MATSIZ = FS_get_N_active(NODE)
      IF( MATSIZ.GT.0 ) THEN
         call FS_INFOG2L(ID,ID,ROOT_NODE,IPQ,JPQ,IDROW,IDCOL)
         if(ID.ne.N)then
            CALL DSTEDC( 'I', MATSIZ, D(ID), E(ID), Q(IPQ,JPQ), LDQ, WORK, LWORK, IWORK, LIWORK, INFO )
            !       CALL DSTEQR( 'I', MATSIZ, D(ID), E(ID), Q(IPQ,JPQ), LDQ, WORK, INFO )
         else
            ! ID = Nは最終行なので副対角要素はない.
            Q(IPQ,JPQ)=1.0;
            INFO=0
         endif

        IF( INFO.NE.0 ) THEN
          GO TO 10
        END IF
      ENDIF
#if TIMER_PRINT
      CALL FS_prof_end(prof, 28)
#endif
!
!     Successively merge eigensystems of adjacent submatrices
!     into eigensystem for the corresponding larger matrix.
!
 60   CONTINUE
      IF( associated(NODE%parent_node) ) THEN

        PARENT_NODE => NODE%parent_node
        ID  = PARENT_NODE%nstart
        N0  = FS_get_N_active(PARENT_NODE)
        N1  = FS_get_N_active(PARENT_NODE%sub_bt_node(1))
        IF( N0.EQ.N1 ) THEN
          NODE => NODE%parent_node
          GO TO 60
        ENDIF
        NB  = FS_get_NB(PARENT_NODE)
        RHO = E(ID+N1-1)

        call FS_get_QTOP(PARENT_NODE, IPQ, JPQ)

#ifdef _DEBUGLOG
  if( FS_MYRANK.eq.0 ) then
    write(*,'(a)') "+---------------------"
    write(*,'(a)')   "FS_PDLAED0 merge loop"
    ! write(*,'(a,i)') "  layer = ", PARENT_NODE%layer
    ! write(*,'(a,i)') "  N     = ", N0
    ! write(*,'(a,i)') "  NB    = ", NB
    write(*,*) "  layer = ", PARENT_NODE%layer
    write(*,*) "  N     = ", N0
    write(*,*) "  NB    = ", NB
  endif
#endif

#if TIMER_PRINT
        CALL FS_prof_init(prof_layer)
#endif
        CALL FS_PDLAED1( N0, N1, D(ID), Q(IPQ,JPQ), LDQ, PARENT_NODE, &
                         RHO, WORK, IWORK, INFO, prof_layer )
#if TIMER_PRINT>2
        CALL FS_prof_finalize(prof_layer, outall=0)
#endif
#if TIMER_PRINT
        CALL FS_prof_add(prof, prof_layer)
#endif
        IF( INFO.NE.0 ) THEN

          GO TO 10
        END IF

        NODE => NODE%parent_node

        GO TO 60

      END IF
#ifdef _DEBUGLOG
  if( FS_MYRANK.eq.0 ) then
    write(*,'(a)') "+---------------------"
  endif
#endif

      ! end of if(FS_COMM_MEMBER)then
      else
        nnod = 1; x_nnod = 1; y_nnod = 1
      endif

#ifdef _DEBUGLOG
  if( FS_MYRANK.eq.0 ) then
    write(*,'(a)') "START BCAST"
  endif
#endif
!
!     Sort eigenvalues and corresponding eigenvectors
!
      call eigen_get_comm(eigen_comm, eigen_x_comm, eigen_y_comm)
      call MPI_Bcast(x_nnod,1,MPI_INTEGER,0,eigen_comm,ierr)
      call MPI_Bcast(y_nnod,1,MPI_INTEGER,0,eigen_comm,ierr)
      
      NBLK = FS_get_NBLK(ROOT_NODE)
      NB   = FS_get_NB(ROOT_NODE)
      ! NP   = (NBLK / ROOT_NODE%x_nnod) * NB
      ! NQ   = (NBLK / ROOT_NODE%y_nnod) * NB
      NP   = (NBLK / x_nnod) * NB
      NQ   = (NBLK / y_nnod) * NB
      LDQ2 = NP
!
      IPQ2 = 1
      ISENDQ = IPQ2 + LDQ2*NQ
      IRECVQ = ISENDQ + LDQ2*NQ
      IBUF   = IRECVQ + LDQ2*NQ
!
      INDROW = 1
      INDCOL = INDROW + N
      INDX   = INDCOL + N
      INDRCV = INDX + N
!
      call eigen_get_procs(eigen_NP, eigen_NPROW, eigen_NPCOL)

      if(nnod.eq.eigen_NP)then
         CALL FS_PDLASRT( N, D, Q, LDQ, ROOT_NODE, WORK(IPQ2), LDQ2, &
              WORK(ISENDQ), WORK(IRECVQ), WORK(IBUF), &
              IWORK(INDROW), IWORK(INDCOL), IWORK(INDX), &
              IWORK(INDRCV), INFO, prof )

      else
         CALL FS2eigen_PDLASRT( N, D, Q, LDQ, ROOT_NODE, &
              WORK,WORK,WORK(max(1,NP*NQ/2+1)), &
              IWORK, INFO, prof )

      endif
!
 10   CONTINUE
!
!     free tree
!
      if(FS_COMM_MEMBER)then
         CALL FS_dividing_free(ROOT_NODE)
      endif
!
#if TIMER_PRINT
      CALL FS_prof_end(prof, 20)
#endif
#ifdef _DEBUGLOG
      if( FS_MYRANK.eq.0 ) then
        write(*,'(a,i0)') "FS_PDLAED0 end. INFO=", INFO
      endif
#endif
      RETURN
!
!     End of FS_PDLAED0
!
end subroutine FS_PDLAED0

